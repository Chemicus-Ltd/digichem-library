import hashlib
import json
import datetime
import periodictable
from uuid import uuid4
import types

from digichem.parse.base import Parser_abc
import digichem.log
from digichem.input.digichem_input import si_from_file



class Crest_parser(Parser_abc):
    """
    Top level class for parsing output from crest data.
    """
    
    def __init__(self, mol_name, program, **kwargs):
        self.program = program
        self.mol_name = mol_name
        super().__init__(**kwargs)

    def _parse(self):
        """
        Extract results from our output files.
        """
        self.data = types.SimpleNamespace()

        # Metadata we can get entirely from our passed in program object.
        self.data.metadata = {
            "name": self.mol_name,
            "user": self.program.calculation.user.username,
            "userID": self.program.calculation.user.userID,
            "jobId": self.program.calculation.job_id,
            "wall_time": [self.program.duration],
            "cpu_time": [self.program.duration * self.program.calculation.performance['num_cpu']],
            "date": datetime.datetime.now(datetime.timezone.utc).timestamp(),
            "package": "CREST",
            #"package_version": None,
            "calculations": ["Optimisation"],
            "success": not self.program.error,
            "methods": [self.program.calculation.method['gfn']['level']],
            "charge": self.program.calculation.charge,
            "multiplicity": self.program.calculation.multiplicity,
            "optimisation_converged": not self.program.error,
            #"temperature": None,
            #"pressure": None,
            # TODO: CREST doesn't actually use orbitals...
            "orbital_spin_type": "restricted",
            "solvent_name": str(self.program.calculation.solution['solvent']) if self.program.calculation.solution['calc'] else None,
            "solvent_model": self.program.calculation.solution['model'] if self.program.calculation.solution['calc'] else None,
            "num_cpu": self.program.calculation.performance['num_cpu'],
            #"memory_used": None,
            "memory_available": self.program.calculation.performance['memory'],
        }

        # Crest calculations are normally used for conformer searching, so we should expect multiple structures as output.
        # We'll use the lowest energy conformer as our 'main' structure.
        main_si = si_from_file(self.program.working_directory / "crest_conformers.xyz", gen3D = False, charge = self.program.calculation.charge, multiplicity = self.program.calculation.multiplicity)
        
        self.data.atomnos = []
        self.data.atomcoords = [[]]
        
        for atom in main_si.atoms:
            self.data.atomnos.append(periodictable.elements.symbol(atom['atom']).number)
            self.data.atomcoords[0].append(
                [
                    atom['x'],
                    atom['y'],
                    atom['z']
                ]
            )
    
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """ 
        super().post_parse()

        try:
            # Try to generate a checksum from metadata.
            self.data._id = hashlib.sha1(json.dumps(self.data.metadata, sort_keys = True, default = str).encode('utf-8')).hexdigest()

        except Exception:
            # No luck, something in metadata must be unhashable.
            digichem.log.get_logger().error("Unable to generate hash ID from calculation metadata, using random ID instead", exc_info = True)
            # TODO: Think of a better way to do this.
            self.data._id = hashlib.sha1(uuid4().hex.encode('utf-8')).hexdigest()
        
            