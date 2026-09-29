import hashlib
import json
import datetime
import periodictable
from uuid import uuid4
import types

from digilib.parse.base import Parser_abc
import digilib.log
from digilib.input.digichem_input import si_from_file


class Censo_parser(Parser_abc):
    """
    Top level class for parsing output from censo data.
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

        # Censo calculations are normally made up of four distinct phases (prescreening, screening, optimisation, and refinement).
        calculations = []
        methods = []
        functional = []
        basis_set = []
        engines = []

        if self.program.calculation.properties['prescreening']['calc']:
            calculations.append("Screening")
            methods.append(self.program.calculation.properties['prescreening']['gfn'])
            methods.append("DFT")
            engines.append(self.program.calculation.properties['prescreening']['engine'])
            functional.append(self.program.calculation.properties['prescreening']['functional'])
            basis_set.append(self.program.calculation.properties['prescreening']['basis_set'])
        
        if self.program.calculation.properties['screening']['calc']:
            calculations.append("Screening")
            methods.append(self.program.calculation.properties['screening']['gfn'])
            methods.append("DFT")
            engines.append(self.program.calculation.properties['screening']['engine'])
            functional.append(self.program.calculation.properties['screening']['functional'])
            basis_set.append(self.program.calculation.properties['screening']['basis_set'])
        
        if self.program.calculation.properties['optimisation']['calc']:
            calculations.append("Optimisation")
            methods.append(self.program.calculation.properties['optimisation']['gfn'])
            methods.append("DFT")
            engines.append(self.program.calculation.properties['optimisation']['engine'])
            functional.append(self.program.calculation.properties['optimisation']['functional'])
            basis_set.append(self.program.calculation.properties['optimisation']['basis_set'])
        
        if self.program.calculation.properties['refinement']['calc']:
            calculations.append("Single Point")
            methods.append(self.program.calculation.properties['refinement']['gfn'])
            methods.append("DFT")
            engines.append(self.program.calculation.properties['refinement']['engine'])
            functional.append(self.program.calculation.properties['refinement']['functional'])
            basis_set.append(self.program.calculation.properties['refinement']['basis_set'])

        engines = ["CENSO"] + list(set(engines))

        # Metadata we can get entirely from our passed in program object.
        self.data.metadata = {
            "name": self.mol_name,
            "user": self.program.calculation.user.username,
            "userID": self.program.calculation.user.userID,
            "jobId": self.program.calculation.job_id,
            "wall_time": [self.program.duration],
            "cpu_time": [self.program.duration * self.program.calculation.performance['num_cpu']],
            "date": datetime.datetime.now(datetime.timezone.utc).timestamp(),
            "package": "/".join(engines),
            #"package_version": None,
            "calculations": calculations,
            "success": not self.program.error,
            "methods": methods,
            "functional": "/".join(functional),
            "basis_set": "/".join(basis_set),
            "charge": self.program.calculation.charge,
            "multiplicity": self.program.calculation.multiplicity,
            "optimisation_converged": not self.program.error,
            #"temperature": None,
            #"pressure": None,
            "orbital_spin_type": "restricted" if self.program.calculation.multiplicity == 1 else "unrestricted",
            "solvent_name": str(self.program.calculation.solution['solvent']) if self.program.calculation.solution['calc'] else None,
            "solvent_model": self.program.calculation.solution['model'] if self.program.calculation.solution['calc'] else None,
            "num_cpu": self.program.calculation.performance['num_cpu'],
            #"memory_used": None,
            "memory_available": self.program.calculation.performance['memory'],
        }

        try:
            coord_file = self.program.next_coords
        
        except Exception:
            coord_file = None
        
        if coord_file:
            main_si = si_from_file(coord_file, gen3D = False, charge = self.program.calculation.charge, multiplicity = self.program.calculation.multiplicity)
            
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
            digilib.log.get_logger().error("Unable to generate hash ID from calculation metadata, using random ID instead", exc_info = True)
            # TODO: Think of a better way to do this.
            self.data._id = hashlib.sha1(uuid4().hex.encode('utf-8')).hexdigest()
        
            