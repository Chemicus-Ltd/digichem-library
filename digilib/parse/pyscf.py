import hashlib
import json
from uuid import uuid4

from cclib.bridge.cclib2pyscf import cclibfrommethods

from digilib.parse.base import Parser_abc
from digilib.parse.std2 import Std2_parser_mixin
import digilib.file.types as file_types
import digilib.log


class Pyscf_parser(Parser_abc, Std2_parser_mixin):
    """
    Top level class for parsing output from pyscf data.
    """

    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
        file_types.xtb_std2: "std2_file"
    }
    
    def __init__(self, mol_name, methods, **kwargs):
        self.methods = methods
        self.mol_name = mol_name
        super().__init__(**kwargs)

    def _parse(self):
        """
        Extract results from our output files.
        """
        self.data = cclibfrommethods(**self.methods)

        # Do manual parsing of excited states.
        if 'std2_file' in self.auxiliary_files:
            self.parse_std2_log()
    
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """ 
        super().post_parse()
        # Set some metadata objects.
        self.data.metadata['name'] = self.mol_name
        self.data._aux = {'methods': self.methods}

        if hasattr(self.data, "etenergies"):
            # Reorder excited states.
            # First, set energies properly, keeping track of each energy's old index.
            energy_index = sorted(
                [(energy, index) for index, energy in enumerate(self.data.etenergies)],
                key=lambda energy_index: energy_index[0],
            )
            self.data.etenergies = [energy for energy, old_index in energy_index]

            # Sort everything else.
            for property in ("etdips", "etmagdips", "etoscs", "etrotats", "etsyms", "etsecs", "etveldips"):
                if hasattr(self.data, property):
                    setattr(self.data, property, [getattr(self.data, property)[old_index] for energy, old_index in energy_index])

        try:
            # Try to generate a checksum from metadata.
            self.data._id = hashlib.sha1(json.dumps(self.data.metadata, sort_keys = True, default = str).encode('utf-8')).hexdigest()

        except Exception:
            # No luck, something in metadata must be unhashable.
            digilib.log.get_logger().error("Unable to generate hash ID from calculation metadata, using random ID instead", exc_info = True)
            # TODO: Think of a better way to do this.
            self.data._id = hashlib.sha1(uuid4().hex.encode('utf-8')).hexdigest()
        
            