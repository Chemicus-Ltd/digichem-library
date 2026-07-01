import numpy as np

from digichem.parse.cclib import Cclib_parser
import digichem.file.types as file_types
import digichem.log


class Xtb_parser(Cclib_parser):
    """
    Top level class for parsing output from Xtb log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
        file_types.xtb_topology: "topology_file",
        file_types.xtb_molden: "molden_file"
    }
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        # Avoid circular import.
        from digichem.input import si_from_file

        # CClib does most of the work for us.
        super()._parse()

        # If atomcoords are missing, add those from topo file.
        if self.data.atomcoords.size == 0 and self.auxiliary_files['topology_file'] is not None:
            # Parse the file.
            si = si_from_file(self.auxiliary_files['topology_file'])

            # Remember that atomcoords has a dimension for opt cycles.
            self.data.atomcoords = np.array([[[atom['x'], atom['y'], atom['z']] for atom in si.atoms]])
