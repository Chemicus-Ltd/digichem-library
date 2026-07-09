import numpy as np
from datetime import timedelta

import digichem.log
from digichem.parse.cclib import Cclib_parser
from digichem.parse.std2 import Std2_parser_mixin
import digichem.file.types as file_types


class Xtb_parser(Cclib_parser, Std2_parser_mixin):
    """
    Top level class for parsing output from Xtb log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
        file_types.xtb_topology: "topology_file",
        file_types.xtb_molden: "molden_file",
        file_types.xtb_density: "density_file",
        file_types.xtb_std2: "std2_file",
        file_types.xtb_xtb4stda: "xtb4stda_file"
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
        
        # Do manual parsing of xtb4stda.
        if 'xtb4stda_file' in self.auxiliary_files:
            try:
                self.parse_xtb4stda_log()
            
            except Exception:
                digichem.log.get_logger().warning("Failed to parse xtb4stda aux file", exc_info = True)
        
        # Do manual parsing of excited states.
        if 'std2_file' in self.auxiliary_files:
            self.parse_std2_log()
    
    def parse_xtb4stda_log(self):
        """
        Parse data from the xtb4stda log file.
        """
        # Timing.
        # This section has variable output sadly...
        # speedup  2.72
        # cpu  time for all    8.28 s
        # wall time for all    3.05 s
        with open(self.auxiliary_files['xtb4stda_file'], "rt") as xtb4stda_file:
            for line in xtb4stda_file:
                if "cpu  time for all" in line or "wall time for all" in line:
                    split_line = line.split()
                    days = 0
                    hours = 0
                    minutes = 0
                    seconds = 0
                    if len(split_line) > 10:
                        days = float(split_line[-8])
                    
                    if len(split_line) > 8:
                        hours = float(split_line[-6])
                    
                    if len(split_line) > 6:
                        minutes = float(split_line[-4])
                    
                    seconds=float(split_line[-2])
                    
                    if "cpu" in line:
                        if 'cpu_time' not in self.data.metadata:
                            self.data.metadata['cpu_time'] = []
                            
                        self.data.metadata['cpu_time'].append(
                            timedelta(
                                days = days,
                                hours = hours,
                                minutes = minutes,
                                seconds = seconds
                            )
                        )
                    else:
                        if 'wall_time' not in self.data.metadata:
                            self.data.metadata['wall_time'] = []

                        self.data.metadata['wall_time'].append(
                            timedelta(
                                days = days,
                                hours = hours,
                                minutes = minutes,
                                seconds = seconds
                            )
                        )