import numpy as np

from digichem.parse.cclib import Cclib_parser
import digichem.file.types as file_types
import digichem.log
from digichem.result import Result_object


class Xtb_parser(Cclib_parser):
    """
    Top level class for parsing output from Xtb log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
        file_types.xtb_topology: "topology_file",
        file_types.xtb_molden: "molden_file",
        file_types.xtb_density: "density_file",
        file_types.xtb_std2: "std2_file"
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
        
        # Do manual parsing of excited states.
        if 'std2_file' in self.auxiliary_files:
            self.parse_std2_log()
    
    def parse_std2_log(self):
        """
        Parse excited states from an std2.log file.
        """
        with open(self.auxiliary_files['std2_file'], "rt") as std2_file:
            # Setup variables.
            symmetry =  ""
            etenergies = []
            etoscs = []
            etsyms = []
            etsecs = []

            for line in std2_file:
                if "triplet" in line:
                    " triplet                       :  T"
                    if line.split()[-1] == "T":
                        symmetry = "Triplet"
                    else:
                        symmetry = "Singlet"

                if " excitation energies, transition moments and" in line:
                    #  excitation energies, transition moments and TDA amplitudes
                    #  state    eV      nm       fL        Rv(corr)
                    #     1    4.979   249.0     0.0000    -0.0000     0.74(  14->  16) -0.67(  15->  17)  0.04(  15->  16)
                    #     2    5.002   247.9     0.0000    -0.0000     1.00(  15->  16) -0.05(  11->  26)  0.03(  15->  17)
                    #     3    5.010   247.5     0.0000     0.0000    -0.74(  15->  17) -0.67(  14->  16)  0.01(  15->  21)
                    # ...
                    line = next(std2_file)
                    while len(split_line := next(std2_file).split()) != 0:
                        index = int(split_line[0])
                        energy = float(split_line[1])
                        wavelength = float(split_line[2])
                        oscillator = float(split_line[3])

                        raw_configs = [split_line[5 + pos:5 + pos + 3] for pos in range(0, len(split_line[5:]), 3)]
                        configurations = []

                        for raw_config in raw_configs:
                            startMO = int(raw_config[1][:-2]) -1
                            endMO = int(raw_config[2][:-1]) -1
                            coef = float(raw_config[0][:-1])

                            configurations.append([
                                (startMO, 0),
                                (endMO, 0),
                                coef
                            ])
                        
                        etenergies.append(Result_object.energy_to_wavenumbers(energy))
                        etoscs.append(oscillator)
                        etsyms.append(symmetry)
                        etsecs.append(configurations)
            
            # Reorder.
            # First, set energies properly, keeping track of each energy's old index.
            energy_index = sorted(
                [(energy, index) for index, energy in enumerate(etenergies)],
                key=lambda energy_index: energy_index[0],
            )

            # Sort everything else.
            etenergies = [energy for energy, old_index in energy_index]
            etoscs = [etoscs[old_index] for energy, old_index in energy_index]
            etsyms = [etsyms[old_index] for energy, old_index in energy_index]
            etsecs = [etsecs[old_index] for energy, old_index in energy_index]

            # Save.        
            self.data.etenergies = np.array(etenergies)
            self.data.etoscs = np.array(etoscs)
            self.data.etsyms = etsyms
            self.data.etsecs = etsecs