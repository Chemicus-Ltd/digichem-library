import numpy as np

from digichem.result import Result_object


class Std2_parser_mixin():
    """Mixin for classes that can calculate excited states with std2"""

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

            # Some orbitals are missing in the sTD-DFT treatment (probably core e- ?).
            orbital_offset = 0

            for line in std2_file:
                if "triplet" in line:
                    " triplet                       :  T"
                    if line.split()[-1] == "T":
                        symmetry = "Triplet"
                    else:
                        symmetry = "Singlet"
                
                if " ordered frontier orbitals" in line:
                    #  ordered frontier orbitals
                    #          eV     # centers
                    #   84    -7.384     3.5
                    #   85    -7.362     3.4
                    #   86    -7.215    12.1
                    #   87    -7.188     8.4
                    #   88    -7.142    13.9
                    #   89    -6.987    13.9
                    #   90    -6.962    12.2
                    #   91    -6.877    10.6
                    #   92    -6.722     9.0
                    #   93    -6.604     9.4
                    #   94    -5.108     8.0

                    #   95    -2.189    10.7
                    #   96    -2.129    10.2
                    #   97    -0.858    16.7
                    #   98    -0.680     6.6
                    #   99    -0.366     8.6
                    #  100    -0.364     8.9
                    #  101    -0.286     9.8
                    #  102    -0.199    10.2
                    #  103     0.135    11.1
                    #  104     0.634    13.0
                    #  105     0.784    11.6
                    line = next(std2_file)

                    while (line := next(std2_file)).strip():
                        last_orbital = int(line.split()[0]) -1
                    
                    # Adjust based on where the HOMO really is.
                    orbital_offset = self.data.homos[0] - last_orbital


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
                                (startMO + orbital_offset, 0),
                                (endMO + orbital_offset, 0),
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