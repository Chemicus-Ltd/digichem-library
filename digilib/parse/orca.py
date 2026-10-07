from digilib.parse.cclib import Cclib_parser
import digilib.file.types as file_types
import digilib.log


class Orca_parser(Cclib_parser):
    """
    Top level class for parsing output from Orca log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
            file_types.orca_gbw_file: "gbw_file",
            file_types.orca_density_file: "density_file",
            file_types.orca_density_info_file: "density_info_file"
        }
    
    def parse_output_line(self, log_file, line):
        """
        Perform custom line-by-line parsing of an output file.
        """
        
        # Spin-orbit coupling.
        # We're looking to populate two or three attributes
        # 'socstates'   : A two-membered list of the singlet and triplet symbols that make up this coupling (eg, ["S(1)", "T(2)"].
        # 'socenergies' : The total spin-orbit coupling value (RSS of socelements).
        # 'socelements' : A three-membered list of the soc values for the triplet state with number +1, 0 and -1.
        if "CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS" in line:
            # Start of the SOC section.
            # The same header is used for SOC in cartesian basis (x,y,z) and for individual triplet states (+1, 0, -1).
            line = next(log_file)
            line = next(log_file)
            line = next(log_file)
            
            soc_type = None
            if "Z" in line and "X" in line and "Y in line":
                # Cartesian SOC
                # In this format we can parse total SOC only.
                soc_type = "cartesian"
            elif "0" in line and "-1" in line and "+1" in line:
                # Triplet SOC.
                # In this format we can parse individual SOC as well as total SOC.
                soc_type = "triplet"
            
            else:
                pass
                digilib.log.get_logger().debug("Unrecognised SOC section started by line '{}'".format(line))
            
            if soc_type is not None:
                # Reset our attributes.
                self.data.socstates = []
                self.data.socenergies = []
                if soc_type == "triplet":
                    self.data.socelements = []
                elif hasattr(self.data, 'socelements'):
                    delattr(self.data, 'socelements')
                
                line = next(log_file)
                line = next(log_file)
                
                while line.strip() != "--------------------------------------------------------------------------------" and \
                    line.strip() != "":
                    # Each line is the coupling between one singlet state and one triplet state.
                    # 1      1    (   0.00 ,    0.00)    (  -0.00 ,   -0.00)    (  -0.00 ,    0.00)
                    split_line = line.split()
                    try:
                        triplet_index = int(split_line[0])
                        singlet_index = int(split_line[1])
                    
                    except Exception:
                        print(line)
                    
                    # Split on brackets to get each xyz/0, -1, +1 element.
                    soc_elements = []
                    soc_element_strings = line.split("(")[1:]
                    
                    for soc_element_string in soc_element_strings:
                        # The last character will be the closing bracket, so we can discard.
                        real, imagine = [float(ele) for ele in soc_element_string.strip()[:-1].split(",")]
                        
                        # We're not interested in the real or imaginary parts, just combine (root of the sum of the squares).
                        soc_elements.append((real**2 + imagine**2)**0.5)
                    
                    # We now have everything we need.
                    self.data.socstates.append(["S({})".format(singlet_index), "T({})".format(triplet_index)])
                    self.data.socenergies.append((soc_elements[0]**2 + soc_elements[1]**2 + soc_elements[2]**2)**0.5)
                    
                    # Only add elements if they are triplet elements (not cartesian).
                    if soc_type == "triplet":
                        # The order of triplet states is different in Orca to what we expect.
                        # We want +1, 0, -1
                        # Orca is 0, -1, +1
                        self.data.socelements.append([soc_elements[2], soc_elements[0], soc_elements[1]])
                
                    line = next(log_file)
                
                