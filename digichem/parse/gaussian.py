# General imports.
from datetime import datetime, timedelta

from digichem.exception.base import Result_unavailable_error

# Digichem imports.
from digichem.parse.cclib import Cclib_parser
import digichem.log
import digichem.file.types as file_types

# Hidden imports.
#import pysoc.io.SOC

class Gaussian_parser(Cclib_parser):
    """
    Top level class for parsing output from Gaussian log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
            file_types.gaussian_chk_file: "chk_file",
            file_types.gaussian_fchk_file: "fchk_file",
            file_types.gaussian_rwf_file: "rwf_file"
        }
    
    # Headers for date strings.
    DATE_HEADER = "Normal termination of"
    ELAPSED_TIME_HEADER = "Elapsed time:"
    CPU_TIME_HEADER = "Job cpu time:"
    CPU_HEADER = "Will use up to"

    def __init__(self, *log_files, rwfdump = "rwfdump", options, ornt = None, ornt_args = (), metadata_defaults = None, **auxiliary_files):
        self.rwfdump = rwfdump
        super().__init__(*log_files, options = options, ornt = ornt, ornt_args = ornt_args, metadata_defaults = metadata_defaults, **auxiliary_files)
    
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        super().parse_metadata()
    
    def pre_parse(self):
        """
        Perform any setup before line-by-line parsing.
        """
        super().pre_parse()
        # Assume we used 1 CPU if not otherwise clear (is this a good idea?)
        self.data.metadata['num_cpus'] = 1
        
        self.wall_time = []
        self.cpu_time = []
    
    def parse_output_line(self, log_file, line):
        """
        Perform custom line-by-line parsing of an output file.
        """
        # Although we only need the last ~5 lines from the (possibly huge) log file, we read all the way through because negative seek()ing is tricky.
        # Look for our key string.
        if self.DATE_HEADER in line:
            # This line looks like: "Normal termination of Gaussian 16 at Sun Dec  6 19:13:09 2020"
            date_str = " ".join(line.split()[-4:])
            self.data.metadata['date'] = datetime.strptime(date_str, "%b %d %H:%M:%S %Y.").timestamp()
            
        elif self.ELAPSED_TIME_HEADER in line:
            # This line looks like: "Elapsed time:       0 days  2 hours 38 minutes 50.9 seconds."
            datey = line.split()[-8:]
            self.wall_time.append(timedelta(days = int(datey[0]), hours = int(datey[2]), minutes = int(datey[4]), seconds = float(datey[6])).total_seconds())
            
        elif self.CPU_TIME_HEADER in line:
            # This line looks like: "Job cpu time:       0 days 20 hours 52 minutes 17.3 seconds."
            datey = line.split()[-8:]
            self.cpu_time.append(timedelta(days = int(datey[0]), hours = int(datey[2]), minutes = int(datey[4]), seconds = float(datey[6])).total_seconds())
            
        elif self.CPU_HEADER in line:
            # This line looks like: "Will use up to   10 processors via shared memory."
            self.data.metadata['num_cpus'] = int(line.split()[4])
            
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """
        super().post_parse()
        
        if 'wall_time' not in self.data.metadata and len(self.wall_time) != 0:
            self.data.metadata['wall_time'] = self.wall_time
        
        if 'cpu_time' not in self.data.metadata and len(self.cpu_time) != 0:
            self.data.metadata['cpu_time'] = self.cpu_time
        
        # Get SOC.
        # Next try and get SOC.
        try:
            self.calculate_SOC()
            
        except Exception:
            digichem.log.get_logger().debug("Cannot calculate spin-orbit-coupling from output file '{}'".format(self.log_file_path), exc_info = True)
    
    def calculate_SOC(self):
        """
        Parse spin-orbit coupling using PySOC.
        """
        try:
            import pysoc.io.SOC
        
        except Exception as e:
            raise Result_unavailable_error("Spin-orbit coupling", "PySOC is not available") from e
        
        # For SOC, we need both .log and .rwf file.
        # No need to check for these tho; pysoc does that for us.
        # We also need etsyms to decide which excited state is which.
        if not hasattr(self.data, "etsyms"):
            raise Result_unavailable_error("Spin-orbit coupling", "There are no excited states available")
        
        # Get a PySOC parser.
        soc_calculator = pysoc.io.SOC.Calculator(self.log_file_path, rwfdump = self.rwfdump, rwf_file_name = self.auxiliary_files['rwf_file'])
        soc_calculator.calculate()
        SOC_table = soc_calculator.soc_td.SOC
        
        # We'll split the SOC table given to use by PySOC to better match the format used by cclib.
        socstates = []
        socelements = []
        
        for SOC_line in SOC_table:
            # Add states.
            socstates.append([SOC_line.singlet_state, SOC_line.triplet_state])
            
            # Add coupling.
            socelements.append([SOC_line.positive_one, SOC_line.zero, SOC_line.negative_one])
                
        # Add to data.
        self.data.socstates = socstates
        self.data.socelements = socelements
            
            