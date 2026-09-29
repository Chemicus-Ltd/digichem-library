from pathlib import Path
import itertools

import digichem.log
from digichem.parse.base import File_parser_abc
from digichem.exception import Digichem_exception
from digichem.misc.io import checksum

# Hidden imports.
#import cclib.io


class Cclib_parser(File_parser_abc):
    """
    ABC for parsers that use cclib to do most of their work for them.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {}
    
    def __init__(self, *log_files, options, ornt = None, ornt_args = (), metadata_defaults = None, **auxiliary_files):
        """
        Top level constructor for calculation parsers.
        
        :param log_files: A list of output file to analyse/parse. The first log_file given will be used for naming purposes.
        :param auxiliary_files: A dictionary of auxiliary input files related to the calculation.
        """
        # Also save our aux files, stripping None.
        self.auxiliary_files = {name: aux_file for name,aux_file in auxiliary_files.items() if aux_file is not None}
        
        super().__init__(*log_files, options = options, ornt = ornt, ornt_args = ornt_args, metadata_defaults = metadata_defaults, profile_file = Path(log_files[0].parent, "../Logs/profile.csv"))
        
    @classmethod
    def from_logs(self, *log_files, hints = None, options, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file(s).
        
        :param given_log_files: Output file(s) to parse or a directory of output files to parse.
        """
        # Have a look for aux. files.
        auxiliary_files = {}

        basename = log_files[0].name if len(log_files) > 0 else ""
        
        for hint in itertools.chain(log_files, hints if hints is not None else []):
            auxiliary_files.update(self.find_auxiliary_files(hint, basename))
            
        # Finally, update our auxiliary_files with kwargs, so any user specified aux files take precedence.
        auxiliary_files.update(kwargs)
        
        return self(*log_files, options = options, **auxiliary_files)
    
    @classmethod
    def find_auxiliary_files(self, hint, basename):
        """
        Find auxiliary files from a given hint.
        
        :param hint: A path to a file to use as a hint to find additional files.
        :returns: A dictionary of found aux files.
        """
        hint = Path(hint)
        
        # Now have a look for aux. input files, which are defined by each parser's INPUT_FILE_TYPES
        auxiliary_files = {}
        for file_type in self.INPUT_FILE_TYPES:
            for extension in file_type.extensions:
                if hint.is_dir():
                    # Peak inside.
                    if Path(hint, basename).with_suffix(extension).exists():
                        auxiliary_files[self.INPUT_FILE_TYPES[file_type]] = Path(hint, basename).with_suffix(extension)

                if hint.with_suffix(extension).exists():
                    auxiliary_files[self.INPUT_FILE_TYPES[file_type]] = hint.with_suffix(extension)
        
        return auxiliary_files
    
    def _parse(self):
        """
        Extract results from our output files.
        
        This default implementation will parse the log file with cclib and save the result to self.data.
        """
        import cclib.io
        
        # We start by using cclib to get most of the data we need.
        
        # Output a message (because this is slow).
        digichem.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        
        # Use cclib to open our log files.
        # ccread will accept a list of log files to read, but will sometimes choke if the list contains only one entry,
        # in which case we give it only one file name.
        # ccread will also choke if we give it pathlib Paths.
        file_paths = [str(log_file) for log_file in self.log_file_paths]
        
        # Get data from cclib.
        self.data = cclib.io.ccread(file_paths if len(file_paths) > 1 else file_paths[0])
        if self.data is None:
            raise Digichem_exception("Could not parse any data at all!")
        
        # Get a unique ID (checksum) from the given log files.
        # First, order the list of filenames so we also process in the same order.
        # We do this because not all parsers define a custom sort.
        
        file_paths.sort()
        self.data._id = checksum(*file_paths, hash_func = "sha1")
        
        # Do some setup.
        self.pre_parse()
        
        # Do our own parsing (if any).
        for log_file_path in self.log_file_paths:
            with open(log_file_path, "rt") as log_file:
                for line in log_file:
                    self.parse_output_line(log_file, line)

        # Add profiling data.
        try:
            self.parse_profile_file()
        
        except Exception:
            if self.profile_file.exists():
                digichem.log.get_logger().warning("Could not parse profile.csv file; profiling data will be unavailable", exc_info=True)
            
            else:
                pass
        
    def parse_output_line(self, log_file, line):
        """
        Perform custom line-by-line parsing of an output file.
        
        This method will be called for each line of each log-file given to the parser (although be aware that some implementations may skip some lines during parsing),
        and it allows for data not supported by cclib to be extracted. It is program specific.
        """
        # Do nothing.
        
    def pre_parse(self):
        """
        Perform any setup before line-by-line parsing.
        """
        # Do nothing.
        
    @classmethod
    def au_to_debye(self, au):
        """
        Convert a dipole moment in au to debye.
        """
        return au * 2.541746473
    
    @classmethod
    def bohr_to_angstrom(self, bohr_distance):
        """
        Convert a length in bohr to angstrom.
        """
        return bohr_distance * 0.529177210903