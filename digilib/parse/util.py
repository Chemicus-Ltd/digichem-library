"""Utilities and handy functions for reading and parsing result files."""

# General imports.
from functools import partial
from itertools import filterfalse, zip_longest
from pathlib import Path
import itertools
import shutil
from tempfile import mkdtemp
import collections
import warnings
# IMPORTANT: Do not replace multiprocessing pools with pathos, the latter is too buggy for production ATM (26-05-2023).
import multiprocessing

from configurables.misc import is_iter, is_int
from configurables.defres import Default, defres

# Digichem imports.
from digilib.exception.base import Digichem_exception
from digilib.parse.cclib import Cclib_parser
from digilib.parse.gaussian import Gaussian_parser
from digilib.parse.turbomole import Turbomole_parser
from digilib.parse.orca import Orca_parser
from digilib.parse.xtb import Xtb_parser
from digilib.result.multi import Merged
from digilib.result.result import Result_set
import digilib.log
from digilib.parse.dump import Yaml_multi_parser, Json_multi_parser

# Hidden imports.
#import cclib.io
#import cclib.parser

custom_parsing_formats = {
    "sir": Yaml_multi_parser,
    "yaml": Yaml_multi_parser,
    "sij": Json_multi_parser,
    "json": Json_multi_parser,
}

archive_formats = list(itertools.chain(*[extensions for name, extensions, desc in shutil.get_unpack_formats()]))

def find_log_files_from_hint(hint):
    """
    Find output (log) files from a given hint.
    
    :param hint: A path to a file to use as a hint to find additional log files. hint can optionally be a directory, in which case files inside this directory will be found.
    :returns: A list of found log files.
    """
    hint = Path(hint)
    
    # First, find our parent dir.
    # hint may actually be a dir.
    if hint.is_dir():
        # Look for all .log files.
        # File extensions that we recognise.
        log_types = itertools.chain(["*." + custom_format for custom_format in custom_parsing_formats], ["*.log", ]) # "*.out" disabled for now...
        parent = hint
        log_files = [found_log_file for found_log_file in itertools.chain(*[parent.glob(log_type) for log_type in log_types])]
        
        # Remove any 'digichem.log' files as we know these are not calc log files.
        # We don't actually write 'digichem.log' files anymore either (we use digichem.out instead),
        # but older versions did...
        log_files = [log_file for log_file in log_files if log_file.name not in ["digichem.log", "digichem.out", "silico.log", "silico.out", ".digichem.yaml"]]
    else:
        parent = hint.parent
        log_files = [hint]
    
    # If we have a computational style log file, look for others.
    if hint.suffix not in ["." + custom_format for custom_format in custom_parsing_formats]:
        # Try and find job files.
        # These files have names like 'job.0', 'job.1' etc, ending in 'job.last'.
        for maybe_job_file in parent.glob("job.*"):
            # We only want the numbered files.
            if is_int(".".join(maybe_job_file.name.split(".")[1:])):
                # Looks good.
                log_files.append(maybe_job_file)
        
        # Look for other files.
        for maybe_file_name in ("basis", "control", "mos", "alpha", "beta", "coord", "gradient", "aoforce", "job.last", "numforce/aoforce.out"):
            maybe_file_path = Path(parent, maybe_file_name)
            
            if maybe_file_path.exists():
                # Found it.
                log_files.append(maybe_file_path)
            
    # Make sure we only have unique log files.
    # We also now reverse our ordering, so that files given earlier by the user have priority.
    if len(log_files) > 0:
        return log_files
            
    # If we have no log files, and there's a directory called Output or Result that we can use, try again using that as the hint.
    elif hint.is_dir():
        # Try the Output dir first.
        if Path(hint, "Output").is_dir():
            log_files = find_log_files_from_hint(Path(hint, "Output"))
            
        # If we still have nothing, try the Results dir.
        if len(log_files) == 0 and Path(hint, "Results").is_dir():
            log_files = find_log_files_from_hint(Path(hint, "Results"))
            
    return log_files


def find_log_files(*hints):
    """
    Find log files from a number of given hints.
    
    Each hint should point to an existing log file to parse, or a directory containing such log files.
    Each log file given (and found) should refer to the same calculation. 
    """
    # Get a list of found log files.
    log_files = [found_log for hint in hints for found_log in find_log_files_from_hint(hint)]
    
    # Make sure we only have unique log files.
    # We also now reverse our ordering, so that files given earlier by the user have priority.
    log_files = list(reversed(list(dict.fromkeys([log_file.resolve() for log_file in log_files]))))

    return log_files
    

def class_from_log_files(*log_files, format_hint = "auto"):
    """
    Get a parser class based on some calculation log files.
    
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'cclib' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    """
    import cclib.io
    import cclib.parser
    
    if format_hint in custom_parsing_formats:
        return custom_parsing_formats[format_hint]
    
    elif format_hint == "auto" and len(log_files) > 0 and log_files[0].suffix[1:] in custom_parsing_formats:
        return custom_parsing_formats[log_files[0].suffix[1:]]
    
    if format_hint not in ["cclib", "auto"]:
        raise ValueError("Unrecognised format hint '{}'".format(format_hint))
    
    # We'll use cclib to guess the file type for us.
    try:
        log_file_type = type(cclib.io.ccopen([str(found_log_file) for found_log_file in log_files]))
        
    except Exception as e:
        if isinstance(e, FileNotFoundError):
            raise

        # cclib couldn't figure out the file type, it probably wasn't a .log file.
        raise Digichem_exception("Could not determine file type of file(s): '{}'; are you sure these are computational log files?".format(", ".join((str(log_file) for log_file in log_files)))) from e
    
    # Either get a more complex parser if we have one, or just return the base parser.
    if log_file_type == cclib.parser.gaussianparser.Gaussian:
        return Gaussian_parser
    elif log_file_type == cclib.parser.turbomoleparser.Turbomole:
        return Turbomole_parser
    elif log_file_type == cclib.parser.orcaparser.ORCA:
        return Orca_parser
    elif log_file_type == cclib.parser.xtbparser.XTB:
        return Xtb_parser
    else:
        return Cclib_parser

def from_log_files(*log_files, format_hint = "auto", parser_options = {}, **auxiliary_files):
    """
    Get an output file parser of appropriate type.
    
    This is a convenience function.
    
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    """
    found_log_files = find_log_files(*log_files)

    try:
        return class_from_log_files(*found_log_files, format_hint = format_hint).from_logs(*found_log_files, hints = log_files, **parser_options, **auxiliary_files)
     
    except Exception:
        if len(found_log_files) == 0:
            raise ValueError("There are no log files at '{}'".format(log_files[0] if len(log_files) == 1 else log_files)) from None
         
        else:
            raise
    
def parse_calculation(*log_files, options, ornt = None, ornt_args = (), parse_all = False, format_hint = "auto", keep_archive = False, parser_options = {}, **auxiliary_files):
    """
    Parse a single calculation result.
    
    :param log_files: A number of calculation result files corresponding to the same calculation.
    :param options: A Digichem options nested dictionary containing options to control parsing.
    :param parse_all: Whether to parse all results in the given log file. If True, a list of result sets will be returned, if False, only the first result will be returned if there are multiple.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    :param auxiliary_files: Optional auxiliary calculation files corresponding to the calculation.
    :return: A single Result_set object.
    """        
    # Handle aux files.
    # Auxiliary files are files associated with a calculation but that do not contain calculation output directly (they are not log files).
    # Often, these files are written in a program-dependent binary format, and may be used for eg, program restarting, post-calculations, or image generation.
    # For example, common aux files for Gaussian include: .chk, .fchk and .rwf.
    # Auxiliary files to associate with a log file(s) can be given by the auxiliary_files key-word argument, but this is cumbersome in some instances.
    # Alternatively, aux files can be given as normal log_files by prefixing the file name with 'aux':, where aux is the type of file.
    # For example, a chk file could be given as "chk:calc/output.chk".
    real_log_files = []
    for maybe_log_file in log_files:
        maybe_log_file = str(maybe_log_file)
        found = False
        # Loop through all known auxiliary_files:
        for file_type, aux_file_name in itertools.chain(Gaussian_parser.INPUT_FILE_TYPES.items(), Turbomole_parser.INPUT_FILE_TYPES.items(), Orca_parser.INPUT_FILE_TYPES.items()):
            code_len = len(file_type.short_code) +1
            if maybe_log_file[:code_len] == file_type.short_code + ":":
                auxiliary_files[aux_file_name] = maybe_log_file[code_len:]
            
                # Done.
                found = True
                break
        
        if not found:
            # Take care of 'log:' files (which we support in-case a normal log file happened to start with 'chk:' or similar).
            if maybe_log_file[:4] == "log:":
                real_log_files.append(maybe_log_file[4:])
                
            else:
                # Take as is.
                real_log_files.append(maybe_log_file)
            
    log_files = real_log_files
        
    # Open files for reading (handles archives for us).
    archive = open_for_parsing(*log_files, auxiliary_files = auxiliary_files)
    
    try:
        open_log_files, open_aux_files = archive.open()
        
        if parse_all:
            results = from_log_files(*open_log_files, format_hint = format_hint, parser_options = parser_options, options = options, ornt = ornt, ornt_args = ornt_args, **open_aux_files).process_all()
        
        else:       
            results = from_log_files(*open_log_files, format_hint = format_hint, parser_options = parser_options, options = options, ornt = ornt, ornt_args = ornt_args, **open_aux_files).process()
        
    finally:
        if not keep_archive:
            archive.cleanup()
            
    if keep_archive:
        # We've been asked not to close the archive, return it.
        return (results, archive)
    
    else:
        # The caller isn't interested in the archive.
        return results

def multi_parser(log_files, auxiliary_files, *, options, ornt = None, ornt_args = (), format_hint = "auto", keep_archive = False, parser_options = {},):
        """
        The inner function which will be called in parallel to parse files.
        """
        # If the given 'log_file' is actually already a result object, then there's nothing for us to do.
        # This is allowed to support merging calculation results between previously parsed results and new log files.
        if isinstance(log_files, Result_set):
            # Nothing we need to do.
            return [log_files]
        
        # Next, decide if this is a single log file path, or is actually a list of multiple paths.
        # Regardless of whether this is a single file or multiple, all given files should correspond to the same calculation.
        if not isinstance(log_files, str) and is_iter(log_files):
            # Multiple paths.
            logs = log_files
        
        else:
            # Single path.
            logs = (log_files,)
        
        try:    
            return parse_calculation(*logs, options = options, ornt = ornt, ornt_args = ornt_args, parse_all = True, format_hint = format_hint, keep_archive = keep_archive, parser_options = parser_options, **auxiliary_files)
            
        except Exception:
            digilib.log.get_logger().warning("Unable to parse calculation result file '{}'; skipping".format(logs[0] if len(logs) == 1 else logs), exc_info = True)
            return None

def parse_multiple_calculations(*log_files, auxiliary_files = None, options, ornt = None, ornt_args = (), parser_options = {}, pool = None, init_func = None, init_args = None, format_hint = "auto", processes = 1, keep_archive = False):
    """
    Parse a number of separate calculation results in parallel.
    
    If the argument 'pool' is not given, a multiprocessing.pool object will first be created using the arguments init_func and num_cpu.
    
    :param log_files: A number of calculation result files corresponding to different calculations. Each item can optionally be a list itself, to specify files from the same calculation but which are spread across multiple files.
    :param auxiliary_files: A list of dicts of aux files. The ordering of the dicts should correspond to that of log_files.
    :param options: A Digichem options nested dictionary containing options to control parsing.
    :param pool: An optional subprocessing.pool object to use for parallel parsing.
    :param init_func: An optional function to call to init each newly created process.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    :param processes: The max number of processes to create the new pool object with; if the number of given log_files is less than processes, then len(log_files) will be used instead.
    """
    if len(log_files) == 0:
        # Give up now.
        return []
    
    if auxiliary_files is None:
        auxiliary_files = [{}] * len(log_files)
    
    if len(log_files) < processes:
        processes = len(log_files)
    
    # Sort out our pool if we need to.
    own_pool = False
    if pool is None:
        own_pool = True
        pool = multiprocessing.Pool(processes, initializer = init_func, initargs = init_args if init_args is not None else [])
    
    # Do some parsing.
    try:
        result_lists = list(
            filterfalse(lambda x: x is None,
                pool.starmap(partial(multi_parser, options = options, ornt = ornt, ornt_args = ornt_args, format_hint = format_hint, keep_archive = keep_archive, parser_options = parser_options), zip_longest(log_files, auxiliary_files, fillvalue = {}))
                #pool.map(partial(multi_parser, options = options, format_hint = format_hint, keep_archive = keep_archive), *transpose(list(zip_longest(log_files, auxiliary_files, fillvalue = {})), 2))
            )
        )
        
        if keep_archive:
            return [(result, archive) for results, archive in result_lists for result in results]
        
        else:
            return [result for result_list in result_lists for result in result_list]
    
    finally:
        # Do some cleanup if we need to.
        if own_pool:
            pool.__exit__(None, None, None)
    
def parse_and_merge_calculations(*log_files, auxiliary_files = None, options, ornt = None, ornt_args = (), parser_options = {}, format_hint = "auto", inner_pool = None, keep_archive = False):
    """
    Get a single result object by parsing a number of computational log files.
    
    Multiple log files can be given both from the same calculation, or from multiple different calculations.
    If multiple different calculations are given, the individually parsed results will be merged together (which may give bizarre results if the calculations are unrelated, eg if they are of different molecules).
    
    Example:
        parse_and_merge_calculations(['calc1/primary.log', 'calc1/secondary.log'], 'calc2/calc.log', 'calc3/calc.log')
    Would parse three separate calculations (calc1, calc2 and calc3), of which the first is contained in two output files (primary.log and secondary.log), merging the result sets together.
    
    :param log_files: A list of paths to computational chemistry log files to parse. If more than one file is given, each is assumed to correspond to a separate calculation in which case the parsed results will be merged together. In addition, each given 'log file' can be an iterable of log file paths, which will be considered to correspond to an individual calculation.
    :param options: A Digichem options nested dictionary containing options to control parsing.: An alignment class to use to reorientate each molecule.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    :param auxiliary_files: A list of dictionaries of auxiliary files. The ordering of auxiliary_files should match that of log_files.
    :return: A single Result_set object (or child thereof).
    """
    parsed_results = parse_multiple_calculations(*log_files, options = options, ornt = ornt, ornt_args = ornt_args, parser_options = parser_options, format_hint = format_hint, pool = inner_pool, auxiliary_files = auxiliary_files, keep_archive = keep_archive)
    
    # If we asked for archives as well, unpack.
    if keep_archive:
        parsed_results, archives = list(map(list, zip(*parsed_results)))
    
    # If we have more than one result, merge them together.
    if len(parsed_results) > 1:
        parsed_results = Merged.from_results(*parsed_results, options = options, ornt = ornt, ornt_args = ornt_args)
    
    elif len(parsed_results) == 0:
        parsed_results = None
    
    else:
        parsed_results = parsed_results[0]
        
    if keep_archive:
        return (parsed_results, archives)
    
    else:
        return parsed_results
            
def multi_merger_parser(log_files, auxiliary_files, *, options, ornt = None, ornt_args = (), parser_options = {}, format_hint = "auto" , inner_pool = None, keep_archive = False):
        """
        The inner function which will be called in parallel to parse files.
        """
        try:
            return parse_and_merge_calculations(*log_files, options = options, ornt = ornt, ornt_args = ornt_args, parser_options = parser_options, format_hint = format_hint, inner_pool = inner_pool, auxiliary_files = auxiliary_files, keep_archive = keep_archive)
            
        except Exception:
            digilib.log.get_logger().warning("Unable to parse and merge calculation results '{}'; skipping".format(", ".join([str(log_file) for log_file in log_files])), exc_info = True)
            return None

def parse_and_merge_multiple_calculations(*multiple_results, options, ornt = None, ornt_args = (), parser_options = {}, format_hint = "auto", init_func = None, init_args = None, processes = None, auxiliary_files = None, keep_archive = False):
    """
    Parse a number of separate calculation results in parallel, merging some or all of the results into combined result sets.
    
    :param multiple_results: A list of two dimensions, where the first dimension is a list of separate results to process, and the second dimension is a list of results that should be merged together.
    :param options: A Digichem options nested dictionary containing options to control parsing.
    :param format_hint: A hint as to the format of the given log files. Either 'auto' (to guess), 'log' (calc log file), 'sir' (digichem result file) or 'sid' (digichem database file).
    :param init_func: An optional function to call to init each newly created process.
    :param processes: The max number of processes to create the new pool object with.
    :param auxiliary_files: An optional list of lists of dicts of auxiliary files. Each item in auxiliary_files should match the corresponding log file in multiple_results.
    :return: A single Result_set object (or child thereof).
    """
    if auxiliary_files is None:
        auxiliary_files = [None] * len(multiple_results)
    
    pool = None
    # Do some parsing.
    # TODO: This parallelization isn't ideal, currently we process each group of to-be merged calcs separately, meaning processes can be wasted.
    try:
        pool = multiprocessing.Pool(processes, initializer = init_func, initargs = init_args if init_args is not None else [])
        
        result_lists = list(
            filterfalse(lambda x: x is None,
                map(partial(multi_merger_parser, options = options, ornt = ornt, ornt_args = ornt_args, parser_options = parser_options, format_hint = format_hint, inner_pool = pool, keep_archive = keep_archive), multiple_results, auxiliary_files)
            )
        )
        
        return result_lists
    
    finally:
        if pool:
            pool.__exit__(None, None, None)
    

class open_for_parsing():
    """
    A context manager for opening log files for parsing. Returns a list of pathlib.Path objects suitable for parsing.
    
    Currently, the main purpose of this context manager is to intelligently handle unpacking of archives (.zip, .tar etc) for parsing.
    """
    
    def __init__(self, *log_files, auxiliary_files = Default(None)):
        log_files = [Path(log_file).resolve() for log_file in log_files]
        
        # Check we haven't been given any duplicate log files.
        # This is just for convenience, if duplicates have been given the user has probably made a mistake.               
        duplicates = [path for path, number in collections.Counter(log_files).items() if number > 1]
        
        for duplicate in duplicates:
            warnings.warn("Ignoring duplicate log file: ".format(duplicate))
        
        # Remove duplicates but retain order.
        self.log_files = list(dict.fromkeys(log_files).keys())

        # Should we worry about duplicate aux files?
        self.has_aux_files = not isinstance(auxiliary_files, Default)
        self.auxiliary_files = {aux_type: Path(auxiliary_file).resolve() for aux_type, auxiliary_file in auxiliary_files.items()} if defres(auxiliary_files) is not None else {}
        
        # A list of tempfile.TemporaryDirectory objects that should be closed when we are finished.
        #self.temp_dirs = []

        self.archive_dirs = {}

        # Keep track of past work to prevent unpacking duplicates.
        self.done_archives = []
    
    @classmethod
    def get_archive_formats(self):
        """
        Get a list of supported archive formats.
        """
        return archive_formats
    
    @property
    def archive_formats(self):
        return archive_formats
        
    def __enter__(self):
        """
        'Open' files for reading.
        """
        return self.open()
    
    @classmethod
    def is_archive(self, path):
        """
        Determine (based on a given file extensions) whether a path points to an archive.
        """
        # We can't simply join all the extensions together, what if the file is named something like "Benzene.log.zip"
        # Also be careful for files like Benzene.tar.gz, this is .tar.gz not .gz
        for pos in range(0, len(path.suffixes)):
            if "".join(path.suffixes[pos:]) in self.get_archive_formats():
                # This is an archive format.                
                return True
        
        return False
    
    def find(self, hint):
        """
        """
        # First, open the file if it is an archive.
        open_files = []

        if hint.exists() and self.is_archive(hint):
            open_files = self.extract(hint)
        
        else:
            open_files = [hint]

        return open_files
        
    def open(self):
        """
        'Open' files for reading.
        """
        new_log_files = []
        
        # First, unpack any explicitly specified aux files.
        # We do this first because we won't extract the same file twice, and we want to make sure
        # we definitely capture the file here.
        new_aux_files = {}
        for aux_type, path in self.auxiliary_files.items():
            if self.is_archive(path):
                files = list(self.extract(path))

                # If the archive contains multiple files, complain (because we don't know which one the user wants).
                if len(files) == 0:
                    raise Digichem_exception("Cannot extract auxiliary file archive '{}'; this archive is empty".format(path))

                elif len(files) > 1:
                    raise Digichem_exception("Cannot extract auxiliary file archive '{}'; this archive contains multiple files".format(path))

                new_aux_files[aux_type] = files[0]
            
            else:
                new_aux_files[aux_type] = path
        
        # Next, build a list of files and folders to check.
        for log_file in self.log_files:
            found_child_archive = False
            
            if log_file.is_dir():
                parent_dir = log_file
                
                # If 'log_file' is a directory, check for an archive inside called 'Output.xxx'.
                for archive_format in self.archive_formats:
                    child_archive = Path(log_file, "Output" + archive_format)
                    if child_archive.exists():
                        if not found_child_archive:
                            # Found an Output dir archive, use this instead.
                            new_log_files.extend(self.find(child_archive))
                            found_child_archive = child_archive
                        
                        else:
                            # For now, only care about the first.
                            warnings.warn("Ignoring subsequent Output archive '{}'; already found '{}'".format(child_archive, found_child_archive))
                
            else:
                # The hint is not a directory, either it is a normal log file, or an archive.
                # If it is an archive, it could either be:
                #  - An archive of a single file (eg, Output.log.zip).
                #  - An archive of an Output dir.
                #  - An archive of a calculation dir (inside of which is Output).
                new_files = list(self.find(log_file))
                new_log_files.extend(new_files)
                parent_dir = log_file.parent
                
                if not all([file.is_file() for file in new_files]):

                #if "Output" in [file.name for file in new_files] or len(list(itertools.chain(*[file.glob("Output") for file in new_files]))) > 0:
                    found_child_archive = True
                
                elif self.is_archive(log_file):
                    # If the hint is a single file archive, also add the parent dir (incase not all of the files are archives).
                    new_log_files.append(log_file.parent)


            if not found_child_archive:
                # If there's not an Output.zip type file, also look for individually zipped output files.
                for sub_file in parent_dir.iterdir():
                    if self.is_archive(sub_file):
                        new_log_files.extend([found.parent for found in self.find(sub_file)])
                        
                new_log_files.extend(self.find(log_file))
        
        if self.has_aux_files:
            return new_log_files, new_aux_files
        
        else:
            return new_log_files
    
    def extract(self, file_name):
        """
        Extract an archive and return the contained log files.
        """
        file_name = file_name.resolve()
        if file_name in self.done_archives:
            digilib.log.get_logger().debug("Skipping duplicate archive '{}'".format(file_name))
            return []

        else:
            self.done_archives.append(file_name)
        
        # Get a temp dir to extract to.
        # We can't use TemporaryDirectory here, because these are auto deleted on program exit. This is not compatible with multi-processing.
        tempdir = mkdtemp()
        self.archive_dirs[file_name] = tempdir
        
        # Extract to it.
        digilib.log.get_logger().info("Extracting archive '{}'...".format(file_name))
        shutil.unpack_archive(file_name, tempdir)
        
        # Add any files/directories that were unpacked.
        return Path(tempdir).glob("*")
    
    def cleanup(self):
        """
        Perform cleanup of any open files.
        
        Alias for __exit__().
        """
        return self.__exit__()
    
    def __exit__(self, exc_type = None, exc_value = None, traceback = None):
        """
        'Close' any open files.
        """
        for tempdir in self.archive_dirs.values():
            shutil.rmtree(tempdir, ignore_errors = True)
        
        