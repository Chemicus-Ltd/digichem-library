import warnings
import logging.handlers
import sys
import textwrap
from subprocess import CalledProcessError

from digichem.exception import Digichem_exception

# The handler object digichem uses for logging.
LOGGING_HANDLER = None

# The name of the digichem logger, can be passed to logging.get_logger() to get the digichem logging object.
LOGGER_NAME = "digichem"


def init_logger(file_name = None, time = False):
    """
    Init the package wide logger.
    
    :param file_name: Optional file to write to. If not given, messages are logged to stderr.
    :param time: Whether to include the time in the logging message.
    """
    global LOGGING_HANDLER
    
    logging.captureWarnings(True)
    
    logger = logging.getLogger(LOGGER_NAME)
    warnings_logger = logging.getLogger("py.warnings")
    
    # Choose our handler.
    if file_name is None:
        # STDERR logging.
        # The console handler, where we'll print most messages.
        LOGGING_HANDLER = Handler(sys.stderr)
    
    else:
        # File logging.
        LOGGING_HANDLER = logging.handlers.RotatingFileHandler(file_name, maxBytes = 1024, backupCount = 5)
        
    # Handle everything.
    LOGGING_HANDLER.setLevel(logging.DEBUG)
    # Set its formatter.
    var_formatter = Variable_formatter(logger, show_time = time, default_warning_formatter = warnings.formatwarning)
    LOGGING_HANDLER.setFormatter(var_formatter)
    
    # Remove old handlers.
    loggers = (logger, warnings_logger)
    for log in loggers:
        while len(log.handlers) > 0:
            log.removeHandler(log.handlers[0])
            
        log.addHandler(LOGGING_HANDLER)
    
    # Add the handler.
    warnings.formatwarning = var_formatter.formatWarning
    
    
def set_logging_level(log_level, verbose = None):
    """
    Set the logging level of the digichem logger object.
    
    :param log_level: The base logging level as a string.
    :param verbose: An integer which specifies how much to increase the logging level by. If verbose is zero (or None) then the logging level is determined only by log_level.
    """
    logger = get_logger()
    
    # Set from log_level first.
    if log_level == "OFF":
        logger.setLevel(60)
    else:
        logger.setLevel(log_level)
    
    # Now adjust with verbosity.
    if verbose is not None:
        # Set from verbosity.
        new_level = logger.level - verbose * 10
        
        # Don't allow us to reach 0 (because this is actually 'UNSET').
        if new_level <= 0:
            new_level = 10
        
        # And set.
        logger.setLevel(new_level)
        
#     # Adjust warnings if necessary.
#     if logger.level == 10:
#         warnings.simplefilter('always', DeprecationWarning)


def get_logger():
    """
    Get the logger used by all parts of digichem.
    """
    return logging.getLogger(LOGGER_NAME)


class Handler(logging.StreamHandler):
    """
    """


class Variable_formatter(logging.Formatter):
    """
    The logging formatter used by digichem, the format changes depending on a logger's log_level.
    """
    
    # Different formatters for printing the message.
    DEFAULT_FORMATTER = '%(handle)s: %(levelname)s: %(message)s'
    WHEN_FORMATTER = '%(handle)s: %(asctime)s: %(levelname)s: %(message)s'
    
    # Format string for printing the date/time.
    DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    
    def __init__(self, logger, show_time = False, *, default_warning_formatter):
        super().__init__(
            #fmt = self.WHEN_FORMATTER if show_time else self.DEFAULT_FORMATTER,
            datefmt = '%Y-%m-%d %H:%M:%S',
            style = '%'
        )
        # Save our logger.
        self.show_time = show_time
        self.logger = logger
        self.default_warning_formatter = default_warning_formatter
    
    def formatMessage(self, record):
        #return logging.PercentStyle(self.WHEN_FORMATTER if self.show_time else self.DEFAULT_FORMATTER).format(record)
        if record.name in ['digichem.silico', 'py.warnings']:
            handle = 'digichem'
        
        else:
            handle = record.name.split('.')[-1]

        msg = handle + ": "
        
        if self.show_time:
            msg += record.asctime + ": "

        target_len = 16

        msg += "{: >{}}".format(getattr(record, 'stream_name', record.levelname) + ": ", target_len - len (msg))

        msg += record.getMessage()

        return msg

    
    # def format(self, record):
    #     ""
        
        
    #     return super().format(record)
        
    def formatWarning(self, message, category, filename, lineno, line=None):
        """
        Format a warning (from the warnings module).
        
        :return: The formatted warning message.
        """
        if self.logger.getEffectiveLevel() > logging.DEBUG:
            return "{}: {}".format(category.__name__, message)
        else:
            return self.default_warning_formatter(message, category, filename, lineno, line=None)
                
    def formatException(self, exc_info):
        """
        Format an exception.
        
        In Variable_formatter, how we do this depends on our log level.
        
        :return: The formatted exception.
        """
        if isinstance(exc_info[1], Digichem_exception) and self.logger.getEffectiveLevel() > logging.DEBUG:
            return self.exception_to_str(exc_info[1])
        else:
            return textwrap.indent(super().formatException(exc_info), "  ")
        
    @classmethod
    def process_output_to_str(self, process_exception):
        """
        Get additional descriptive text from a CalledProcessError exception.
        
        
        """
        output = ""
        if process_exception.stdout is not None and process_exception.stdout != "\n":
            output += "\n" + process_exception.stdout.strip()
        if process_exception.stderr is not None and process_exception.stderr != "\n":
            output += "\n" + process_exception.stderr.strip()
        return output
        
    @classmethod
    def exception_to_str(self, exception):
        """
        Recursively get a string representation of an exception.
        
        :return: A string representation of exception and any prior exceptions.
        """
        excstr = textwrap.indent("{}: {}".format(type(exception).__name__, exception), "\t")
        
        # If the exception is a CalledProcessError, we'll also append the stdout/stderr.
        if isinstance(exception, CalledProcessError):
            excstr += textwrap.indent(self.process_output_to_str(exception), "\t\t")
        
        if exception.__cause__ is not None:
            excstr += "\n{}".format(self.exception_to_str(exception.__cause__))
        elif exception.__context__ is not None and not exception.__suppress_context__:
            excstr += "\n{}".format(self.exception_to_str(exception.__context__))
        
        return excstr