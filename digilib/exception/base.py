# Exceptions and other errors.

class Digichem_exception(Exception):
    """
    General digichem exception.
    """
    
    def __init__(self, message):
        self.message = message
        
    def __str__(self, *args, **kwargs):
        return self.message
    

class Result_unavailable_error(Digichem_exception):
    """
    Exception for when a requested result is not available (because it could not be found in the calculation results for example).
    """
    
    def __init__(self, result_name, reason = "result could not be found"):
        """
        Constructor for Result_unavailable_error objects.
        
        :param result_name: The name of the result that is unavailable.
        :param reason: Optional message explaining why the result is unavailable. If not given, a default message will be used.
        """
        self.result_name = result_name
        self.reason = reason
        
    def __str__(self, *args, **kwargs):
        """
        Stringify this error.
        """
        return "'{}' is not available; {}".format(self.result_name, self.reason)
        
class File_maker_exception(Digichem_exception):
    """
    Exception for when a file cannot be made/rendered for whatever reason.
    """
    
    def __init__(self, file_maker, reason = ""):
        """
        Constructor for File_maker_exception objects.
        
        :param file_maker: The file_maker object where the exception occurred.
        :param reason: Optional string describing why the exception occurred.
        """
        self.file_maker = file_maker
        self.reason = reason
        
    def __str__(self, *args, **kwargs):
        """
        Stringify this error.
        """
        return "Error making '{}' file '{}'; {}".format(type(self.file_maker).__name__, self.file_maker.output, self.reason)
    
class Unknown_file_type_exception(Digichem_exception):
    """
    Exception for when a file is given but its type cannot be determined.
    """
    
    def __init__(self, file_path, expected = None):
        """
        Constructor for Unknown_file_type_exception objects.
        
        :param file_path: String-like path of the file that is unrecognised.
        :param expected: An optional string-like representing the type of file that was expected. 
        """
        self.file_path = file_path
        self.expected = expected
        
    def __str__(self, *args, **kwargs):
        """
        Stringify this error.
        """
        err_str = "Unknown file type '{}'".format(self.file_path)
        if self.expected is not None:
            err_str = "{}; expected file of type '{}'".format(err_str, self.expected)
        return err_str
                
class Submission_error(Digichem_exception):
    """
    Exceptions for when an error occurs during calculation submission.
    """
    
    def __init__(self, calculation, reason, file_name = None):
        """
        Constructor for Submission_error exception objects.
        
        :param calculation: The calculation that was in process of being submitted when the error occurred. This can be any one of a Method_target, Program_target or Calculation_target.
        :param reason: String describing why the error occurred.
        """        
        # Do some quick type checking.
        if calculation.meta['TYPE'] == "destination":
            # 'Calculation' is actually a Method_target.
            calculation = calculation.program.calculation
        elif calculation.meta['TYPE'] == "program":
            # 'Calculation' is actually a Program_target.
            calculation = calculation.calculation
        
        # Decide on file name.
        if file_name is None:
            file_name = calculation.molecule_name
        self.file_name = file_name
        
        self.calculation = calculation
        self.reason = reason
        
    def __str__(self, *args, **kwargs):
        """
        Stringify this error.
        """
        return "Error submitting file '{}' to '{}'; {}".format(self.file_name, self.calculation.meta['name'], self.reason)
    
    
class Format_error(Digichem_exception):
    """
    Exceptions for when an error occurs during result formatting.
    """
    
    def __init__(self, format, reason):
        """
        Constructor for Format_error exception objects.
        """
        self.format = format
        self.reason = reason
        
    def __str__(self):
        """
        Stringify this error.
        """
        return "{} ({}); {}".format(type(self.format).__name__, ", ".join(getattr(self.format, 'CLASS_HANDLE', [])), self.reason)
    