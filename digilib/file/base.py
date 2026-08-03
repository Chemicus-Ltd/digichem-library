# General imports.
from pathlib import Path
from timeit import default_timer as timer
import datetime
import warnings

from digilib import misc
from digilib.exception.base import File_maker_exception, Digichem_exception

# Digichem imports.
import digilib.log


class File_maker_ABC():
    """
    ABC for file maker classes.
    
    These classes make file(s) on the fly, they are mostly used for image rendering (and making the prerequisite files for image rendering).
    
    Do not try to init an instance of this class yourself; use a real one below (or elsewhere).
    """
    #####################################
    # To be implemented in sub classes. #
    #####################################
    
    def __init__(self, output):
        """
        :param output:  The filename/path to the new file (this path doesn't need to point to a real file yet; we will use this path to write to).
        """
        self.output = output
        self.file_path = {'file': self.output}
        
        # A flag of whether we have already created our output file.
        self.done_file_creation = False
    
    def get_file(self, name = 'file'):
        """
        Get the path to one of the files that this class represents, converting from our input file and writing to file first if necessary.
        """
        raise NotImplementedError("get_file() is not implemented in this ABC")
    
    def delete(self, lazy = False):
        """
        Delete the file(s) described by this object.
        """
        raise NotImplementedError("delete() is not implemented in this ABC")
    
    def get_parallel(self, name = 'file', cpus = 1):
        """
        Generate the file represented by this object in a parallel context.

        get_parallel() will be called by a higher level function as an argument to ThreadPoolExecutor.map() or similar,
        to generate many files simultaneously. This is useful for slow operations (such as cube generation) that are
        difficult to parallelise individually, but easy to parallelise across multiple files.

        :param name: The file to generate.
        :param cpus: The number of CPUs this operation should use, nearly always 1.
        """
        # This default implementation does nothing.
        pass
    
    #################################
    # Implemented in this sub-class #
    #################################
        
    def safe_get_file(self, name = 'file'):
        """
        Get the path to one of the files that this class represents, converting from our input file and writing to file first if necessary.
        
        The functioning of this method is controlled by the dont_modify & use_existing flags.
        
        You can also use the normal python attribute mechanism (either through getattr() or dot notation) to get these paths.
        
        :raises KeyError: If name is not the name of one of the files this class represents.
        :param name: The name of a file to get. Depends on the paths in self.file_path.
        :return: A pathlib Path object pointing to the file represented by name, or None if no file could be created.
        """
        try:
            return self.get_file(name)
        
        except Exception:
            digilib.log.get_logger().error("Unable to get file '{}'".format(self.file_path[name]), exc_info = True)
            
            return None
    
    def relative_path(self, file_name = 'file', *, output_base):
        """
        Get the relative path to the file (or one of the files) represented by this object.
        
        The path is relative to the directory where we write our file, called the output_base.
        
        :param file_name: The name of one of the files that we represent to get the path of, this can be excluded if we only represent one file.
        :param output_base: Pathlib Path to the base directory to get a relative path to.
        :return: A relative Path object to the file.
        """
        # First get the file's real path (this could trigger rendering).        
        file_path = self.safe_get_file(file_name)
        
        # Return None if that doesn't work (eg, because dont_modify = True and use_existing = False) to match the rest of the class.
        if file_path is None:
            return None
                        
        # And return a relative path (this can still raise exceptions).
        return file_path.relative_to(output_base)

    
    def __str__(self):
        return str(self.safe_get_file())


class Dummy_file_maker(File_maker_ABC):
    """
    A dummy file maker that will always fail to make any files.
    
    This class is useful where it is not possible to construct a real file maker object, where instead an instance of this class can be seamlessly used instead.
    """

    def __init__(self, output, message):
        """
        :param output: The path to the file this class is not creating.
        :param message: The reason why this dummy file maker is being used.
        """
        super().__init__(output)
        self.message = message
        
    def get_file(self, name = 'file'):
        """
        Get the path to one of the files that this class represents, converting from our input file and writing to file first if necessary.
        """
        raise File_maker_exception(self, self.message)
    
    def delete(self, lazy = False):
        """
        Delete the file(s) described by this object.
        """
        return False
    

class File_maker(File_maker_ABC):
    """
    Superclass for classes that automatically create files when needed.
    """
    
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "output"
    
    def __init__(self, output, existing_file = None, dont_modify = False, use_existing = False, full_path_names = False):
        """
        Constructor for File_converter objects.
        
        :param output:  The filename/path to the new file (this path doesn't need to point to a real file yet; we will use this path to write to).
        :param existing_file: An optional existing file of the type we're converting to. If this is given, then no conversion is done. This option exists so the user can specify files they have already converted themselves.
        :param dont_modify: Flag that prevents modifying the file on disk. If True, no new file will be written even if one does not already exist. dont_modify is automatically set to True if existing_file is not None (to prevent over-writing whatever file the user gave us).
        :param use_existing: Flag that modifies how file conversion works. If True, existing files will be preferentially used if available (set to False to force overwriting existing files).
        :param full_path_names: Whether to print full path names in log messages.
        """
        # If we've been given an existing_file explicitly, check it exists.
        if existing_file is not None:
            # Convert to path.
            existing_file = Path(existing_file)
            
            # Force set dont_modify = True so we don't overwrite the input that the user has given us.
            dont_modify = True
            digilib.log.get_logger().debug("Setting dont_modify == True to prevent overwrite of user specified {} file '{}'".format(self.output_file_type, existing_file))
            
            # Check
            if not existing_file.exists():
                warnings.warn("The given {} file '{}' does not appear to exist".format(self.output_file_type, existing_file))
                # Continue anyway.
            
            # Set our output appropriately, and continue as normal.
            output = existing_file
        else:
            # Convert to path.
            output = Path(output)
        
        self.dont_modify = dont_modify
        self.use_existing = use_existing
        
        self.failed_file_creation = False
        
        self.full_path_names = full_path_names
        
        super().__init__(output)
    
    def get_file(self, name = 'file'):
        """
        Get the path to one of the files that this class represents, converting from our input file and writing to file first if necessary.
        
        The functioning of this method is controlled by the dont_modify & use_existing flags.
                
        :raises KeyError: If name is not the name of one of the files this class represents.
        :raises File_maker_exception: If the file could not be created for some reason.
        :param name: The name of a file to get. Depends on the paths in self.file_path.
        :return: A pathlib Path object pointing to the file represented by name.
        """
        # We first need to see if we've already made our files or not.
        # Although we can actually represent multiple files, they are all created at the same time.
        if self.done_file_creation:
            return self.file_path[name]
        else:
            # Our files haven't been created. What happens next depends on the options given to us at __init__().
            if self.use_existing:
                # We've been asked to use existing files if they exist, see if they do.
                # There is obviously a race condition here, but we're not making any guarantee to the calling function that the file really exists (or even what it is), we're only checking to see if we need to create a new file or not.
                if self.file_path[name].exists():
                    # File exists, return its path.
                    #digilib.log.get_logger().debug("Using existing {} file '{}'".format(self.output_file_type, self.file_path[name]))
                    return self.file_path[name]
            
            # There are no files we can use. If we're allowed, write new files.
            if not self.dont_modify:
                # Try and make.
                try:
                    self.make()
                    return self.file_path[name]
                
                except File_maker_exception as e:
                    # Don't raise a new exception if we've already got one as this just pollutes the log.
                    self.failed_file_creation = True
                    raise
                
                except Exception as e:
                    # Know that a KeyError could be thrown here if name is not valid, but I think we'd want that to propagate anyway.
                    self.failed_file_creation = True
                    raise File_maker_exception(self, "Unable to create file") from e
            else:
                raise File_maker_exception(self, "Not creating file because dont_modify is True")
    
    @property
    def creation_message(self):
        """
        A short message that may (depending on log-level) be printed to the user before make_files() is called.
        """
        return "Generating {} file '{}'".format(self.output_file_type, self.output if self.full_path_names else self.output.name)
    
    def check_can_make(self):
        """
        Check whether it is feasible to try and create the files(s) that we represent.
        
        Reasons for making not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
        
        This method returns nothing, but will raise an File_maker_exception exception if the rendering is not possible.
        """
        if self.failed_file_creation:
            # We've already tried to make this file once before and failed, don't try again.
            raise File_maker_exception(self, "Not creating file due to previous errors")
    
    def make(self):
        """
        Make the file(s) described by this object.
        
        This method is a wrapper, used to call other methods that should be written by the user.
        
        Most inheriting classes won't need to modify this.
        
        This method should not be called directly, and doing so can overwrite files that we've been asked not to overwrite.
        """
        # First check if we can render. This method will throw an exception to stop us if we can't proceed.
        self.check_can_make()

        # Print a debug message (because lots can go wrong next and this step an be quite slow).    
        digilib.log.get_logger().info(self.creation_message)
        
        start_timer = timer()
        
        # Make our parent folder(s).
        try:
            self.output.parent.mkdir(parents = True)
        except FileExistsError:
            # This will happen when the dir already exists, which is fine.
            pass
            
        
        # Call the real workhorse, which is implementation specific.
        self.make_files()
        
        # Set our flag.
        self.done_file_creation = True
        
        # Print timing info.
        self.duration = datetime.timedelta(seconds = timer() - start_timer)
        digilib.log.get_logger().info("File generation duration: {} ({} total seconds)".format(misc.timedelta_to_string(self.duration), self.duration.total_seconds()))
        
    def delete(self, lazy = False):
        """
        Delete the file(s) described by this object.
        
        This method will try to delete all files represented by this object, even if exceptions are raised during deletion of an earlier file.
        If lazy is not True, note that it is the last caught exception (from the last file deletion that went wrong) that is re-raised.
        
        :param lazy: If true, no exceptions are raised if file deletion is not possible.
        :return: True if all files were deleted successfully, False otherwise.
        """
        # We'll keep hold of any exceptions for later. 
        exc = None
        
        # Delete each of our files.
        for file_name, file_path in self.file_path.items():
            try:
                # First check we're allowed to delete the file.
                if self.dont_modify:
                    raise Digichem_exception("Unable to delete file '{}'; dont_modify is True".format(file_path))
                
                # DELETE.
                file_path.unlink()
                
                # Unset the done_file_creation flag, so we'll re-write any files if they are called for again.
                self.done_file_creation = False
            except Exception as e:
                # Hold onto our exception. This might discard a previous exception, but we can only raise one exception anyway (?).
                exc = e
        
        # If we've been asked to report exceptions, do so.
        if not lazy and exc is not None:
            raise exc
        
        # Also return whether anything bad happened.
        return exc is None
        
    def make_files(self):
        """
        Make the files(s) described by this object. Inheriting classes should write their own implementation.
        """
        pass


class File_converter(File_maker):
    """
    Superclass for classes that automatically convert between different file formats.
    """
    
    # Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
    input_file_type = "input"
    
    
    def __init__(self, output, input_file = None, existing_file = None, dont_modify = False, use_existing = True):
        """
        Constructor for File_converter objects.
        
        :param output:  The filename/path to the new file (this path doesn't need to point to a real file yet; we will use this path to write to).
        :param input_file: The input file that will be converted. Not that while this option is formally optional, conversion will normally fail if no input is available.
        :param existing_file: An optional existing file of the type we're converting to. If this is given, then no conversion is done. This option exists so the user can specify files they have already converted themselves.
        :param dont_modify: Flag that modifies how file conversion works. If True, no new file will be written.
        :param use_existing: Flag that modifies how file conversion works. If True, no new file will be written. If True, existing files will be preferentially used if available (set to False to force overwriting existing files).
        """
        # Call our parent.
        super().__init__(output, existing_file, dont_modify, use_existing)
        
        # And save our input file.
        self.input_file = Path(input_file) if isinstance(input_file, str) else input_file
    
    
    def check_can_make(self):
        """
        Check whether it is feasible to try and create the files(s) that we represent.
        
        Reasons for making not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
        
        This method returns nothing, but will raise an File_maker_exception exception if the rendering is not possible.
        """
        super().check_can_make()
        
        if self.input_file is None:
            raise File_maker_exception(self, "No {} file is available".format(self.input_file_type))
        try:
            if self.input_file.safe_get_file() is None:
                raise File_maker_exception(self, "No {} file is available".format(self.input_file_type))
        except AttributeError:
            # Input file does not have a safe_get_file() method.
            pass
    
    @property
    def creation_message(self):
        """
        A short message that may (depending on log-level) be printed to the user before make_files() is called.
        """
        return "Converting {} file '{}' to {} file '{}'".format(self.input_file_type, self.input_file if self.full_path_names else Path(str(self.input_file)).name, self.output_file_type, self.output if self.full_path_names else Path(str(self.output)).name)        
