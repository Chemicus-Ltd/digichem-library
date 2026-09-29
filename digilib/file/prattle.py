# General imports.
import subprocess
from subprocess import CalledProcessError
import re
import os
import copy
from pathlib import Path
import warnings
import json

# Digichem imports.
from digilib.exception.base import Digichem_exception
import digilib.log


class Openprattle_converter():
    """
    Provides an interface to oprattle.

    Openprattle provides a library interface, but because of the GPL we cannot use it.
    Calling the oprattle command is fine, however.
    """

    def __init__(self, *, input_file = None, input_file_buffer = None, input_file_path = None, input_file_type = None, executable = "oprattle"):
        """
        Constructor for the OpenBabel converter.
        
        Input files can be specified in one of three ways:
         - As an open file descriptor (input_file and input_file_type)
         - As an in-memory buffer, most probably a string or bytes-like object (input_file_buffer and input_file_type)
         - As a file path (input_file_path and optionally input_file_type)
        
        :param input_file: An open file descriptor in the format given by input_file_type that should be converted.
        :param input_file_buffer: Alternatively, a buffer (unicode string or bytes) in the format given by input_file_type that should be converted.
        :param input_file_path: Alternatively, a path to a file that should be converted.
        :param input_file_type: A shortcode identifying the format of the input file. If not given but input_file_path is given, then this will be determined automatically.
        :param executable: Path or command name to the oprattle executable.
        """
        self.input_file = input_file
        self.input_file_buffer = input_file_buffer if input_file is None else input_file.read()
        self.input_file_path = input_file_path
        self.input_file_type = input_file_type
        self.executable = executable
        
    @property
    def input_name(self):
        """
        A descriptive name of the file we are converting. Works even if converting from memory.
        """
        if self.input_file_path is not None:
            return self.input_file_path
        else:
            return "(file loaded from memory)"
        
    @classmethod
    def get_cls(self, input_file_type):
        """
        This function is deprecated, and does nothing.
        """
        warnings.warn("get_cls() is deprecated, use Openprattle_converter directly instead", DeprecationWarning)
        return self
        
    @classmethod
    def from_file(self, input_file_path, input_file_type = None, **kwargs):
        """
        This function is deprecated, and does nothing.
        """
        warnings.warn("from_file() is deprecated, use the Openprattle_converter constructor directly instead", DeprecationWarning)
        return self.get_cls(None)(input_file_path = input_file_path, input_file_type = input_file_type, **kwargs)
    
    @classmethod
    def type_from_file_name(self, input_file_name, allow_none = False):
        """
        Get the type of a file based on its file name.
        
        This method largely uses the file extension (.com, .tmol etc), with a few other simple rules.
        
        :param input_file_name: The name of the file to check.
        :param allow_none: If the type of the file cannot be determined and allow_none is False (the default), an exception is raised. Otherwise, None is returned.
        """
        try:
            input_file_name = Path(input_file_name)
        
        except TypeError:
            if allow_none:
                return None
            
            else:
                # No file given.
                raise ValueError("Could not automatically determine file format; no file name was given") from None
        
        # Get file extension (removing the dot character).
        extension = input_file_name.suffix[1:]
        
        if extension != "":
            # All done.
            return extension.lower()
        
        elif input_file_name.name == "coord":
            # This is a turbomole file.
            return "tmol"
        
        else:
            if allow_none:
                return None
            
            else:
                # Don't recognise the file format.
                raise ValueError("Could not determine file format of file '{}'; the file does not have an extension and is not recognised".format(input_file_name))
        
    def convert(self, output_file_type = None, output_file = None, *, gen3D = None, charge = None, multiplicity = None):
        """
        Convert the input file wrapped by this class to the designated output_file_type.
         
        :param output_file_type: The file type to convert to.
        :param output_file: Optional file name to write to. If not given, the converted file will be returned as a string (or binary string depending on format).
        :param gen3D: If True and the loaded molecule does not have 3D coordinates, these will be generated (this will scramble atom coordinates).        
        :param charge:  Optional charge of the output format.
        :param multiplicity: Optional multiplicity of the output format.
        :return: The converted file, or None if output_file is not None.
        """
        output_file = str(output_file) if output_file is not None else None
        
        # The signature we'll use to run oprattle.
        sig = [
            str(self.executable),
            "--json"
        ]
        
        # Add the input path if we're reading from file.
        if self.input_file_buffer is None:
            sig.append(str(self.input_file_path))
            
        # Now add the input and output switches.
        if output_file_type:
            sig.extend([
                "-o", output_file_type
            ])
        if self.input_file_type:
            sig.extend([
                "-i", self.input_file_type
            ])
        
        # Add gen3D command if we've been asked to.
        if gen3D is True:
            sig.extend(["--gen3D", "True"])
        elif gen3D is False:
            sig.extend(["--gen3D", "False"])

        # Add charge and mult.
        if charge is not None:
            sig.extend(["--charge", str(charge)])
        
        if multiplicity is not None:
            sig.extend(["--multiplicity", str(multiplicity)])
            
        # If a file to write to has been given, set it.
        if output_file is not None:
            sig.extend(['-O', output_file])
        
        # Give our input_file as stdin if we're not reading from file.
        inputs = self.input_file_buffer

        # Encode strings.
        if isinstance(inputs, str):
            inputs = inputs.encode()
        
        # GO.
        done_process = subprocess.run(
            sig,
            input = inputs,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
        )

        # This can throw exceptions.
        self.handle_logging(done_process.stderr.decode())

        if done_process.returncode != 0:
            raise Digichem_exception("prattle subprocess returned code {}".format(done_process.returncode))
        
        # Return our output.
        return done_process.stdout.decode() if output_file is None else None
    
    def handle_logging(self, raw_output):
        """
        """
        exceptions = []
        for raw_message in raw_output.split("\n"):
            if raw_message == "":
                # Nothing returned, nothing to do.
                continue
            
            # Each message should be in JSON, but check.
            try:
                message = json.loads(raw_message)
                message_text = message['message']
                if message['exception']:
                    exceptions.append(Exception(message['exception']))
                    continue
                    #message_text += "\n" + message['exception']
                digilib.log.get_logger().log(
                    message['levelno'],
                    message_text
                )

            except Exception:
                digilib.log.get_logger().error("Unexpected output from oprattle: '{}'".format(raw_message), exc_info=False)

        if len(exceptions) > 0:
            raise exceptions[0]
    


class Oprattle_formats():
    """
    Class for retrieving the supported file formats from oprattle.
    """
    
    def __init__(self, executable = "oprattle"):
        """
        """
        self.executable = executable
    
    def run(self, readwrite):
        """
        """
        if readwrite not in ["read", "write"]:
            raise ValueError("readwrite must be one of either 'read' or 'write'")
        
        # The signature we'll use to run oprattle.
        sig = [
            self.executable,
            "--json",
            "--readable" if readwrite == "read" else "--writable"
        ]
        
        # GO.
        done_process = subprocess.run(
             sig,
             stdout = subprocess.PIPE,
             stderr = subprocess.PIPE,
             universal_newlines = True,
             check = True,
         )
        
        formats = json.loads(done_process.stdout)
        
        return formats

    
    def read(self, refresh = False):
        """
        Retrieve supported input (read) formats.
        """
        if refresh:
            try:
                del self._read
            
            except AttributeError:
                pass
                
        try:
            return self._read
        
        except AttributeError:
            # Cache miss.
            self._read = self.run("read")
            return self._read
        
    def write(self, refresh = False):
        """
        Retrieve supported input (read) formats.
        """
        if refresh:
            try:
                del self._write
            
            except AttributeError:
                pass
        
        try:
            return self._write
        
        except AttributeError:
            # Cache miss.
            self._write = self.run("write")
            return self._write
    
    