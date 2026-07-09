# Functions and classes for determining file types etc.
from pathlib import Path
from itertools import accumulate

from digichem.exception.base import Unknown_file_type_exception

class File_type():
    """
    Class for representing a known file type.
    """
    
    def __init__(self, name, family = None, extensions = None, short_code = None, exact = None):
        """
        Constructor for File_type objects.
        
        :param name: The name of the file (eg, "JPEG").
        :param family: The name of the family/group that the file type belongs to (eg, "image").
        :param extensions: An iterable of known file extensions (this is case-insensitive) (eg, [".jpg", ".jpeg"]).
        :param exact: An iterable of the file names this file can appear as.
        """
        self.name = name
        self.family = family
        self.extensions = [extension.lower() for extension in extensions] if extensions is not None else []
        self.exact = exact or []
        
        if short_code is None and len(self.extensions) > 0:
            # Take one of our extensions (without the dot).
            self.short_code = self.extensions[0][1:]
            
        else:
            self.short_code = short_code
    
    def check(self, file_path):
        """
        Determine whether a given file is of this file type.
        
        :param file_path: A string-like path to the file to check.
        :return: True if the file is of this type, false otherwise.
        """
        # Get a Path object.
        file_path = Path(file_path)

        # If we match exactly, then nothing more to do.
        for exact in self.exact:
            if file_path.name == exact:
                return True
        
        # Check to see if the file extension matches one of our file extensions.
        # We recurse through all the given file's extensions in case one of our extensions contains multiple parts (".tar.gz" for example).
        for extension in accumulate(reversed(file_path.suffixes), lambda a, b: b+a):
            if extension.lower() in self.extensions:
                return True
        
        # No match found.
        return False
    
    @property
    def file_type(self):
        """
        Get a string uniquely identifying this file type.
        """
        if self.family is not None:
            return "{}/{}".format(self.family, self.name)
        else:
            return self.name
        
    def __str__(self):
        """
        Stringify this file type.
        """
        return self.file_type
    
    
# Known file types.
log_file = File_type("log", "general", [".log"])
gaussian_chk_file = File_type("checkpoint", "gaussian", [".chk"])
gaussian_NTO_chk_file = File_type("NTO-checkpoint", "gaussian", [".chk"])
gaussian_fchk_file = File_type("formatted-checkpoint", "gaussian", [".fchk"])
gaussian_rwf_file = File_type("read-write", "gaussian", [".rwf"])
gaussian_cube_file = File_type("cube", "gaussian", [".cub", ".cube"])

orca_gbw_file = File_type("gbw", "orca", [".gbw"])
orca_density_file = File_type("density", "orca", [".densities"])
orca_density_info_file = File_type("density-info", "orca", [".densitiesinfo"])

xtb_topology = File_type("topo", "xtb", exact = ["xtbtopo.mol"], short_code = "xtbtopo.mol")
xtb_molden = File_type("molden", "xtb", exact = ["molden.input"], short_code = "molden.input")
xtb_density = File_type("density", "xtb", exact = ["density.cub"], short_code = "density.cub")
xtb_std2 = File_type("std2", "xtb", [".std2.log", ".std2.out"])
xtb_xtb4stda = File_type("xtb4stda", "xtb", [".xtb4stda.log", ".xtb4stda.out"])


# A list of all our known types.
known_types = [log_file, gaussian_chk_file, gaussian_fchk_file, gaussian_rwf_file, gaussian_cube_file, orca_gbw_file]

# TODO: This seems to be unused?
def get_file_type(file_path):
    """
    Get the type of a file from a list of known file types.
    
    :raises Unknown_file_type_exception: If the given file is of an unknown type.
    :param file_path: String-like path to the file to check.
    :return: A File_type object.
    """
    # Iterate through our list of known file types.
    for file_type in known_types:
        if file_type.check(file_path):
            return file_type
        
    # Unknown file type.
    raise Unknown_file_type_exception(file_path)

