from pathlib import Path

class Input_file():
    """
    ABC for classes that support various input format types.
    
    In most cases the concrete class digilib.input.Digichem_coords is what you are likely looking for.
    """
    
    def __init__(self, charge = None, multiplicity = None, name = None, file_name = None, history = None):
        """
        Abstract constructor for Input_file classes.
        
        :param charge: An explicit molecular charge (as an integer). If not given, a charge may be interpreted from other properties of the input file.
        :param multiplicity: An explicit molecular multiplicity (as an integer). If not given, a multiplicity may be interpreted from other properties of the input file.
        :param name: An explicit molecule name. If not given, a name may be interpreted from other properties of the input file.
        :param file_name: The path to the file from which this object was created. If this Input_file was not created from a file (eg, from memory instead), this should be None.
        :param history: The SHA of the previous calculation (if applicable) from which these coordinates were taken.
        """
        if charge is not None and not isinstance(charge, int):
            raise TypeError("Charge must be an integer (or None)")
        if multiplicity is not None and not isinstance(multiplicity, int):
            raise TypeError("Multiplicity must be an integer (or None)")
        
        self.charge = charge
        self.multiplicity = multiplicity
        self.name = name
        # TOOD: This should be called file_path
        self.file_name = Path(file_name) if file_name is not None else None
        self.history = history

    def dump(self):
        """
        Get this input file as a dict.
        """
        return {
            'name': self.name,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'history': self.history,
        }
    
    @property
    def implicit_name(self):
        """
        A more intelligent name for this molecule, taking into account our old file name if necessary.
        """
        # If a real name wasn't given, but a file name was, use it.
        if self.name is None and self.file_name is not None:
            name = self.file_name.with_suffix("").name
            
        elif self.name is not None:
            name = self.name
        
        else:
            name = "molecule"
            
        return name

    @property
    def implicit_charge(self):
        """
        The charge of the molecule/system, accounting for cases where no explicit charge is set.
        """ 
        if self.charge is None:
            return 0
        else:
            return self.charge
        
    @property
    def implicit_multiplicity(self):
        """
        The multiplicity (as an integer) of the molecule/system, accounting for cases where no explicit multiplicity is set.
        """ 
        if self.multiplicity is None:
            return 1
        else:
            return self.multiplicity