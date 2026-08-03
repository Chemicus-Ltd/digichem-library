"""This module contains translation tables between various different, but equivalent, terms expected by different calculation programs."""

import csv

from configurables.misc import is_int

import digilib.log
from digilib.datas import get_resource

# Hidden imports.
#import basis_set_exchange as bse


# Load functional conversion table from file.
functional_db = {}
with open(get_resource('data/functionals.csv')) as csv_file:
    _reader = csv.reader(csv_file)
    headers = None
    for line_num, line in enumerate(_reader):
        if line_num == 0:
            headers = line
        
        else:
            line = [line_part if line_part != "" else None for line_part in line]
            
            # The various parts of the definition.
            functional_name = line[0]
            aliases = line[1].split(",") if line[1] is not None else []
            functional = {"name": functional_name, "aliases": aliases, "turbomole": line[2], "gaussian": line[3]}
            
            # Get all the names we can refer to this functional by.
            all_names = [functional_name] + aliases
            for program in ("turbomole", "gaussian"):
                # Make sure we don't overwrite existing values with program specific ones.
                if functional[program] is not None and functional[program] not in functional_db:
                    all_names.append(functional[program])
            
            # Set into the DB, using all possible names.
            for ref_name in all_names:
                functional_db[ref_name] = functional


# Load solvent conversion table from file.
solvent_db = {}
with open(get_resource('data/solvents.csv')) as csv_file:
    _reader = csv.reader(csv_file)
    headers = None
    for line_num, line in enumerate(_reader):
        if line_num == 0:
            headers = line
        
        else:
            line = [line_part if line_part != "" else None for line_part in line]
            
            # The various parts of the definition.
            solvent = {
                "name": line[0],
                "aliases": [line[1]] if line[1] is not None else [],
                "gaussian": line[2],
                "epsilon": float(line[3]) if line[3] is not None else None,
                "refractive": float(line[4]) if line[4] is not None else None,
                "xTB": line[5],
            }
            
            solvent_db[solvent['name']] = solvent


class Translate():
    """
    ABC for translate classes.
    """
    
    def __init__(self, value):
        """
        Constructor for Translate classes.
        """
        if isinstance(value, type(self)):
            self.value = value.value
        
        else:
            self.value = value
    
    def translate(self, program):
        """
        Translate into a name appropriate for a given program.
        """
        raise NotImplementedError("Implement in subclass")
    
    def to_turbomole(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        return self.translate("turbomole")
    
    def to_gaussian(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        return self.translate("gaussian")
    
    def to_orca(self):
        """
        Translate into a name appropriate for ORCA.
        """
        return self.translate("orca")


class Basis_set(Translate):
    """
    A class for converting between the names of basis sets.
    """
    
    @property
    def basis_set(self):
        return self.value
    
    @classmethod
    def build_db(self):
        """
        Build a table of basis set names.
        
        :returns: A dictionary, where each key is the name of a basis set and each value is a BSE metadata dict.
        """
        import basis_set_exchange as bse
        # Build a database of basis set metadata, adding the display name and alternative names as additional keys for each basis set.
        # BSE's database has some semi-duplicate entries (such as '6-31g(d,p)' and '6-31g_st__st_').
        # Although these basis sets are identical, they have reversed metadata (ie, for '6-31g(d,p)', 'display_name' is set to '6-31G(d,p)' and 'other_names'
        # contains '6-31G**', while for '6-31g_st__st_' 'display_name' is '6-31G**' and 'other_names' is '6-31G(d,p)'.
        # This is a problem for us because we need to predictably choose from these names (ie the star form '6-31G**' is suitable for both
        # Gaussian and Turbomole), but the keys referring to these same names are different.
        #
        # Firstly, build a new dict using basename as the key, which appears to be the same even for duplicates.
        db = {value["basename"]: value for key, value in bse.get_metadata().items()}
        
        new_db = {}
        
        for basis_key, basis_set_meta in db.items():
            # First, add the normal key.
            basis_set = {
                "key": basis_set_meta['basename'],
                "name": basis_set_meta['display_name'],
                "aliases": basis_set_meta['other_names'],
                "meta": basis_set_meta
            }
            
            # Decide on our program specific names.
            basis_set["gaussian"] = basis_set['name']
            basis_set["turbomole"] = basis_set['name']
            
            # Turbomole doesn't like pople sets with polarisation function spelled out,
            # it prefers the star format.
            if basis_set['name'][-5:] == "(d,p)":
                basis_set["turbomole"] = basis_set['name'][:-5] + "**"
                
            # Gaussian has a strange, contracted style naming scheme for Karlsruhe,
            # and a misleading/incorrect name for def2-SV(P).
            #print(basis_set['name'])
            if basis_set['name'] == "def2-SV(P)":
                basis_set['gaussian'] = "def2SVPP"
                
            elif basis_set['name'][:5] == "def2-":
                basis_set['gaussian'] = basis_set['name'][:4] + basis_set['name'][5:]
            
            new_db[basis_key] = basis_set
            
            # TODO: Need to support aux basis sets (which end in things like -rifit, -jfit, -jkfit etc.
            # Orca, for example, uses /C, /J and /JK instead...
            # Turbomole uses nothing.
        
        # Basis set exchange seems to be missing some non-polarized Karlsruhe basis sets (or else Gaussian made them up).
        # Add them manually so we can still convert from Gaussian's weird representation.
        new_db.update({
            "def2-TZV": {"key": "def2-TZV", "name": "def2-TZV", "aliases": [], "meta": {}, "gaussian": "def2TZV"},#, "turbomole": "def2-TZV"},
            "def2-QZV": {"key": "def2-QZV", "name": "def2-QZV", "aliases": [], "meta": {}, "gaussian": "def2QZV"}#, "turbomole": "def2-TZV"}
        })
                
        return new_db
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find the entry for a basis set in the basis set exchange.
        
        :raises ValueError: If the given hint could not be found.
        :param hint: The name of the basis set to search for.
        :return: A dictionary of BSE metadata.
        """
        # Get the basis set DB.
        db = self.build_db()
        
        # First, try and see if we can just use hint as an exact match.
        try:
            return db[hint]
            
        except KeyError:
            pass
        
        # Try again, this time ignoring case and looking through all possible names.
        str_hint = "".join(str(hint).lower().split())
        for basis_set in db.values():
            names = [basis_set.get(name, "") for name in ("key", "name", "gaussian", "turbomole")] + basis_set['aliases']
            if str_hint in [name.lower() for name in names if name is not None]:
                return basis_set
        
        # No luck.
        raise ValueError("Could not find basis set definition for '{}'".format(hint))
    
    def translate(self, program = "name"):
        """
        Translate the name of this basis set into one appropriate for a calculation program.
        
        Basis sets are currently handled the same for all programs.
        """
        try:
            basis_def = self.find_in_db(self.basis_set)
            return basis_def[program]
        
        except ValueError:
            # Just return as is.
            if self.basis_set != "auto":
                digilib.log.get_logger().debug("Could not find basis set with name '{}' in the basis set exchange; using name unmodified".format(self.basis_set))
            return self.basis_set
    
        except KeyError:
            return basis_def['name']
    
    def __str__(self):
        return str(self.translate())


class Functional(Translate):
    """
    A class for converting between the names of DFT functionals.
    """
    
    @property
    def functional(self):
        return self.value
            
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a functional in the internal library.
        
        :raises ValueError: If the functional could not be found.
        :param hint: The name of the functional to look for.
        :returns: The corresponding functional dict.
        """
        # First try using exact name.
        try:
            return functional_db[hint]
            
        except KeyError:
            pass
        
        # Try again, this time ignoring case.
        functional_names = [functional_name.upper() for functional_name in functional_db]
        
        try:
            return functional_db[list(functional_db.keys())[functional_names.index(hint.upper())]]
        
        except IndexError:
            pass
        
        # No luck.
        raise ValueError("Could not find functional definition for  '{}'".format(hint))
        
    
    def translate(self, program):
        """
        Translate this functional name to one recognised by a given program.
        
        :param program: The program to translate for.
        """
        # Try and get a definition for the functional from our db.
        try:
            func_def = self.find_in_db(self.functional)
                
            # If there's an explicit value for our program, use that.
            if program in func_def and func_def[program] is not None:
                return func_def[program]
            
            # Otherwise, use the main common name.
            return func_def['name']
            
        except ValueError:
            pass
        
        # No definition found, use as is.
        return self.functional
    
    def to_turbomole(self):
        """
        Translate into a name appropriate for Turbomole.
        """
        func_name = self.translate("turbomole")
        # All turbomole functionals are lower case, except for the word 'Gaussian'.
        return func_name.lower().replace("gaussian", "Gaussian")
    
    def __str__(self):
        return str(self.translate("name"))


class Solvent(Translate):
    """A class for converting between different representations of solvents."""
    
    @property
    def solvent(self):
        return self.value
    
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a solvent in the internal library.
        
        :raises ValueError: If the multiplicity could not be found.
        :param hint: The name, number or symbol of the multiplicity to look for.
        :returns: The corresponding multiplicity dict.
        """
        # First, just try to lookup exactly.
        try:
            return solvent_db[hint]
            
        except KeyError:
            # No luck, look through for a match.
            str_hint = "".join(str(hint).lower().split())
            for solvent in solvent_db.values():
                if hint == solvent['epsilon'] or str_hint in [name.lower() for name in solvent['aliases'] + [solvent['name']] + [solvent['gaussian']] if name is not None]:
                    return solvent
        
        # No luck.
        raise ValueError("Could not find solvent definition for '{}'".format(hint))
    
    @classmethod
    def epsilon_to_name(self, epsilon, threshold = 0.0001):
        """
        Find the name of a solvent based on its epsilon value.
        """
        close = {}
        for solvent_def in solvent_db.values():
            if solvent_def['epsilon'] == epsilon:
                # Exact match, stop.
                return solvent_def['name']
            
            elif abs(solvent_def['epsilon'] - epsilon) < threshold:
                # Close match.
                close[abs(solvent_def['epsilon'] - epsilon)] = solvent_def['name']
                
        # Sort values that were close based on how close they were.
        close_names = [close[closeness] for closeness in sorted(close)]
        
        try:
            return close_names[0]
        
        except IndexError:
            raise ValueError("No solvent definitions within {} of {}".format(threshold, epsilon)) from None
    
    @property
    def epsilon(self):
        return self.find_in_db(self.solvent)['epsilon']
    
    @property
    def refractive_index(self):
        return self.find_in_db(self.solvent)['refractive']
    
    def translate(self, program):
        """
        Translate this solvent name to one recognised by a given program.
        
        :param program: The program to translate for.
        """
        if program == "turbomole":
            program = "epsilon"
            
        elif program == "orca":
            program = "name"
        
        try:
            solvent_def = self.find_in_db(self.solvent)
            
            if program in solvent_def and solvent_def[program] is not None:
                return solvent_def[program]
            
            else:
                return solvent_def["name"]
        
        except ValueError:
            return self.solvent
        
    def __str__(self):
        return str(self.translate("name"))


class SCF_convergence(Translate):
    """
    A class for converting between different shortcuts for SCF convergence thresholds.
    """
    
    table = [
        {"name": "Loose",       "energy": -5,       "density": -3,      "orca": "LooseSCF"},   # Note energy stops changing here.
        {"name": "Weak",        "energy": -5,       "density": -4,      "orca": "SloppySCF"},                     
        {"name": "Medium",      "energy": -6,       "density": -5,      "orca": "MediumSCF"},
        {"name": "Strong",      "energy": -7,       "density": -6,      "orca": "StrongSCF"},
        {"name": "Tight",       "energy": -8,       "density": -7,      "orca": "TightSCF"},   # Typically a sensible default.
        {"name": "VTight",      "energy": -9,       "density": -8,      "orca": "VeryTightSCF"},
        {"name": "VVTight",     "energy": -10,      "density": -9,      "orca": "VeryVeryTightSCF"},
        {"name": "Extreme",     "energy": -14,      "density": -14,     "orca": "ExtremeSCF"},
    ]
    
    @classmethod
    def names(self):
        return [row['name'] for row in self.table]
    
    @classmethod
    def choices(self):
        return [self(row['name']) for row in self.table]

    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry in the internal library.
        
        :raises ValueError: If the value could not be found.
        :param hint: The name to look for
        :returns: The corresponding value dict.
        """
        for row in self.table:
            name, energy, density, orca = row.values()
            
            if str(hint).upper() == name.upper() or str(hint).upper() == orca.upper():
                return row
            
            if str(hint) == str(energy) or str(hint).upper() == "ENERGY{}".format(energy) or str(hint).upper() == "DENSITY{}".format(density):
                return row
            
        # No luck.
        raise ValueError("Could not find convergence definition for '{}'".format(hint))
    
    def translate(self, to_type):
        """
        Translate into a name appropriate for a given program.
        """
        try:
            return self.find_in_db(self.value)[to_type]
        
        except ValueError:
            # Couldn't find in conversion table.
            return self.value
        
    def __str__(self):
        return str(self.translate("name"))
    
    def __eq__(self, other):
        return str(self) == str(other)
    

class Cube_grid_points(Translate):
    """A class for converting between cube grid sizes."""
    
    table = [
        {"name": "Default",     "points": 100,      "gaussian": 0,       "turbomole": "m3",      "orca": 100},   # Default depends on the calc program. Gaussian uses a special value of 'zero' to select a default algorithm.
        {"name": "Tiny",        "points": 25,                            "turbomole": "1",                  },
        {"name": "Small",       "points": 50,                            "turbomole": "2",                  },
        {"name": "Medium",      "points": 100,                           "turbomole": "3",                  },
        {"name": "Large",       "points": 200,                           "turbomole": "4",                  },
        {"name": "Huge",        "points": 500,                           "turbomole": "5",                  },
        
    ]

    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry in the internal library.
        
        :raises ValueError: If the value could not be found.
        :param hint: The name to look for
        :returns: The corresponding value dict.
        """
        for row in self.table:
            if str(hint).upper() == row['name'].upper():
                return row
            
        # No luck.
        raise ValueError("Could not find grid point definition for '{}'".format(hint))
    
    def translate(self, to_type):
        """
        Translate into a name appropriate for a given program.
        """
        try:
            grid_def = self.find_in_db(self.value)
            return grid_def[to_type]
        
        except KeyError:
            # Just return the equivalent number of points.
            return grid_def['points']
        
        except ValueError:
            # Couldn't find in conversion table.
            return self.value
        
    def __str__(self):
        return str(self.translate("name"))


# TODO: Use this.
class Multiplicity(Translate):
    """A class for converting between different representations of multiplicity."""
    
    table = [
        {"name": "no multiplicity", "number": 0, "symbol": "?"},
        {"name": "Singlet", "number": 1, "symbol": "S"},
        {"name": "Doublet", "number": 2, "symbol": "D"},
        {"name": "Triplet", "number": 3, "symbol": "T"},
        {"name": "Quartet", "number": 4, "symbol": "Q"}
    ]
    
    @property
    def multiplicity(self):
        return self.value
        
    @classmethod
    def find_in_db(self, hint):
        """
        Try and find an entry for a multiplicity in the internal library.
        
        :raises ValueError: If the multiplicity could not be found.
        :param hint: The name, number or symbol of the multiplicity to look for.
        :returns: The corresponding multiplicity dict.
        """
        for row in self.table:
            name, number, symbol = row.values()
            if hint == number or str(hint).upper() in (str(number), name.upper(), symbol.upper()):
                return row
            
        # No luck.
        raise ValueError("Could not find multiplicity definition for '{}'".format(hint))
    
    @property
    def symbol(self):
        """
        Get this multiplicity as a symbol.
        """
        # Get a shorthand symbol if we can.
        try:
            return self.find_in_db(self.multiplicity)['symbol']
            
        except ValueError:
            if self.multiplicity % 1 == 0:
                # Multiplicity is an integer, so return as a stringy whole number.
                return str(int(self.multiplicity))
            else:
                return str(self.multiplicity)
            
    @property
    def string(self):
        """
        Get this multiplicity as a string.
        """
        return self.translate("name")
    
    def __str__(self):
        return self.string
    
    @property
    def number(self):
        """
        Get this multiplicity as a number.
        """
        try:
            return self.find_in_db(self.multiplicity)['number']
        
        except ValueError:
            # No pre-defined number, see if it is an int.
            if is_int(self.multiplicity):
                return int(self.multiplicity)
            
            else:
                return float(self.multiplicity)
    
    def translate(self, to_type):
        """
        Translate into a name appropriate for a given program.
        """
        if to_type == "gaussian":
            to_type = "name"
        
        elif to_type == "turbomole":
            to_type = "number"
        
        
        try:
            return self.find_in_db(self.multiplicity)[to_type]
        
        except ValueError:
            # Couldn't find in conversion table.
            return self.multiplicity
