# Classes for representing MOs.

# General imports.
from itertools import filterfalse, zip_longest, chain
import warnings
import math

from digilib.exception import Result_unavailable_error
from digilib.result import Result_container
from digilib.result import Result_object
from digilib.result import Floatable_mixin


class Molecular_orbital_list(Result_container):
    """
    Class for representing a group of MOs.
    """
    
    # A warning issued when attempting to merge non-equivalent orbital lists.
    MERGE_WARNING = "Attempting to merge lists of orbitals that are not identical; non-equivalent orbitals will be ignored"
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    @property
    def HOMO_energy(self):
        """
        Get the energy of the highest occupied orbital in this list.
        
        Use get_orbital(HOMO_difference = 0) to retrieve the HOMO as an object.
        
        :raises Result_unavailable_error: If the HOMO is not available.
        :return The orbital energy (in eV).
        """
        return self.get_orbital(HOMO_difference = 0).energy
    
    @property
    def LUMO_energy(self):
        """
        Get the energy of the lowest unoccupied orbital in this list.
        
        Use get_orbital(HOMO_difference = 1) to retrieve the LUMO as an object.
        
        :raises Result_unavailable_error: If the LUMO is not available.
        :return The orbital energy (in eV).
        """
        return self.get_orbital(HOMO_difference = 1).energy
        
    @property
    def HOMO_LUMO_energy(self):
        """
        Get ΔE HOMO-LUMO; the energy difference between the HOMO and LUMO.
        
        Depending on how the orbitals are occupied, it might not make sense to calculate the HOMO-LUMO energy.
        
        :raises Result_unavailable_error: If the HOMO or LUMO is not available.
        :return: The HOMO-LUMO energy gap (in  eV).
        """
        # Get out orbitals.
        HOMO = self.get_orbital(HOMO_difference = 0)
        LUMO = self.get_orbital(HOMO_difference = 1)
        # Return the difference.
        return float(LUMO) - float(HOMO)
    
    @property
    def occupied(self):
        """
        Return a new Molecular_orbital_list containing only the occupied orbitals of this list.
        """
        return type(self)([orbital for orbital in self if orbital.is_occupied])
    
    @property
    def virtual(self):
        """
        Return a new Molecular_orbital_list containing only the virtual (unoccupied) orbitals of this list.
        """
        return type(self)([orbital for orbital in self if not orbital.is_occupied])

    @property
    def spin_type(self):
        """
        Get the spin type (alpha, beta etc.) of the orbitals in this list.
        
        :raises Result_unavailable_error: If there are no orbitals in this list.
        :return: The spin type, one of either 'alpha' or 'beta' for unrestricted calcs, 'none' for restricted calcs, or 'mixed' if multiple spin types are present in this list.
        """
        # Unique spin types in our list.
        spin_types = list(set([orbital.spin_type for orbital in self]))
        
        if len(spin_types) == 1:
            return spin_types[0]
        elif len(spin_types) > 1:
            return "mixed"
        else:
            raise Result_unavailable_error("Orbital Spin", "There are no orbitals")
    
    def get_orbital(self, criteria = None, *, label = None, HOMO_difference = None, level = None):
        """
        Retrieve an orbital based on some property.
        
        Only one of the criteria should be specified.
        
        :deprecated: Use find() instead.
        :raises Result_unavailable_error: If the requested MO could not be found.
        :param criteria: A string describing the orbital to get. The meaning of criteria is determined automatically based on its content. If criteria begins with '+' or '-' then it is used as a HOMO_difference. Otherwise, if criteria is a valid integer, then it is used as a level. Otherwise, criteria is assumed to be an orbital label.
        :param label: The label of the orbital to get.
        :param HOMO_difference: The distance from the HOMO of the orbital to get.
        :param level: The level of the orbital to get.
        :return: The Molecular_orbital object.
        """
        warnings.warn("get_orbital is deprecated, use find() instead", DeprecationWarning)
        return self.find(criteria, label = label, HOMO_difference = HOMO_difference, level = level)
    
    def find(self, criteria = None, *, label = None, HOMO_difference = None, level = None):
        """
        Retrieve an orbital based on some property.
        
        Only one of the criteria should be specified.
        
        :raises Result_unavailable_error: If the requested MO could not be found.
        :param criteria: A string describing the orbital to get. The meaning of criteria is determined automatically based on its content. If criteria begins with '+' or '-' then it is used as a HOMO_difference. Otherwise, if criteria is a valid integer, then it is used as a level. Otherwise, criteria is assumed to be an orbital label.
        :param label: The label of the orbital to get.
        :param HOMO_difference: The distance from the HOMO of the orbital to get.
        :param level: The level of the orbital to get.
        :return: The Molecular_orbital object.
        """
        # Use the search() method to do our work for us.
        return self.search(criteria, label = label, HOMO_difference = HOMO_difference, level = level, allow_empty = False)[0]
    
    def search(self, criteria = None, *, label = None, HOMO_difference = None, level = None, allow_empty = True):
        """
        Attempt to retrieve a number of orbitals based on some property.
        
        Only one of the criteria should be specified.
        
        :param criteria: A string describing the orbitals to get. The meaning of criteria is determined automatically based on its content. If criteria begins with '+' or '-' then it is used as a HOMO_difference. Otherwise, if criteria is a valid integer, then it is used as a level. Otherwise, criteria is assumed to be an orbital label.
        :param label: The label of the orbitals to get.
        :param HOMO_difference: The distance from the HOMO of the orbitals to get.
        :param level: The level of the orbitals to get.
        :param allow_empty: If False and no matching orbitals can be found, raise a Result_unavailable_error.
        :return: A (possibly empty) Molecular_orbital_list object.
        """
        # If we've been given a generic criteria, decide what it is actually talking about.        
        if criteria is not None:
            try:                
                # If our string start with a sign (+ or -) then it's a HOMO_difference.
                if criteria[:1] == "+" or criteria[:1] == "-":
                    HOMO_difference = int(criteria)
                # If its and integer (or a string that looks like one), then its a level.
                elif criteria.isdigit() or isinstance(criteria, int):
                    level = int(criteria)
                # Otherwise we assume it a label.
                else:
                    label = criteria
                
            except Exception:
                # We couldn't parse criteria, get upset.
                raise ValueError("Unable to understand given search criteria '{}'".format(criteria))
        
        # Now get our filter func.
        if label is not None:
            filter_func = lambda mo: mo.label != label
        elif HOMO_difference is not None:
            filter_func = lambda mo: mo.HOMO_difference != HOMO_difference
        elif level is not None:
            filter_func = lambda mo: mo.level != level
        else:
            raise ValueError("Missing criteria to search by; specify one of 'criteria', 'label', 'HOMO_difference' or 'level'")
        
        # Now search.
        found = type(self)(filterfalse(filter_func, self))
        
        # Possibly panic if the list is empty.
        if not allow_empty and len(found) == 0:
            # Work out what we searched for to give better error.
            if label is not None:
                criteria_string = "label = '{}'".format(label)
            
            elif HOMO_difference is not None:
                criteria_string = "HOMO difference = '{}'".format(HOMO_difference)
            
            elif level is not None:
                criteria_string = "index = '{}'".format(level)
            
            raise Result_unavailable_error("Orbital", "could not find orbital with {}".format(criteria_string))
        
        return found
    
    def ordered(self):
        """
        Return a copy of this list of MOs that is ordered in terms of energy and removes duplicate MOs.
        """
        ordered_list = type(self)(set(self))
        #ordered_list.sort(key = lambda mo: mo.level)
        ordered_list.sort()
        return ordered_list
        
    @classmethod
    def from_parser(self, parser, cls = None):
        """
        Construct a Molecular_orbital_list object from an output file parser.
        
        :param parser: An output file parser.
        :param cls: Optional class of objects to populate this list with, should inherit from Molecular_orbital. Defaults to Molecular_orbital if only one set of orbitals are available, or Alpha_orbital if both alpha and beta are available (in which case you should call Molecular_orbital_list.from_cclib() again with cls = Beta_orbital to get beta as well).
        :returns: The new Molecular_orbital_list object. The list will be empty if no MO data is available.
        """
        try:
            # Set our default class if we've not been given one.
            if cls is None:
                # Check to see if we have only 'alpha' or beta as well.
                if len(parser.data.moenergies) == 1:
                    cls = Molecular_orbital
                else:
                    cls = Alpha_orbital
            # Get our list.
            return self(cls.list_from_parser(parser))
        except AttributeError:
            return self()
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        # Decide which class to use.
        if data['spin_type'] == "none":
            cls = Molecular_orbital
        elif data['spin_type'] == "alpha":
            cls = Alpha_orbital
        elif data['spin_type'] == "beta":
            cls = Beta_orbital
        elif data['spin_type'] is None and len(data['values']) == 0:
            # No orbitals, no type defined.
            return self()
        
        return self(cls.list_from_dump(data['values'], result_set, options))
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dict = {
            "dE(HOMO-LUMO)": {
                "value": self.HOMO_LUMO_energy if self.safe_get("HOMO_LUMO_energy") else None,
                "units": "eV"
            },
            "num_occupied": len(self.occupied),
            "num_virtual": len(self.virtual),
            "values": super()._dump_(digichem_options, all),
            "spin_type": self.safe_get("spin_type")
        }
#         # Add HOMO and LUMO
#         for label in ["HOMO", "LUMO"]:
#             try:
#                 orbital = self.find(label)
#                 dump_dict[label] = {
#                     "value": float(orbital.energy),
#                     "units": "eV"
#                 }
#                 
#             except Result_unavailable_error:
#                 pass
        
        return dump_dict
    
    def find_common_level(self, *other_lists, HOMO_difference):
        """
        Find either:
            The orbital with the lowest level that has no less than the given negative HOMO_difference.
        or:
            The orbital with the highest level that has no more than the given positive HOMO_difference.
        Across one or more orbital lists.
        
        The method is useful for determining which orbitals to traverse between two limits from the HOMO-LUMO gap.
        
        :raises Result_unavailable_error: If all of the given orbital_lists (including this one) are empty.
        :param *other_lists: Optional lists to search. If none are given, then only this orbital_list is search.
        :param HOMO_difference: The distance from the HOMO to search for. Negative values indicate HOMO-n, positive values indicate LUMO+(n-1). The LUMO should be at +1 by definition.
        :return: The level (as an integer) of the matching orbital.
        """
        # Our list of orbitals that match our criteria.
        found_orbitals = []
        
        # Our list of orbital_list objects to look through.
        orbital_lists = list(other_lists)
        orbital_lists.append(self)
        
        # Cant think of a better way to do this...
        if HOMO_difference <= 0:
            search_func = lambda orbital: orbital.HOMO_difference < HOMO_difference
        else:
            search_func = lambda orbital: orbital.HOMO_difference > HOMO_difference
        
        # Loop through each list and search.
        for orbital_list in orbital_lists:
            # Now search each list for orbitals that match our criteria.
            matching_orbitals = list(filterfalse(search_func, orbital_list))
            
            # We can just add all the orbitals we find because we'll only look at the lowest/highest anyway.
            found_orbitals.extend(matching_orbitals)
            
        # Get a list of orbital levels that match our criteria.
        orbital_levels = [orbital.level for orbital in found_orbitals]
        
        # Now either return the smallest or largest orbital level, depending on what we were asked for.
        try:
            if HOMO_difference <= 0:
                # The lowest orbital below HOMO.
                return min(orbital_levels)
            else:
                # The highest orbital above LUMO.
                return max(orbital_levels)
        except ValueError:
            # Min/Max couldn't find anything, this should only happen if all orbital_lists are completely empty.
            raise Result_unavailable_error("Common orbital level", "there are no orbitals")
    
    def assign_total_level(self, other_list):
        """
        Assign total levels to the orbitals in this list and another orbital list.
        
        'Total levels' are the index +1 of each orbital out of both alpha and beta orbitals.
        
        :param: other_list: If this list contains alpha orbitals, other_list should contain beta_orbitals (vice versa).
        """
        for index, orbital in enumerate(sorted(chain(self, other_list), key = lambda orbital: orbital.energy)):
            orbital.total_level = index +1
            
    @classmethod
    def merge_orbitals(self, molecular_orbital_lists, beta_orbital_lists):
        """
        """
        if len(molecular_orbital_lists) != len(beta_orbital_lists):
            # Panic.
            raise TypeError("molecular_orbital_lists and beta_orbital_lists must be of the same length")
        
        MOs = molecular_orbital_lists[0]
        betas = beta_orbital_lists[0]
        
        # Check all other lists are the same.
        for MO_list, beta_list in zip(molecular_orbital_lists[1:], beta_orbital_lists[1:]):
            # If this list has atoms and our current doesn't, use this as our base list.
            # We only do this if both lists are empty, so we don't have alpha and beta orbitals from two different results.
            if len(MOs) == 0 and len(betas) == 0:
                if len(MO_list) > 0:
                    MOs = MO_list
                if len(beta_list) > 0:
                    betas = beta_list
                
            else:
                MOs.check_equivalent(MOs, MO_list)
                betas.check_equivalent(betas, beta_list)
                
        return (MOs, betas)
    
    @classmethod
    def check_equivalent(self, MOs, other_MO_list):
        """
        Check that two MO lists are equivalent, issuing merge warnings if not.
        """
        for index, MO in enumerate(other_MO_list):
            try:
                other_MO = MOs[index]
                if not self.are_items_equal(MO, other_MO):
                    warnings.warn(self.MERGE_WARNING)
            except IndexError:
                warnings.warn(self.MERGE_WARNING)
        
            
    @classmethod
    def merge(self, *multiple_lists):
        """
        Merge multiple lists of MOs into a single list.
        """
        raise NotImplementedError("Molecular orbital lists do not implement the merge() method; use merge_orbitals() instead")
    
    @classmethod
    def are_items_equal(self, MO, other_MO):
        """
        A method which determines whether two items are the same for the purposes of merging.
        """
        return math.isclose(MO.energy, other_MO.energy) and MO.level == other_MO.level and MO.spin_type == other_MO.spin_type
    

class Molecular_orbital(Result_object, Floatable_mixin):
    """
    Class representing a molecular orbital.
    """
    
    # True MOs don't have a spin.
    spin_type = "none"
    
    def __init__(self, level, HOMO_difference, energy, symmetry = None, symmetry_level = None):
        """
        Constructor for MOs.
        
        :param level: The 'level' of this MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
        :param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
        :param energy: The energy of this MO (in eV).
        :param symmetry: The symmetry of this MO.
        :param symmetry_level: The ordered (by energy) index of this orbital in terms of orbitals that share the same symmetry.
        """
        super().__init__()
        self.level = level
        self.HOMO_difference = HOMO_difference
        self.symmetry = symmetry
        self.symmetry_level = symmetry_level
        self.energy = energy
                
        # This needs to be set outside of this constructor, because it relies on multiple lists of orbitals being completed.
        # Total level is the index +1 of this orbital out of both alpha and beta orbitals (for restricted this will == level).
        self.total_level = None
        
    @property
    def is_occupied(self):
        """
        Determine whether this orbital is occupied.
        """
        return self.HOMO_difference <= 0
            
    def __float__(self):
        return float(self.energy)
    
    @property
    def HOMO_level(self):
        """
        The level of the HOMO in the collection of orbitals of which this orbital is a member.
        """
        return self.level - self.HOMO_difference
    
    @property
    def LUMO_level(self):
        """
        The level of the LUMO in the collection of orbitals of which this orbital is a member.
        """
        return self.HOMO_level +1        
    
    @property
    def label(self):
        """
        A label describing this MO in terms of its proximity to the HOMO and LUMO.
        
        :return: A string label, of the form either HOMO-n or LUMO+n.
        """
        # The label we return depends on how close to the HOMO we are.
        if self.level == self.HOMO_level:
            # We are the HOMO.
            label = "HOMO"
        elif self.level < self.HOMO_level:
            # We are below the HOMO (and presumably occupied).
            label = "HOMO{}".format(self.level - self.HOMO_level)
        elif self.level == self.LUMO_level:
            # We are the LUMO.
            label = "LUMO"
        else:
            # We are above the LUMO (and presumably unoccupied).
            label = "LUMO{0:+}".format(self.level - self.LUMO_level)
            
        return label
    
    @property
    def irrep(self):
        """
        A unique description of this orbital combining symmetry and energy, as used by Turbomole.
        """
        return "{}{}".format(self.symmetry_level, self.symmetry.lower()) if self.symmetry_level is not None and self.symmetry is not None else None
    
    @property
    def sirrep(self):
        """
        A unique description of this orbital combining symmetry, energy and spin (alpha/beta), as used by Turbomole.
        """
        if self.spin_type == "none":
            return self.irrep
        else:
            if self.spin_type == "alpha":
                spin_tag = "a"
            elif self.spin_type == "beta":
                spin_tag = "b"
            else:
                spin_tag = "?"
            return "{}_{}".format(self.irrep, spin_tag)
    
    def __eq__(self, other):
        """
        Equality operator between MOs.
        """
        return self.label == other.label
    
    def __hash__(self):
        """
        Hash operator.
        """
        return hash(tuple(self.label))
    
    # The index used to access data from cclib (which always has two lists, one for alpha one for beta).
    ccdata_index = 0
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "index": self.level,
            "label": self.label,
            "homo_distance": int(self.HOMO_difference),
            "energy": {
                "value": float(self.energy),
                "units": "eV"
            },
            "symmetry": self.symmetry,
            "symmetry_index": self.symmetry_level
        }
        
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return [self(orbital_dict['index'], orbital_dict['homo_distance'], orbital_dict['energy']['value'], orbital_dict['symmetry'], orbital_dict['symmetry_index']) for orbital_dict in data]
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of Molecular_orbital objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of Molecular_orbital objects. The list will be empty if no MO is available.
        """
        try:
            # Get and zip required params.
            if hasattr(parser.data, "mosyms"):
                # We have symmetries.
                symmetry = parser.data.mosyms[self.ccdata_index]
            else:
                # No symmetry.
                symmetry = []
                
            # Don't catch this exception; if we don't have MO energies there's nothing we can do.
            energy = parser.data.moenergies[self.ccdata_index]
            
            # Check we don't have more symmetries than we do energies.
            if len(energy) < len(symmetry):
                diff = len(symmetry) - len(energy)
                warnings.warn("Parsed more MO symmetries than MO energies; excess symmetries will be discarded ('{}' symmetries, '{}' energies)".format(len(symmetry), len(energy)))
                # Dump the excess symmetry.
                symmetry = symmetry[:-diff]
            
            # Keep a track of all the symmetries we've seen.
            symmetries = {}

            orbitals = []
            for index, (symmetry, energy) in enumerate(zip_longest(symmetry, energy, fillvalue = None)):
                if symmetry is not None:
                    try:
                        # Add one to the level.
                        symmetries[symmetry] += 1
                        symm_level = symmetries[symmetry]
                    except KeyError:
                        # We haven't seen this mult before.
                        symmetries[symmetry] = 1
                        symm_level = symmetries[symmetry]
                else:
                    symm_level = None
                    
                orbitals.append(
                    self(
                        index +1,
                        index - parser.data.homos[self.ccdata_index],
                        energy,
                        symmetry,
                        symm_level
                    )
                )
                 
            return orbitals       
                
        except (AttributeError, IndexError):
            return []
    
class Unrestricted_orbital(Molecular_orbital):
    """
    Top-level class for unrestricted orbitals.
    """
    
    def __init__(self, level, HOMO_difference, energy, spin_type, symmetry = None, symm_level = None):
        """
        Constructor for MOs.
        
        :param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
        :param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
        :param symmetry: The symmetry of the MO.
        :param energy: The energy of the MO (in eV).
        :param spin_type: The spin of this spin-orbital (either alpha or beta).
        """
        # Call parent first.
        super().__init__(level, HOMO_difference, energy, symmetry, symm_level)
        self.spin_type = spin_type
    
    @property
    def label(self):
        # Get the base of the label first.
        label = super().label
        # Append our spin type.
        label = "{} ({})".format(label, self.spin_type)
        return label

class Alpha_orbital(Unrestricted_orbital):
    """
    An alpha spin orbital (these types of orbitals are only singly occupied, electrons are spin-up).
    """
    
    def __init__(self, level, HOMO_difference, energy, symmetry, symm_level):
        """
        Constructor for alpha MOs.
        
        :param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
        :param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
        :param energy: The energy of the MO (in eV).
        :param symmetry: The symmetry of the MO.
        """
        super().__init__(level, HOMO_difference, energy, "alpha", symmetry, symm_level)
        
class Beta_orbital(Unrestricted_orbital):
    """
    A beta spin orbital (these types of orbitals are only singly occupied, electrons are spin-down).
    """
    
    # Beta orbitals use the other list in cclib.
    ccdata_index = 1
    
    def __init__(self, level, HOMO_difference, energy, symmetry, symm_level):
        """
        Constructor for beta MOs.
        
        :param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
        :param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
        :param energy: The energy of the MO (in eV).
        :param symmetry: The symmetry of the MO.
        """
        super().__init__(level, HOMO_difference, energy, "beta", symmetry, symm_level)
    