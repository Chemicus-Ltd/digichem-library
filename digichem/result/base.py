import scipy.constants
import warnings
import itertools
import math

from digichem.exception.base import Result_unavailable_error


class Result_object():
    """
    Top level class for objects that are designed to hold calculation results.
    """
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge lists of results that are not identical; non-equivalent results will be ignored"
    
    def safe_get(self, *attr_names, default = None):
        """
        Access an attribute of this object, returning default if the attribute is None or is not available (raises Result_unavailable_error).
        """
        # DO YOU LIKE MY B-TEAM SPAGHETTI CODE!?
        
        # Get the first name, which we'll be getattr'ing.
        first_name = attr_names[0]
        remaining_names = attr_names[1:]
        
        # Get.
        try:
            attr = getattr(self, first_name)
            
            if attr is not None and len(remaining_names) > 0:
                if isinstance(attr, Result_object):
                    # We have remaining names to resolve, call attr's safe_get().
                    attr = attr.safe_get(*remaining_names, default = default)
                else:
                    # We have remaining names to resolve, but our attr is not a Result_set (so we can't use safe_get()).
                    for remaining_name in remaining_names:
                        attr = getattr(attr, remaining_name)
        except Result_unavailable_error:
            return default
        
        # And done.
        return attr
    
    @classmethod
    def wavelength_to_energy(self, wavelength):
        """
        Convert a wavelength (in nm) to energy (in eV).
        """
        # e = (c * h) / λ
        return ((scipy.constants.speed_of_light * scipy.constants.Planck) / (wavelength / 1000000000)) / scipy.constants.eV
    
    @classmethod
    def energy_to_wavelength(self, energy):
        """
        Convert an energy (in eV) to wavelength (in nm).
        """
        # λ = (c * h) / e
        return ((scipy.constants.speed_of_light * scipy.constants.Planck) / (energy * scipy.constants.eV)) * 1000000000
    
    @classmethod
    def energy_to_wavenumbers(self, energy):
        """
        Convert an energy (in eV) to wavenumbers (cm-1)
        """
        return 10000000 / self.energy_to_wavelength(energy)
    
    @classmethod
    def wavenumbers_to_energy(self, wavenumbers):
        """
        Convert wavenumbers (in cm-1) to energy (in eV).
        """
        return self.wavelength_to_energy((1 / wavenumbers) * 10000000)
    
    @classmethod
    def merged_attr(self, name, objects):
        """
        Merge a number of str like attributes.
        
        :param name: The name of the attribute to merge.
        :param objects: A list of objects to merge from.
        :return: A single string containing the merged attributes.
        """
        attributes = list(dict.fromkeys([getattr(obj, name) for obj in objects]))
        attributes = [attribute for attribute in attributes if attribute is not None]
        return ", ".join(attributes) if len(attributes) > 0 else None
    
    
    @classmethod
    def merge(self, *multiple_objects):
        """
        Merge multiple Result objects of the same type into a single Result object.
        
        Note that this default implementation is intended for objects that cannot actually be merged;
        truly mergable objects should implement their own merge() method.
        """
        # Check all items are the same.
        if not all(obj == multiple_objects[0] for obj in multiple_objects if obj is not None):
            warnings.warn(self.MERGE_WARNING)
            
        return multiple_objects[0]
    
    def calculate(self, item, digichem_options):
        """
        Retrieve/calculate an on-demand value, caching the value if it has not already been retrieved.

        This function is a wrapper around _get_dump_(), caching the calculation result to avoid
        expensive recalculation.

        :param item: The key corresponding to an item in this object's _get_dump_() dict.
        :param digichem_options: Digichem options that will be passed to _get_dump_() (unused on a cache hit).
        """
        if not hasattr(self, "_dump_cache"):
            self._dump_cache = {}

        if item in self._dump_cache:
            # Cache hit.
            return self._dump_cache[item]

        else:
            # Cache miss.
            # Generate (and cache) the data.
            self._dump_cache[item] = self.generate_for_dump()[item](digichem_options)

            # And return of course.
            return self._dump_cache[item]
    
    def generate_for_dump(self):
        """
        Method used to get a dictionary used to generate on-demand values for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
        """
        warnings.warn("generate_for_dump() is deprecated, use calculate() or _get_dump_() instead", DeprecationWarning)
        return self._get_dump_()
    
    def dump(self, digichem_options, all = False):
        # First, get simple key:value pairs.
        base_dict = self._dump_(digichem_options, all)

        # Then extend with generator ones if we have any.
        if all:
            generators = self._get_dump_()
            for key in generators:
                base_dict[key] = generators[key](digichem_options)
        
        return base_dict
    
    def _dump_(self, digichem_options, all):
        """
        Abstract function that is called to dump the value of the result object to a primitive type, suitable for serializing with yaml.
        
        For real objects, this will normally be a dict with (at least) the following items:
         - 'value': The value of the result (for example, 34.5.
         - 'units': The units of the result (for example, k m^-s).
        """
        raise NotImplementedError("Implement in subclass")
    
    def _get_dump_(self):
        """
        Method used to get a dictionary used to generate on-demand values for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
        """
        return {}


class Floatable_mixin():
    """
    A mixin class for result objects that have an unambigious numerical representation.
    
    This mixin class expects the implementing class to define __float__().
    """
    
    def __lt__(self, other):
        """
        Is this class less than another?
        """
        return float(self) < float(other)
    
    def __le__(self, other):
        """
        Is this class equal or less than another?
        """
        return math.isclose(float(self), float(other)) or float(self) < float(other)
    
    def __eq__(self, other):
        """
        Is this class equal to another?
        """
        #return float(self) == float(other)
        return math.isclose(float(self), float(other))
    
    def __ge__(self, other):
        """
        Is this class greater or equal to another?
        """
        return math.isclose(float(self), float(other)) or float(self) > float(other)
    
    def __gt__(self, other):
        """
        Is this class greater than another?
        """
        return float(self) > float(other)
    

class Unmergeable_container_mixin():
    """
    Mixin for Result_containers that cannot be merged.
    
    Classes which inherit from this mixin should also set a class attribute with the name MERGE_WARNING, which should be a message to be displayed when attempting to merge this container with another.
    """
    
    @classmethod
    def false_merge(self, *multiple_lists, **kwargs):
        """
        'Merge' a number of containers of the same type.
        
        As Unmergeable_container objects cannot be merged, this method will instead check that all given containers are equivalent, issuing warnings if this is not the case.
        
        :param multiple_lists: A number of container objects to 'merge'.
        :param kwargs: Key-word arguments to pass to the returned container.
        """
        items = self(multiple_lists[0], **kwargs)
        
        # Check all other lists are the same.
        for item_list in multiple_lists[1:]:
            # If this list has atoms and our current doesn't, use this as our base list.
            if len(item_list) > 0 and len(items) == 0:
                items = item_list
            else:
                for index, item in enumerate(item_list):
                    try:
                        other_item = items[index]
                        if not self.are_items_equal(item, other_item):
                            warnings.warn(self.MERGE_WARNING)
                    except IndexError:
                        warnings.warn(self.MERGE_WARNING)
                
        # Return the 'merged' list.
        return items
    
    @classmethod
    def are_items_equal(self, item, other_item):
        """
        A method which determines whether two items are the same for the purposes of merging.
        
        This default implementation simply checks for equivalence between the two items.
        """
        return item == other_item
        
class Result_container(list, Result_object):
    """
    Top level class for Result_objects that hold a list of results.
    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        Result_object.__init__(self)
    
    def assign_levels(self):
        """
        (Re)assign total levels of the objects in this list.
        
        The contents of this list will be modified in place.
        """
        total_level = 0
        
        for state in self:            
            # Increment total level.
            total_level += 1
                
            # Set level
            state.level = total_level
            
    @classmethod
    def merge_default(self, *multiple_lists, **kwargs):
        """
        The default implementation of the merge method.
        """
        # First, get a new list containing all objects.
        merged_list = self(itertools.chain(*multiple_lists), **kwargs)
        
        # Now sort.
        # This works because the contained objects are all Floatables (see Floatable_mixin) and therefore have well-defined sort mechanics,
        # or otherwise define their own __eq__ and related methods.
        merged_list.sort()
        
        # Next, we need to reassign levels of the contained objects.
        merged_list.assign_levels()
        
        # All done.
        return merged_list
            
    @classmethod
    def merge(self, *multiple_lists, **kwargs):
        """
        Merge multiple lists of of the same type into a single, ordered list.
        
        Note that this method will return a new list object, but will modify the contained objects (eg, object.level) in place.
        Inheriting classes may safely override this method.
        """
        if hasattr(self, "false_merge"):
            return self.false_merge(*multiple_lists, **kwargs)
        else:
            return self.merge_default(*multiple_lists, **kwargs)
        
    def _dump_(self, digichem_options, all):
        return [item.dump(digichem_options, all) for item in self]
    