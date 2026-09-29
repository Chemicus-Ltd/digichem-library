import numpy


class Memory():
    """
    Class for storing memory-amounts.
    """
    
    # The units that we know about.
    UNITS = {
        'TB': 1000000000000,
        'GB': 1000000000,
        'MB': 1000000,
        'KB': 1000,
         'B': 1
        }
    
    # When outputting units, whether to separate the number and unit with a space.
    SPACE_UNIT = False
    
    def __init__(self, value = None, print_decimal = False, round_decimal = True):
        """
        """
        self.value = None
        self.auto = value
        self.print_decimal = print_decimal
        self.round_decimal = round_decimal
    
    @classmethod
    def from_units(self, value, units):
        """
        Create a memory object from a given amount of memory and a specified unit.
        
        :param value: The amount of memory.
        :param units: The units of memory.
        """
        return self("{}{}".format(value, units))
    
    @property
    def auto_units(self):
        """
        Get the amount of memory, as a tuple, using an automatically determined suffix.
        """
        # First, get an ordered list from our known units.
        ordered_units = sorted(self.UNITS.items(), key = lambda item: item[1])
        
        # Now see where our number best fits (thanks numpy).
        suffix_index = numpy.digitize((abs(self.value),), [ordered_unit[1] for ordered_unit in ordered_units])[0]
        
        # Wrap if we're out of bounds and convert to proper index.
        if suffix_index == 0:
            suffix_index += 0
        elif suffix_index == len(ordered_units):
            suffix_index -= 1
        else:
            suffix_index -= 1
        
        # Now check to see if we have a decimal part which is non zero.
        while not self.print_decimal and suffix_index > -1 and (self.value / ordered_units[suffix_index][1]) % 1 != 0:
            # Got a decimal part, use the next smallest unit.
            suffix_index -= 1
            
        if suffix_index < 0:
            # We ran out of units and couldn't remove our fraction.
            if self.round_decimal:
                # We're allowed to round our decimal away.
                suffix_index = 0
            else:
                raise ValueError("Unable to represent memory value '{}'B as non decimal".format(self.value))
        
        # Now get our value in the correct units.
        value = self.value / ordered_units[suffix_index][1]
        if not self.print_decimal:
            value = int(value)
            
        return (value, ordered_units[suffix_index][0])
    
    @property
    def auto(self):
        """
        Get the amount of memory, as a string, using an automatically determined suffix.
        """
        value, units = self.auto_units
        return "{}{}{}".format(value, " " if self.SPACE_UNIT else "", units)
            
        
    @auto.setter
    def auto(self, value):
        """
        Set the amount of memory, automatically determining the unit suffix (KB, MB, GB etc).
        """
        if isinstance(value, str):
            # First, determine the suffix (if any).
            suffix = None
            for unit, amount in  sorted(self.UNITS.items(), key = lambda item: item[1], reverse = True):
                if value.lower().endswith(unit.lower()):
                    suffix = unit
                    # This is a bit of a hack because all suffixes contain 'B' at the end.
                    break
                    
            # Now remove the suffix (if there is one), convert to float and multiply by the value of suffix (to get bytes).
            if suffix is not None:
                value = float(value[:-len(suffix)]) * self.UNITS[suffix]
        
        # Convert to int (because not sure it makes sense to represent fractions of bytes) and store.
        self.value = int(value)
    
    @property    
    def MiB(self):
        """
        The value of this memory in MiB.
        """
        return round(self.value/1048576, None)
    
    @property    
    def MB(self):
        """
        The value of this memory in MB.
        """
        return round(self.value/1000000, None)
    
    @property
    def KiB(self):
        """
        The value of this memory in KiB.
        """
        return round(self.value/1024, None)
    
    @property
    def KB(self):
        """
        The value of this memory in KB.
        """
        return round(self.value/1000, None)
    
    def __float__(self):
        """
        Floatify this memory amount (returning the number of bytes).
        """
        return float(self.value)
    
    def __int__(self):
        """
        Intify this memory amount (returning the number of bytes).
        
        In most cases, the int and float equivalents should be the same (because fractional bytes are not supported).
        """
        return int(self.value)
    
    def __str__(self):
        """
        Stringify this memory amount (returning the number of bytes) with an appropriate suffix.
        """
        return self.auto
    
    def __eq__(self, other):
        return int(self) == other
    
    @classmethod
    def is_memory(self, value):
        """Convenience method to determine whether a value is a valid memory amount."""
        try:
            self(value)
            return True
        
        except Exception:
            return False


class Turbomole_memory(Memory):
    """
    A class for representing memory quantities in formats suitable for turbomole.
    """
    
    # When outputting units, whether to separate the number and unit with a space.
    SPACE_UNIT = False