import math

import digilib.config


class Angle():
    """
    A class for representing angles (either in radians or degrees).
    """
    
#     # The default units to print angles in.
#     _default_angle_units = "deg"
#     
#     @classmethod
#     def set_default_angle_units(self, angle_units):
#         """
#         Change the default angle units used by subsequent objects.
#         """
#         if angle_units != "rad" and angle_units != "deg":
#             self._raise_angle_unit_error(angle_units)
#         else:
#             self._default_angle_units = angle_units

    @property
    def _default_angle_units(self):
        return  digilib.config.get_config()['angle_units']
            
    
    def __init__(self, angle, angle_units = None, output_units = "def"):
        """
        Constructor for angle objects.
        
        :param angle: The numerical value of the angle. Note that this will be stored internally in radians under the 'radians' attribute, and as degrees under 'degrees'.
        :param angle_units: The units of the angle, either "rad" for radians or "deg" for degrees. If angle is an Angle object, angle_units is ignored and determined automatically. Otherwise and if angle_units is None, "rad" is assumed as the default.
        :param output_units: The units to use when accessing the angle, either "def" to use the class default 'default_angle_units', "rad" for radians or "deg" for degrees. This can be changed post init by setting the 'units' property.
        """
        if isinstance(angle, Angle):
            angle_units = angle.units
            angle = float(angle)
        
        if angle_units is None:
            angle_units = "rad"
        
        # Check our units are valid and save our angle.
        if angle_units == "rad":
            self._angle = angle
        elif angle_units == "deg":
            self._angle = self.deg_to_rad(angle)
        else:
            # Get upset.
            self._raise_angle_unit_error(angle_units)
        
        # Save our output units (this is type checked for us).
        if output_units == "def":
            # Use our default (which can be set at the module level).
            self.units = self._default_angle_units
        else:
            self.units = output_units
            
    @classmethod
    def _raise_angle_unit_error(self, angle_unit):
        """
        Convenience function that raises an exception.
        """
        raise ValueError('\'{}\' is not a valid angle unit. Accepted values are "rad" or "deg"'.format(angle_unit))
            
    @property
    def radians(self):
        """
        The angle in radians.
        """
        return self._angle
    
    @property
    def degrees(self):
        """
        The angle in degrees.
        """
        return self.rad_to_deg(self._angle)
    
    @property
    def angle(self):
        if self._units == "rad":
            return self._angle
        elif self._units == "deg":
            return self.rad_to_deg(self._angle)
        else:
            self._raise_angle_unit_error(self._units)
        
    @property
    def units(self):
        """
        The units of the angle, either "deg" for degrees or "rad" for radians.
        """
        return self._units
    
    @units.setter
    def units(self, value):
        if value == "rad" or value == "deg":
            self._units = value
        else:
            self._raise_angle_unit_error(self._units)
    
    @property
    def pretty_units(self):
        """
        The units of the angle, currently either "°" for degrees or "rad" for radians.
        
        Note that the units returned by this method may change without warning so they are not suitable for type checking, use 'units' for actually checking the type of angle.
        """
        return self.units_to_pretty_units(self.units)
        
    @classmethod
    def units_to_pretty_units(self, units):
        """
        Convert a string describing angle units (either 'deg' or 'rad') to an appropriate symbol, currently either "°" for degrees or "rad" for radians.
        
        :param units: The units as a string.
        :return: An appropriate unit symbol as a string.
        """
        if units == "deg":
            return "°"
        elif units == "rad":
            return "rad"
        else:
            return ""
        
    
    def __float__(self):
        """
        Get the numerical value of this angle, according to our current 'units'.
        
        Use either the 'radians' or 'degrees' properties to use a specific unit.
        """
        return self.angle
    
    def __str__(self):
        """
        String representation of this angle, currently in the format: "angle" "units".
        """
        return "{} {}".format(self.angle, self.pretty_units)
            
            
    @classmethod
    def deg_to_rad(self, deg):
        return math.radians(deg)
    
    @classmethod
    def rad_to_deg(self, rad):
        return math.degrees(rad)
            
            
        