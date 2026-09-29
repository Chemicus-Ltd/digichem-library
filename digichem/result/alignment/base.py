import math
from timeit import default_timer as timer
import datetime

from configurables.parent import Dynamic_parent

from digichem.result.atom import Atom_list


class Alignment(Atom_list, Dynamic_parent):
    """
    A class that carries out a series of transformations to realign the Cartesian axes of a set of atoms in a particular manner.
    """
    
    def __init__(self, atoms, *args, charge = None, **kwargs):
        """
        Constructor for this alignment class.
        
        :param atoms: The list of atoms to transform (can be any object that provides a 'coords' attribute which is a tuple of (x, y, z) coordinates.) Note that these coordinates will be transformed in place, so make a copy of your atom list if you want to preserve your original coordinates.
        """
        # Call our parent first.
        super().__init__(atoms, *args, charge = charge, **kwargs)
        
        # Keep track of the transformation we make so we can apply them later.
        # The translations applied to all atoms. 
        self.translations = (0, 0, 0)
        # The rotations in radians applied to all atoms. This is a list of tuples, of the form (axis, angle), where axis is 0 = X, 1 = Y, 2 = Z.
        self.rotations = []
        # And transform (if we have some atoms).
        self.duration = None
        if len(self) > 0:
            start_timer = timer()
            self.align_axes()
            end_timer = timer()
            self.duration = datetime.timedelta(seconds = end_timer - start_timer)
            
        #self.debug_print()
        #exit()
    
    @property
    def method_type(self):
        return self.CLASS_HANDLE[1]
    
    @property
    def human_method_type(self):
        return self.CLASS_HANDLE[0]

    
    def align_axes(self):
        """
        The 'main' method of this alignment class; executes the necessary transformations to align the given atoms.
        
        The inheriting subclass should provide a concrete implementation of this method.
        """
        pass
    
    @classmethod
    def translate_coords(self, coords, factor):
        """
        Translate a given set of coordinates by a constant amount.
        
        :param coords: A tuple of the coords (x, y, z) to translate.
        :param factor: A tuple of the amount to translate the coord by (x, y, z).
        :return: The translated coords.
        """
        return (coords[0] + factor[0],
                coords[1] + factor[1],
                coords[2] + factor[2])
    
    def translate(self, factor):
        """
        Translate all the atoms of this set by a constant amount.
        
        :param factor: A tuple of the amount to translate each atom's (x, y, z) coord by.
        """
        for atom in self:
            # Translate.
            atom.coords = self.translate_coords(atom.coords, factor)
                    
        # Update our param which keeps track of transformations.
        self.translations = self.translate_coords(self.translations, factor)
        
    @classmethod
    def rotate_coords_XY(self, coords, theta):
        """
        Rotate a set of coordinates in the XY plane (down the Z-axis).
        
        :param coords: The (X, Y, Z) coordinates to rotate.
        :param theta: The angle (in radians) to rotate by.
        :return: The rotated coords.
        """
        #TODO: Some of these rotation functions rotate in unexpected directions. This works for the alignment classes written so far, but is unexpected and confusing.
        # Rotate (z stays the same).
        return    ((coords[0] * math.cos(theta) + coords[1] * math.sin(theta)), 
                ((-coords[0]) * math.sin(theta) + coords[1] * math.cos(theta)),
                (coords[2]))
    
    @classmethod
    def rotate_coords_XZ(self, coords, theta):
        """
        Rotate a set of coordinates in the XZ plane (down the Y-axis).
        
        :param coords: The (X, Y, Z) coordinates to rotate.
        :param theta: The angle (in radians) to rotate by.
        :return: The rotated coords.
        """
        # Rotate (y stays the same).
        return    ((coords[0] * math.cos(theta) + coords[2] * math.sin(theta)),
                (coords[1]),
                ((-coords[0]) * math.sin(theta) + coords[2] * math.cos(theta)))
    
    @classmethod
    def rotate_coords_YZ(self, coords, theta):
        """
        Rotate a set of coordinates in the YZ plane (down the X-axis).
        
        :param coords: The (X, Y, Z) coordinates to rotate.
        :param theta: The angle (in radians) to rotate by.
        :return: The rotated coords.
        """
        # Rotate (x stays the same).
        #return ((coords[0]),
        #        (),
        #        ())
        return    ((coords[0]),
                (coords[1] * math.cos(theta) + coords[2] * math.sin(theta)),
                ((-coords[1]) * math.sin(theta) + coords[2] * math.cos(theta)))
        
    @classmethod
    def axis_to_index(self, axis):
        if axis == 'X' or axis == 0:
            return 0
        elif axis == 'Y' or axis == 1:
            return 1
        elif axis == 'Z' or axis == 2:
            return 2
        else:
            raise ValueError("Axis '{}' is out of bounds. Possible  values are 0 (X), 1 (Y) or 2 (Z)")

    @classmethod 
    def rotate_coords(self, coords, axis, theta):
        """
        Rotate a set of coordinates around an axis.
        
        :param coords: Tuple of (X, Y, Z) coords to rotate.
        :param axis: The axis to rotate around, either ('X' or 0), ('Y' or 1) or ('Z' or 2).
        :param theta: The angle (in radians) to rotate by.
        :return: The rotated coords.
        """
        # First decide which function to use.
        axis = self.axis_to_index(axis)
        if  axis == 0:
            func = self.rotate_coords_YZ
        elif axis == 1:
            func = self.rotate_coords_XZ
        elif axis == 2:
            func = self.rotate_coords_XY
        
        # Do the rotation.
        return func(coords, theta)
        
    def rotate(self, axis, theta):
        """
        Rotate all the atoms of this set around a given axis by a given angle.
        
        :param axis: The axis to rotate around, either ('X' or 0), ('Y' or 1) or ('Z' or 2). Coordinates in this axis will not be changed.
        :param theta: The angle (in radians) to rotate by.
        """
        # Rotate each atom.
        for atom in self:
            atom.coords = self.rotate_coords(atom.coords, axis, theta)
            
        # And to our track of all our rotations.
        self.rotations.append((self.axis_to_index(axis), theta))
        
            
    def rotate_XY(self, theta):
        """
        Rotate all the atoms of this set in the xy plane (down the z-axis).
        
        :param theta: The angle (in radians) to rotate by.
        """
        return self.rotate('Z', theta)
        
    def rotate_XZ(self, theta):
        """
        Rotate all the atoms of this set in the xz plane (down the y-axis).
        
        :param theta: The angle (in radians) to rotate by.
        """
        return self.rotate('Y', theta)
            
    def rotate_YZ(self, theta):
        """
        Rotate all the atoms of this set in the yz plane (down the x-axis).
        
        :param theta: The angle (in radians) to rotate by.
        """
        return self.rotate('X', theta)
    
    @classmethod
    def get_theta(self, opposite, adjacent):
        """
        Determine the angle needed (theta) to rotate a point.
        
        This function does atan(opposite / adjacent), watching out for div by 0.
        :param opposite: the coordinate of the point along the axis opposite to the angle to calculate. After rotation, the point will have this coord == 0.
        :param adjacent: The coordinate of the point along the axis adjacent to the angle to calculate. After rotation, the point will be on this axis.
        :return: The angle (in radians).
        """
        try:
            return math.atan(opposite / adjacent)
        except (FloatingPointError, ZeroDivisionError):
            # Think it's safe to assume the angle is 90 degrees in this instance.
            return math.pi /2
    
    def get_coordinate_list(self):
        """
        Get the coordinates of all the atoms of our set as a tuple of lists of the form ([X1, X2...], [Y1, Y2...], [Z1, Z2...]).
        
        This is transposed compared to the normal format, which is [(X1, Y1, Z1), (X2, Y2, Z2)...].
        :return: A tuple of lists of X, Y and Z coordinates. Each list is guaranteed to be the same length.
        """
        return list(map(list, zip(* [(atom.coords[0], atom.coords[1], atom.coords[2]) for atom in self])))
    
    def apply_transformation(self, coords):
        """
        Apply the transformations that have been used to align this atom set to a set of coordinates.
        
        This method is useful if you have objects that you would like to re-align to this atom set, but you don't want to influence the actual alignment process (eg, dipole moments).
        
        :param coords: A list-like set of 3 coordinates (X, Y, Z).
        :return: The transformed coordinates as a tuple.
        """
        # First, translate.
        coords = self.translate_coords(coords, self.translations)
        # Now, rotate.
        for axis, theta in self.rotations:
            coords = self.rotate_coords(coords, axis, theta)

        # Return as a tuple.
        return coords
    
    def debug_print(self):
        
        for atom in self:
            print("{}, {}, {}, {}".format(atom.element, atom.coords[0], atom.coords[1], atom.coords[2]))
            
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dict = super()._dump_(digichem_options, all)
        dump_dict['alignment_method'] = self.human_method_type
        dump_dict['alignment_duration'] = {
            "value": self.duration.total_seconds(),
            "units": "s"
        }
        return dump_dict
    
    
#     @classmethod
#     def merge(self, *multiple_lists):
#         """
#         Merge multiple lists of atoms into a single list.
#         
#         Note that it does not make logical sense to combine different list of atoms into one; hence the method only ensures that all given lists are the same and then returns the first given.
#         If the atom lists are not equivalent, a warning will be issued.
#         If any of the alignment methods are not the same, a warning will be issued.
#         """
#         alignment = multiple_lists[0]
#         
#         # Check all other lists are the same.
#         for atom_list in multiple_lists[1:]:
#             if type(alignment) != type(atom_list):
#                 warnings.warn("")
#                 
#         # Return the 'merged' list.
#         return alignment
    
    
class Axis_swapper_mix():
    """
    A class mixin for automatically re-assigning the axes so that X, Y & Z axes are in decreasing length.
    """
    
    def reassign_axes(self):
        """
        Automatically swap the axes of this class so that X, Y & Z axes are in decreasing length.
        """
        if self.X_length < self.Y_length and self.Y_length > self.Z_length:
            # Swap X and Y.
            self.rotate_XY(-math.pi/2)

        elif self.X_length < self.Z_length:
            # Swap X and Z.
            self.rotate_XZ(-math.pi/2)
            
        if self.Y_length < self.Z_length:
            # Swap. Y and Z.
            self.rotate_YZ(math.pi/2)

            
class Minimal(Alignment, Axis_swapper_mix):
    """
    A basic alignment class, only checks to ensure that the X, Y & Z axes are in decreasing length.
    """
    
    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Minimal", "MIN"]
    
    def align_axes(self):
        """
        The 'main' method of this alignment class; executes the necessary transformations to align the given atoms.
        """
        # All we do is swap axes.
        self.reassign_axes()
        
    
class Fixed(Alignment):
    """
    An alignment method that does nothing at all.
    """

    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Fixed", "FIX"]

    def align_axes(self):
        """
        The 'main' method of this alignment class; executes the necessary transformations to align the given atoms.
        """
        pass