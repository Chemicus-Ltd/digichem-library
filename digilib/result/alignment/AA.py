from digilib.result.alignment.FAP import Furthest_atom_pair
from statistics import  mean
import math
import statistics as stats
from digilib.result.alignment import Axis_swapper_mix

class Average_angle(Furthest_atom_pair, Axis_swapper_mix):
    """
    An enhancement to the Kebab skewer method for estimating molecule linearity, featuring a different method to determine a molecule's long axis.
    """
    
    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Average Angle", "AA", "Kebab+"]
    
    def align_axes(self):
        """
        Realign the axes of our coordinate system, so they have the following meaning:
        
        X-axis: The long axis (the kebab skewer), which we define as passing through the pair of atoms with the greatest separation in our set.
        Y-axis: The middle axis, which we define as perpendicular to the the X-axis (obviously) and passing through the atom that is furthest from the X-axis. Note that the Y-axis only has to pass through one atom; there may not be a corresponding atom on the other side (but there will be if the molecule is symmetrical about the X-axis).
        Z-axis: The short axis, defined as perpendicular to both the X and Y-axes (so we have no choice where this goes).
        
        :return: Nothing. The atoms are rearranged in place.
        """
        self.align_X()
        self.align_Y()
        self.align_Z()
    
    def get_average_angle(self, thetas):
        """
        Get the average angle from a list of angles.
        """
        average_sin = stats.mean([math.sin(angle) for angle in thetas])
        average_cos = stats.mean([math.cos(angle) for angle in thetas])
        
        raw_angle = math.atan(average_sin / average_cos)
        
        # Adjust.
        if average_sin > 0 and average_cos > 0:
            return raw_angle
        elif average_cos < 0:
            return raw_angle + math.pi
        else:
            return raw_angle + math.pi * 2
        
    
    def align_X(self):
        """
        Align the X axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        
        Unlike normal Kebab, we don't align our X axis along the line drawn between our two most separated atoms. Instead, our X axis is defined so as to pass as close to as many points as possible, more accurately passing through the molecule.
        """
        # First get a list of all our coordinates.
        coords = self.get_coordinate_list()
        
        # Set our origin so as to be in the middle of all our points.
        self.translate((-mean(coords[0]) , -mean(coords[1]), -mean(coords[2])))
        
        # Get an average angle from all our points.
        # Rotate xy coords (along Z axis).
        
        theta = self.get_average_angle(([math.atan2(atom.coords[1], atom.coords[0]) for atom in self]))
        self.rotate_XY(theta)
        
        # Rotate xz coords (along Y axis).
        # Get our angle again.
        
        theta = self.get_average_angle([math.atan2(atom.coords[2], atom.coords[0]) for atom in self])
        self.rotate_XZ(theta)
        
    def align_Y(self):
        """
        Align the Y axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        # Rotate xz coords (along Y axis).
        # Get our angle.
        theta = self.get_average_angle([math.atan2(atom.coords[2], atom.coords[1]) for atom in self])
        self.rotate_YZ(theta)
        
    def align_Z(self):
        """
        Align the Z axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        
        In kebab+, this method can reorientate the molecule, as it ensures X is the longest axis, Y is the second longest, and Z is the shortest.
        """
        self.reassign_axes()
        
#         # If Z is bigger than Y, rotate.
#         if (self.Z_length > self.Y_length):
#             self.rotate_YZ(math.pi/2)
#         
#         # If Y is bigger than X, rotate.
#         if (self.Y_length > self.X_length):
#             self.rotate_XY(math.pi/2)
#         
#         return
        
        
        
        
        
#         # Get the length of each axis.
#         axes = [(0, self.get_axis_length(0)), (1, self.get_axis_length(1)), (2, self.get_axis_length(2))]
#         # Sort them in terms of length.
#         axes.sort(key = lambda axis: axis[1], reverse = True)
#         new_axes = (axes[0][0], axes[1][0], axes[2][0])
#         # Rearrange axes.
#         self.swap_axes(new_axes)