from digilib.result.alignment.AA import Average_angle
from statistics import  mean

class Adjusted_average_angle(Average_angle):
    """
    A further enhancement to average angle, this method uses a second set of rotations to improve alignment.
    """
    
    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Adjusted Average Angle", "AAA"]    
    
    def align_axes(self):
        """
        Realign the axes of our coordinate system, so they have the following meaning:
        
        X-axis: The long axis (the kebab skewer), which we define as passing through the pair of atoms with the greatest separation in our set.
        Y-axis: The middle axis, which we define as perpendicular to the the X-axis (obviously) and passing through the atom that is furthest from the X-axis. Note that the Y-axis only has to pass through one atom; there may not be a corresponding atom on the other side (but there will be if the molecule is symmetrical about the X-axis).
        Z-axis: The short axis, defined as perpendicular to both the X and Y-axes (so we have no choice where this goes).
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # First complete a total AA alignment.
        super().align_X()
        super().align_Y()
        super().align_Z()
        
        # Now align using our methods.
        self.align_X()
        self.align_Y()
        self.align_Z()
    
    def align_X(self):
        """
        Align the X axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        # This is important because average angle always produces the same results regardless of the initial rotation of the molecule, which is not true of the rotation we are about to perform. Hence it is vital that we always start with the same coordinates regardless of initial orientation. AA does that for us.
        
        # Now determine the mean (note that this is not a true 'averaged angle' as used in the AA method, this is your bog-standard average).
        mean_angle = mean([self.get_theta(atom.coords[1], atom.coords[0]) for atom in self]) 
        self.rotate_XY(mean_angle)
        
        # Do the same along the Y axis.
        mean_angle = mean([self.get_theta(atom.coords[2], atom.coords[0]) for atom in self]) 
        self.rotate_XZ(mean_angle)
        
                
    def align_Y(self):
        """
        Align the Y axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        
        # Rotate xz coords (along Y axis).
        # Get our angle.
        theta = mean([self.get_theta(atom.coords[2], atom.coords[1]) for atom in self])
        self.rotate_YZ(theta)
        
        
