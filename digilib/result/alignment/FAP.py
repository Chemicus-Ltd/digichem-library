# Methodologies for determining molecule linearity.
import math

from statistics import pstdev
from digilib.result.alignment import Alignment

class Furthest_atom_pair(Alignment):
    """
    The 'kebab' (skewer) method for estimating molecule linearity. 
    """
        
    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Furthest Atom Pair", "FAP", "Kebab"]
        
    def align_axes(self):
        """
        Realign the axes of our coordinate system, so they have the following meaning:
        
        X-axis: The long axis (the kebab skewer), which we define as passing through the pair of atoms with the greatest separation in our set.
        Y-axis: The middle axis, which we define as perpendicular to the the X-axis (obviously) and passing through the atom that is furthest from the X-axis. Note that the Y-axis only has to pass through one atom; there may not be a corresponding atom on the other side (but there will be if the molecule is symmetrical about the X-axis).
        Z-axis: The short axis, defined as perpendicular to both the X and Y-axes (so we have no choice where this goes).
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # First align our X axis.
        self.align_X()
        
        # Next we want to orientate our secondary (Y) axis.
        self.align_Y()
        
        # This does nothing (but might not in the future).
        self.align_Z()
        
    def align_X(self):
        """
        Align the X axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        # First get our most separated atoms. This is probably inefficient.
        # This is a tuple of (distance, atom1, atom2).
        furthest = self.get_furthest_atom_pair()
                    
        # Our new origin is halfway along these two points.
        new_origin = ( (furthest[1].coords[0] + furthest[2].coords[0]) /2, (furthest[1].coords[1] + furthest[2].coords[1]) /2, (furthest[1].coords[2] + furthest[2].coords[2]) /2)
        # Translate to new origin.
        self.translate((-new_origin[0], -new_origin[1], -new_origin[2]))
        
        # Now we need to rotate, because we're going to set these two atoms as points on the x axis.
        # We rotate twice, once to set y = 0, once to set z = 0.
        
        # Rotate xy coords (along Z axis).
        # First determine angle (we should be able to use either atom as reference because we rotate about their midpoint).
        theta = self.get_theta(furthest[1].coords[1], furthest[1].coords[0])
        self.rotate_XY(theta)
        
        # Rotate xz coords (along Y axis).
        # Get our angle again.
        theta = self.get_theta(furthest[1].coords[2], furthest[1].coords[0])
        self.rotate_XZ(theta)
        
    def align_Y(self):
        """
        Align the Y axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        # Note that we're not just looking for the atom furthest away, we're looking for the greatest distance across.
        
        # Get the most separated pair of atoms in this plane.
        furthest = self.get_furthest_atom_pair_in_YZ()
            
        # Work out the angle between the two atoms.
        theta = self.get_theta(furthest[1].coords[2] - furthest[2].coords[2], furthest[1].coords[1] - furthest[2].coords[1])
        
        # Rotate so the Y axis is parallel to a line drawn between our two furthest atoms in this plane.
        self.rotate_YZ(theta)
        
    def align_Z(self):
        """
        Align the Z axis.
        
        You do not need to call this method yourself; it is called as part of align_axes().
        """
        # We don't need to align our Z axis, as it is automatically defined once we've aligned out X and Y.
        pass
    
    def get_deviation_from_X(self):
        # First get the distance of each atom from the X-axis.
        distance = [math.sqrt( (atom.coords[1])**2 + (atom.coords[2])**2 ) for atom in self]
        # Get the deviation.
        return pstdev(distance)
                
    def get_furthest_atom_pair(self):
        """
        Get the pair of atoms in our set that are separated by the greatest distance.
        
        :return: A tuple of the form (distance, atom1, atom2) where 'distance' is the straight line distance between 'atom1' and 'atom2'. The order of 'atom1' vs 'atom2' is essentially random and is irrelevant.
        """
        # This is a tuple of (distance, atom1, atom2).
        furthest = (0, None, None)
        
        for atom in self:
            for foreign_atom in self:
                # Compute distance.
                distance = (atom.distance(foreign_atom), atom, foreign_atom)
                # Check to see if we're further.
                if distance[0] >= furthest[0]:
                    # This distance is greater, so update.
                    furthest = distance
        
        # And return.
        return furthest
    
    def get_furthest_atom_from_X(self):
        """
        Get the atom that is furthest from the X-axis.
        
        :return: A tuple of the form (distance, atom) where 'distance' is the straight line distance from 'atom' to the X-axis.
        """
        furthest = (0, None)
        for atom in self:
            # Calculate distance
            distance = math.sqrt((atom.coords[1])**2 + (atom.coords[2])**2)
            # Check.
            if distance >= furthest[0]:
                # Set.
                furthest = (distance, atom)
                
        # And return
        return furthest
    
    def get_furthest_atom_pair_in_YZ(self):
        """
        Get the pair of atoms in our set that are separated by the greatest distance in the YZ plane.
        
        :return: A tuple of the form (distance, atom1, atom2) where 'distance' is the  distance between 'atom1' and 'atom2'. The order of 'atom1' vs 'atom2' is essentially random and is irrelevant.
        """
        furthest = (0, None, None)
        for atom in self:
            for foreign_atom in self:
                # Get the greatest distance in this plane.
                distance = (math.sqrt( (atom.coords[1] - foreign_atom.coords[1])**2 + (atom.coords[2] - foreign_atom.coords[2])**2 ), atom, foreign_atom)
                if distance[0] >= furthest[0]:
                    furthest = distance
        # And return.
        return furthest
        