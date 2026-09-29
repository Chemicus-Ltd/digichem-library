# General imports
import math

# Digichem imports
from digichem.result import Result_object
from digichem.result.angle import Angle


class Dipole_moment_ABC(Result_object):
    """
    ABC that represents a dipole moment.
    
    See inheriting classes for concrete implementations (depending on whether they are permanent or transition, and electric or magnetic).
    """
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of DPMs that are not identical; non-equivalent DPMs will be ignored"
    
    @property
    def magnitude(self):
        """
        The dipole moment total (root-sum-square of the dipole moment vector).
        
        :deprecated: 'magnitude' is a somewhat misleading name for this attribute, use 'total' instead.
        """
        return self.total
    
    @property
    def total(self):
        """
        The dipole moment total (root-sum-square of the dipole moment vector).
        """
        return math.sqrt( self.vector_coords[0] ** 2 + self.vector_coords[1] ** 2 + self.vector_coords[2] ** 2)
    
    def __init__(self, origin_coords, vector_coords, atoms = None):
        """
        Constructor for Dipole_moment objects.
        
        :param origin_coords: Tuple of (x, y, z) coordinates of the origin of this dipole moment.
        :param vector_coords: Tuple of (x, y, z) coordinates of this dipole  moment.
        :param atoms: Optional atoms atom list object which is used to calculate addition data about this dipole moment. 
        """
        super().__init__()
        # The start of our vector, normally (0,0,0).
        self.origin_coords = origin_coords
        # The end of our vector.
        self.vector_coords = vector_coords
        
        # Save our atoms object.
        self.atoms = atoms if atoms is not None else []
        
        # Save a name describing which dipole we are (permanent vs transition etc).
        self.dipole_type = "permanent"
    
    @property
    def name(self):
        """
        Name that describes this dipole moment.
        """
        return "PDM"
    
    @property
    def origin_coords(self):
        """
        The origin coords of this vector as a tuple of (x, y, z). origin_coords is automatically realigned by the atoms alignment object of this dipole, use _origin_coords if you do not want this behaviour.
        """
        if len(self.atoms) > 0:
            return self.atoms.apply_transformation(self._origin_coords)
        else:
            return self._origin_coords
        
    @origin_coords.setter
    def origin_coords(self, value):
        """
        Set the origin coords of this dipole.
        """
        self._origin_coords = value
        
    @property
    def vector_coords(self):
        """
        The ending coords of this vector as a tuple of (x, y, z). vector_coords is automatically realigned by the atoms alignment object of this dipole, use _vector_coords if you do not want this behaviour.
        """
        if len(self.atoms) > 0:
            return self.atoms.apply_transformation(self._vector_coords)
        else:
            return self._vector_coords
        
    @vector_coords.setter
    def vector_coords(self, value):
        """
        Set the ending coords of this dipole.
        """
        self._vector_coords = value
        
    def __float__(self):
        """
        Floatify this Dipole_moment.
        """
        return self.magnitude
    
    def __eq__(self, other):
        """
        Is this dipole moment equal to another?
        
        A dipole moment is considered equivalent only if all coordinates match.
        """
        return self.origin_coords == other.origin_coords and self.vector_coords == other.vector_coords
    
    @property
    def X_axis_angle(self):
        """
        The angle between this dipole moment and the X axis of the atom set of the molecule.
        """
        if self.total == 0:
            return Angle(0)
        elif len(self.atoms) > 0:
            return Angle(math.fabs(self.atoms.get_X_axis_angle(self.origin_coords, self.vector_coords)), "rad")
        else:
            return None
        
    @property
    def XY_plane_angle(self):
        """
        The angle between this dipole moment and the XY plane of the atom set of the molecule.
        """
        if self.total == 0:
            return Angle(0)
        elif len(self.atoms) > 0:
            return Angle(math.fabs(self.atoms.get_XY_plane_angle(self.origin_coords, self.vector_coords)), "rad")
        else:
            return None
    
    @property
    def gaussian_cgs_vector(self):
        """
        The vector of this dipole moment in Gaussian-cgs units.
        
        Note that the sign of the vector will be reversed in Gaussian-cgs units.
        See J. Phys. Chem. Lett. 2021, 21, 686-695.
        """
        return tuple(self.to_gaussian_cgs(coord) for coord in self.vector_coords)
    
    @property
    def gaussian_cgs(self):
        """
        The total of this dipole moment in Gaussian-cgs units.
        
        Note that the sign of the vector will be reversed in Gaussian-cgs units.
        See J. Phys. Chem. Lett. 2021, 21, 686-695.
        """
        cgs_vector = self.gaussian_cgs_vector
        return math.sqrt( cgs_vector[0] **2 + cgs_vector[1] **2 + cgs_vector[2] **2 )
    
    def angle(self, other, cgs = True):
        """
        Find the angle between this dipole moment and another.
        
        :param cgs: Whether to find the angle in Gaussian-cgs units or standard units (the direction of magnetic dipole moments is reversed in the former).
        """
        cos_angle = self.cos_angle(other, cgs)
        try:
            angle = math.acos(cos_angle)
        
        except ValueError:
            # Rounding errors can result in slightly exceeding the bounds of acos.
            # TODO: define what 'slightly' means.
            if cos_angle > 1:
                angle = 0

            elif cos_angle < -1:
                 angle = math.pi
                
            else:
                raise

        return Angle(angle)
        
    def cos_angle(self, other, cgs = True):
        """
        Find the cosine of the angle between this dipole moment and another.
        
        :param cgs: Whether to find the angle in Gaussian-cgs units or standard units (the direction of magnetic dipole moments is reversed in the former).
        """
        if cgs:
            vector_a = self.gaussian_cgs_vector
            vector_b = other.gaussian_cgs_vector
        
        else:
            vector_a = self.vector_coords
            vector_b = other.vector_coords
            
        vector_a_total = math.sqrt( vector_a[0] **2 + vector_a[1] **2 + vector_a[2] **2 )
        vector_b_total = math.sqrt( vector_b[0] **2 + vector_b[1] **2 + vector_b[2] **2 )
        
        try:
            return (vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2]) / (vector_a_total * vector_b_total)
        except (FloatingPointError, ZeroDivisionError):
            return 0
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "total": {
                "value": self.total,
                "units": self.units
            },
            "origin": {
                "x": {
                    "value": float(self.origin_coords[0]),
                    "units": "Å",
                },
                "y": {
                    "value": float(self.origin_coords[1]),
                    "units": "Å",
                },
                "z": {
                    "value": float(self.origin_coords[2]),
                    "units": "Å",
                }
            },
            "vector": {
                "x": {
                    "value": float(self.vector_coords[0]),
                    "units": self.units
                },
                "y": {
                    "value": float(self.vector_coords[1]),
                    "units": self.units
                },
                "z": {
                    "value": float(self.vector_coords[2]),
                    "units": self.units
                },
            },
            "x-angle": {
                "value": self.X_axis_angle.angle,
                "units": self.X_axis_angle.units
            },
            "xy-angle": {
                "value": self.XY_plane_angle.angle,
                "units": self.XY_plane_angle.units
            }
        }


class Electric_dipole_moment_mixin():
    """
    Mixin class for electric dipole moments.
    """
    
    @property
    def units(self):
        """
        The units of this dipole moment object.
        """
        # Debye
        return "D"
    
    @property
    def coulomb_meters(self):
        """
        The total of this dipole moment in Coulomb Meters.
        """
        return float(self) * 3.335640952e-30
    
    @classmethod
    def D_to_au(self, value):
        """
        Convert an electric dipole moment in D to a.u.
        """
        return value / 2.541746473
    
    @classmethod
    def to_gaussian_cgs(self, value):
        """
        Convert an electric dipole moment in D to Gaussian-cgs units.
        See J. Phys. Chem. Lett. 2021, 21, 686-695.
        """
        # cgs = au * e * a0
        return self.D_to_au(value) * 4.80320425e-10 * 5.29177210903e-9


class Magnetic_dipole_moment_mixin():
    """
    Mixin class for magnetic dipole moments.
    """
    
    @property
    def units(self):
        """
        The units of this dipole moment object.
        """
        return "a.u."
    
    @property
    def bohr_magnetons(self):
        """
        The total of this magnetic dipole moment in bohr magnetons.
        """
        return self.total * 2
    
    @property
    def ampere_square_meters(self):
        """
        The total of this magnetic dipole moment in SI units.
        """
        return self.total * 1.8548e-23
    
    @classmethod
    def to_gaussian_cgs(self, value):
        """
        Convert a magnetic dipole moment in a.u. to Gaussian-cgs units.
        See J. Phys. Chem. Lett. 2021, 21, 686-695.
        """
        return value * -9.2740100783e-21


class Dipole_moment(Dipole_moment_ABC, Electric_dipole_moment_mixin):
    """
    Class that represents the (permanent) dipole moment of a molecule.
    """
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self((data['origin']['x']['value'], data['origin']['y']['value'], data['origin']['z']['value']), (data['vector']['x']['value'], data['vector']['y']['value'], data['vector']['z']['value']), atoms = result_set.atoms)
    
    @classmethod
    def from_parser(self, parser):
        """
        Construct a Dipole_moment object from an output file parser.
        
        :param parser: An output file parser.
        :result: A single Dipole_moment object, or None if no dipole information is available.
        """
        try:
            return self(parser.data.moments[0], parser.data.moments[1], parser.results.atoms)
        except AttributeError:
            return None
    