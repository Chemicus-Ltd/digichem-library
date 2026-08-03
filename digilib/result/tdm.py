import itertools

from digilib.exception.base import Result_unavailable_error
from digilib.result.dipole_moment import Electric_dipole_moment_mixin, Dipole_moment_ABC, Magnetic_dipole_moment_mixin
from digilib.result.base import Result_object


class Transition_dipole_moment_ABC(Dipole_moment_ABC):
    """
    ABC for electric and magnetic transition dipole moments.
    """
    
    def __init__(self, state_level, *args, **kwargs):
        """
        Electric_transition_dipole_moment constructor.
        
        :param state: The excited state level this dipole transitions to.
        """
        super().__init__(*args, **kwargs)
        
        # An int describing the excited state we belong to. This is always available.
        self.state_level = state_level
        
        # The excited state object that we belong to. This may be None. We can't set this here because then both the Excited_state and Transition_dipole_moment constructors would depend on each other.
        self.excited_state = None
        
        # Save a name describing which dipole we are (permanent vs transition etc).
        self.dipole_type = "transition"
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(
            data['state_level'],
            (data['origin']['x']['value'], data['origin']['y']['value'], data['origin']['z']['value']),
            (data['vector']['x']['value'], data['vector']['y']['value'], data['vector']['z']['value']),
            atoms = result_set.atoms
            )
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        data = {'state_level': self.state_level}
        data.update(super()._dump_(digichem_options, all))
        return data
    
    @property
    def name(self):
        """
        Name that describes this dipole moment.
        """
        return "{} TDM".format(self.excited_state.state_symbol)
    
    def set_excited_state(self, excited_state):
        """
        Set the excited state object that we belong to.
        """
        # Not sure what we should do if the given excited_state has a different level to what we expect.
        if excited_state.level != self.state_level:
            pass
        self.excited_state = excited_state
        
        
class Electric_transition_dipole_moment(Transition_dipole_moment_ABC, Electric_dipole_moment_mixin):
    """
    Class that represents an electric transition dipole moment.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.electromagnetic_type = "electric"
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of TDMs from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of TDM objects. An empty list will be returned if no TDM data is available.
        """
        try:
            return [
                self(state_level = index +1,
                origin_coords = (0,0,0),
                vector_coords = (parser.au_to_debye(tdm[0]), parser.au_to_debye(tdm[1]), parser.au_to_debye(tdm[2])),
                atoms = parser.results.atoms) for index, tdm in enumerate(parser.data.etdips)]
        except AttributeError:
            return []


class Magnetic_transition_dipole_moment(Transition_dipole_moment_ABC, Magnetic_dipole_moment_mixin):
    """
    Class that represents a magnetic transition dipole moment.
    
    Unlike electric dipole moments, the units here are in a.u. (0.5 bohr magnetons).
    """
        
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.electromagnetic_type = "magnetic"
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of TDMs from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of TDM objects. An empty list will be returned if no TDM data is available.
        """
        try:
            return [
                self(state_level = index +1,
                origin_coords = (0,0,0),
                vector_coords = tdm,
                atoms = parser.results.atoms) for index, tdm in enumerate(parser.data.etmagdips)]
        except AttributeError:
            return []


class Transition_dipole_moment(Result_object):
    """
    A compound class that represents both the electric and magnetic components of a transition dipole moment.
    
    This class can also be used in most places where an electric TDM is expected, as it will pass references to TEDM attributes to the actual TEDM class.
    If the object has only a TMDM and no TEDM, references will instead be passed to the TMDM.
    """
    
    def __init__(self, electric = None, magnetic = None):
        """
        Constructor for TDM.
        
        Either or both of the substituent TDM can be None.
        """
        self.electric = electric
        self.magnetic = magnetic
        self.electromagnetic_type = "electromagnetic"
    
    def set_excited_state(self, excited_state):
        """
        Set the excited state object that we belong to.
        """
        if self.electric is not None:
            self.electric.set_excited_state(excited_state)
        
        if self.magnetic is not None:
            self.magnetic.set_excited_state(excited_state)
        
    @classmethod
    def from_parser(self, parser):
        return self.list_from_parser(parser)
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of Transition_dipole_moment from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of Transition_dipole_moment objects. An empty list will be returned if no TDM data is available.
        """
        electric = Electric_transition_dipole_moment.list_from_parser(parser)
        magnetic = Magnetic_transition_dipole_moment.list_from_parser(parser)
        
        return [self(electric_part, magnetic_part) for electric_part, magnetic_part in itertools.zip_longest(electric, magnetic)]
        
    def __getattr__(self, name):
        # This check prevents infinite recursion when pickling, and is probably a sensible check anyway...
        if name != "electric" and name != "magnetic":
            try:
                if self.electric is not None:
                    return getattr(self.electric, name)
                else:
                    return getattr(self.magnetic, name)
            except AttributeError as e:
                raise AttributeError(name) from e
        else:
            raise AttributeError(name)
        
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        # This constructor for tdms is a bit odd in that it does not work anything like the other from_dump() methods.
        # In the dumped format, tdms are not stored as their own list (but rather as a sub dict of the relevant excited state),
        # but in order to avoid recursive imports, a list of tdms does need to be created when excited states are loaded again.
        # This being the case, 'data' here is actually a list of excited states, not tdms...
        return [self.from_dump(excited_state['tdm'], result_set, options) for excited_state in data]
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        if data is None:
            return None
    
        # Build the individual parts.
        if data['electric'] is not None:
            electric = Electric_transition_dipole_moment.from_dump(data['electric'], result_set, options)
        
        else:
            electric = None
        
        if data['magnetic'] is not None:
            magnetic = Magnetic_transition_dipole_moment.from_dump(data['magnetic'], result_set, options)
        
        else:
            magnetic = None
        
        return self(electric = electric, magnetic = magnetic)
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "electric": self.electric.dump(digichem_options, all) if self.electric is not None else None,
            "magnetic": self.magnetic.dump(digichem_options, all) if self.magnetic is not None else None,
            "angle": {
                "value": float(self.angle().angle) if self.electric is not None and self.magnetic is not None else None,
                "units": self.angle().units if self.electric is not None and self.magnetic is not None else None,
            },
            "cos_angle": float(self.cos_angle()) if self.electric is not None and self.magnetic is not None else None,
            "dissymmetry_factor": float(self.g_value) if self.safe_get('g_value') is not None else None
        }
            
    
    def angle(self, cgs = True):
        """
        Return the angle between the electric and magnetic components of this transition dipole moment.
        """
        return self.electric.angle(self.magnetic, cgs)
    
    def cos_angle(self, cgs = True):
        """
        Return the cos of the angle between the electric and magnetic components of this transition dipole moment.
        """
        return self.electric.cos_angle(self.magnetic, cgs)
    
    @property
    def g_value(self):
        """
        Return the 'g value'; the dissymmerty factor of this transition dipole moment.
        
        This calculation is based on J. Phys. Chem. Lett. 2021, 12, 1, 686–695.
        """
        try:
            return (4 * self.electric.gaussian_cgs * self.magnetic.gaussian_cgs * self.cos_angle(True)) / (self.electric.gaussian_cgs **2 + self.magnetic.gaussian_cgs **2)
        
        except (FloatingPointError, ZeroDivisionError):
            return 0
        
        except AttributeError:
            raise Result_unavailable_error("g_value", "transition electric or transition magnetic dipole moment is not available") from None

