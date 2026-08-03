# General imports.
import math

from digilib.exception import Result_unavailable_error
from digilib.result import Result_container
from digilib.result import Result_object
from digilib.result.base import Floatable_mixin


class SOC_list(Result_container):
    """
    An augmented list containing spin orbit coupling.
    """
    
    def find(self, criteria):
        """
        """
        states = criteria.split(",")
        
        # Check we have two states.
        if len(states) != 2:
            raise ValueError("SOC can only be found between exactly 2 states, not {} states.".format(len(states)))
        
        return self.between(*states)
    
    def between(self, state1, state2, **kwargs):
        """
        Return the spin-orbit coupling object between the two states with given symbols.
        
        :param states: Symbols identifying the two states to get SOC between (eg, "S(1)" "T(1)").
        :param default: A default argument that will be returned if the given states cannot be found.
        """
        states = set((state1, state2))
        
        # Iterate through all our SOC.
        for soc in self:
            # Check if the two states match the two we've been asked to find.
            if len(states.intersection((soc.singlet_state.state_symbol, soc.triplet_state.state_symbol))) == 2:
                # A match!
                return soc
            
        # No match.
        # Return the default if given, otherwise panic.
        if "default" in kwargs:
            return kwargs['default']
        else:
            raise Result_unavailable_error("Spin-orbit coupling", "could not find SOC between states with symbols '{}' and '{}'".format(state1, state2))
    
    @classmethod
    def from_parser(self, parser):
        """
        Get a SOC_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: A SOC_list object. The list will be empty if no SOC is available.
        """
        # Decide which SOC class to use.
        if hasattr(parser.data, "socelements"):
            soc_cls = Spin_orbit_coupling
        else:
            soc_cls = Total_spin_orbit_coupling
        
        return self(soc_cls.list_from_parser(parser))

    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self([(Spin_orbit_coupling if 'soc' in soc_dict else Total_spin_orbit_coupling).from_dump(soc_dict, result_set, options) for soc_dict in data])
    
    def sort(self, *, key = None, **kwargs):
        """
        Sort this list in place.
        """
        # If no key is given, we use a custom default key.
        if key is None:
            key = lambda SOC: (SOC.singlet_state.level, SOC.triplet_state.level)
            
        return super().sort(key = key, **kwargs)


class Spin_orbit_coupling(Result_object, Floatable_mixin):
    """
    Class that represents spin-orbit coupling between two states.
    """
    
    def __init__(self, singlet_state, triplet_state, positive_one, zero, negative_one):
        """
        Constructor for SOC objects.
        
        :param singlet_state: The singlet excited state this coupling is between.
        :param triplet_state: The triplet excited state this coupling is between.
        :param positive_one: SOC between the singlet state and triplet sub-state with quantum number +1 in cm-1.
        :param zero: SOC between the singlet state and triplet sub-state with quantum number 0 in cm-1.
        :param negative_one: SOC between the singlet state and triplet sub-state with quantum number -1 in cm-1.
        """
        self.singlet_state = singlet_state
        self.triplet_state = triplet_state
        self.positive_one = positive_one
        self.zero = zero
        self.negative_one = negative_one
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "singlet": self.singlet_state.state_symbol,
            "triplet": self.triplet_state.state_symbol,
            "soc": {
                "+1": {
                    "value": self.positive_one,
                    "units": "c m^-1",
                },
                "0": {
                    "value": self.zero,
                    "units": "c m^-1",
                },
                "-1": {
                    "value": self.negative_one,
                    "units": "c m^-1",
                },
            },
            "rss": {
                "units": "c m^-1",
                "value": self.root_sum_square,
            },
            "mixing_coefficient": self.mixing_coefficient
        }
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(result_set.energy_states.find(data['singlet']), result_set.energy_states.find(data['triplet']), data['soc']['+1']['value'], data['soc']['0']['value'], data['soc']['-1']['value'])
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a SOC object from an output file parser.
        """
        # Go through our data.
        try:
            return [self(
                singlet_state = parser.results.energy_states.get_state(state_symbol = singlet_symbol),
                triplet_state = parser.results.energy_states.get_state(state_symbol = triplet_symbol),
                positive_one = positive_soc,
                zero = zero_soc,
                negative_one = negative_soc
            ) for (singlet_symbol, triplet_symbol), (positive_soc, zero_soc, negative_soc) in zip(parser.data.socstates, parser.data.socelements)]
        except AttributeError:
            # No data.
            return []
        
    @property
    def splitting_energy(self):
        """
        The difference in energy between these two states.
        """
        return max(float(self.singlet_state), float(self.triplet_state)) - min(float(self.singlet_state), float(self.triplet_state))
    
    @property
    def mixing_coefficient(self):
        """
        The first order mixing coefficient (λ), given by Hso / dEst.
        Be aware that λ will == math.inf if dEst == 0.
        """
        try:
            return self.energy / self.splitting_energy
        except (FloatingPointError, ZeroDivisionError):
            return math.inf
    
    @property
    def root_sum_square(self):
        """
        The root sum square of this SOC.
        """
        return math.sqrt(self.positive_one **2 + self.zero **2 + self.negative_one **2)
        
    def __float__(self):
        """
        Floatify this SOC class.
        """
        return float(self.root_sum_square)
    
    @property
    def wavenumbers(self):
        """
        SOC in wavenumbers.
        """
        return self.root_sum_square
    
    @property
    def energy(self):
        """
        SOC in eV.
        """
        try:
            return self.wavenumbers_to_energy(self.root_sum_square)
        except (FloatingPointError, ZeroDivisionError):
            return 0.0
    
    def __str__(self):
        """
        Stringify this SOC class.
        """
        return "<{}|Hso|{}> (RSS, +1, 0, -1) (cm-1): {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(self.singlet_state.state_symbol, self.triplet_state.state_symbol, self.root_sum_square, self.positive_one, self.zero, self.negative_one)


class Total_spin_orbit_coupling(Spin_orbit_coupling):
    """
    A class for representing SOC when only the total is known (and not the individual elements).
    """
    
    def __init__(self, singlet_state, triplet_state, total):
        super().__init__(singlet_state, triplet_state, None, None, None)
        
        self._root_sum_square = total
        
    @property
    def root_sum_square(self):
        return self._root_sum_square
        
    def __str__(self):
        """
        Stringify this SOC class.
        """
        return "<{}|Hso|{}> (cm-1): {:10.5f}".format(self.singlet_state.state_symbol, self.triplet_state.state_symbol, self.root_sum_square)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dic = super()._dump_(digichem_options, all)
        # Remove elements.
        del(dump_dic['soc'])
        return dump_dic
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(result_set.energy_states.find(data['singlet']), result_set.energy_states.find(data['triplet']), data['rss']['value'])
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a SOC object from an output file parser.
        """
        # Go through our data.
        try:
            return [self(
                singlet_state = parser.results.energy_states.get_state(state_symbol = singlet_symbol),
                triplet_state = parser.results.energy_states.get_state(state_symbol = triplet_symbol),
                total = total
            ) for (singlet_symbol, triplet_symbol), total in zip(parser.data.socstates, parser.data.socenergies)]
        except (AttributeError, Result_unavailable_error):
            # No data.
            return []
        