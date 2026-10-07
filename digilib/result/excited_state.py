# Hidden import for speed.
#import colour
# from digilib.result.spectroscopy import Spectroscopy_graph,\
#     Absorption_emission_graph

# General imports.
from itertools import filterfalse, zip_longest
import numpy
import math
import warnings
import itertools

from digilib.misc.text import andjoin
from digilib.exception.base import Result_unavailable_error
from digilib.result import Result_container
from digilib.result import Result_object
from digilib.result.base import Floatable_mixin
import digilib.log


class Excited_state_list(Result_container):
    """
    Class for representing a group of excited states.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    @property
    def singlet_triplet_energy(self):
        """
        Get ΔE st; the energy difference between the T1 and S1 excited states.
        
        The returned value is positive if S1 is higher in energy than T1 (the usual case), negative otherwise.
        
        :return: The singlet-triplet splitting energy (in eV).
        """
        # Get our two states of interest.
        S1 = self.get_state("S(1)")
        T1 = self.get_state("T(1)")
        # Calculate their difference.))
        return float(S1) - float(T1)
    
    @property
    def multiplicity_strings(self):
        """
        Get a string describing the different multiplicities of the excited states of this list.
        """
        group = self.group()
        mults = []
        for excited_states in group.values():
            mults.append(excited_states[0].multiplicity_string)
            
        return andjoin(mults)
    
    def difference(self, state1, state2):
        """
        Get the difference in energy between the states by given labels.
        
        The returned energy difference will always be positive.
        
        :param state1: State symbol of the first state.
        :param state2: State symbol of the second state.
        """
        states = (self.get_state(state_symbol = state1).energy, self.get_state(state_symbol = state2).energy)
        return max(states) - min(states)
    
    def find(self, criteria = None, *, state_symbol = None, level = None, mult_level = None):
        """
        Retrieve a particular excited state from its state symbol (S(1), T(1) etc.).
        
        For multiplicities from 1 -> 4, single capital letters are used ('S', 'D', 'T' or 'Q' respectively). For higher multiplicities, the integer representation is used ('5', '6', '7' and so on).
        The multiplicity level follows in brackets, eg ('S(1)') refers to the first (lowest energy) singlet state.
        
        See Excited_state.state_symbol() for a full description of possible symbols.
        
        :raises ValueError: If the requested symbol could not be found in this list.
        :param criteria: Automatically determine which criteria to search by.
        :param state_symbol: The symbol (a string) to retrieve.
        :param level: The level (an int or string that looks like an int) to retrieve.
        :param mult_level: The multiplicity and multiplicity level (as a tuple, eg '(1,1)' for S(1), '(3,2)' for T(2)) to retrieve.
        :return: The requested excited state object.
        """
        if criteria is not None:
            if isinstance(criteria, tuple):
                mult_level = criteria
            elif criteria.isdigit() or isinstance(criteria, int):
                level = int(criteria)
            else:
                state_symbol = criteria
        
        # Now get our search func.
        if state_symbol is not None:
            filter_func = lambda state: state.state_symbol != state_symbol
        elif level is not None:
            filter_func = lambda state: state.level != level
        elif mult_level is not None:
            filter_func = lambda state: state.multiplicity != mult_level[0] or state.multiplicity_level != mult_level[1]
        else:
            raise ValueError("Missing criteria to search by; specify one of 'criteria', 'state_symbol', 'level' or 'mult_level'")
        
        try:
            return next(filterfalse(filter_func, self))
        except Exception:
            # Determine what we searched by for better error reporting.
            if state_symbol is not None:
                criteria_string = "symbol = '{}'".format(state_symbol)
                
            elif level is not None:
                criteria_string = "level = '{}'".format(level)
                
            else:
                criteria_string = "mult_level = '{}'".format(mult_level)
                
            raise Result_unavailable_error("Excited state", "could not find excited state with {}".format(criteria_string)) from None
        
    def get_state(self, *args, **kwargs):
        warnings.warn("get_state() is deprecated use find() instead", DeprecationWarning)
        return self.find(*args, **kwargs)
    
    def group(self):
        """
        Group the excited states in this list by multiplicity.
        
        :return: A dictionary of grouped excited states. Each key will correspond to a multiplicity (1, 2, 3 etc) and each value will be an Excited_state_list of states with that multiplicity. 
        """
        # Dictionary of grouped states.
        grouped_excited_states = {}
        
        for excited_state in self:
            # Group each ES by multiplicity.
            try:
                # Try and append to our list.
                grouped_excited_states[excited_state.multiplicity].append(excited_state)
            except KeyError:
                # No list exists yet.
                grouped_excited_states[excited_state.multiplicity] = type(self)([excited_state])
        
        # Sorting hack for old version of python that didn't technically support sorted dicts (although the sort order was secretly maintained).
        # Just make a new dict that is sorted from the start:
        sorted_group = {}
        for mult in sorted(grouped_excited_states):
            sorted_group[mult] = grouped_excited_states[mult]
        
        return sorted_group
    
    def group_pairs(self):
        """
        Return all the unique pair combinations of the different multiplicities of this list.
        
        :returns: A two membered tuple, where the first element is a list of pairs of multiplicities [(1, 2), (1, 3), (2, 3)] etc, and the second is a grouped dictionary of multiplicites (see group()).
        """
        group = self.group()
        return (list(itertools.combinations(group, 2)), group)
    
    @property
    def num_singlets(self):
        """
        The number of singlet excited states in this list.
        """
        grouped_states = self.group()
        return len(grouped_states[1]) if 1 in grouped_states else 0
    
    @property
    def num_triplets(self):
        """
        The number of triplet excited states in this list.
        """
        grouped_states = self.group()
        return len(grouped_states[3]) if 3 in grouped_states else 0
        
    @classmethod
    def from_parser(self, parser):
        """
        Create an Excited_state_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: The populated Excited_state_list object.
        """
        return self(Excited_state.list_from_parser(parser))
    
    def assign_levels(self):
        """
        (Re)assign total and multiplicity levels of the excited states in this list.
        
        The contents of this list will be modified in place.
        """
        # A dictionary of multiplicities that we've seen so far.
        multiplicities = {}
        total_level = 0
        
        for state in self:
            # (try) and add this state's multiplicity to our dict (if we get a key error that just means it's the first time we've seen this mult).
            try:
                multiplicities[state.multiplicity] += 1
            except KeyError:
                # First time we've seen this mult, set to 1.
                multiplicities[state.multiplicity] = 1
            
            # Increment total level.
            total_level += 1
                
            # Set mult and total level.
            state.multiplicity_level = multiplicities[state.multiplicity]
            state.level = total_level

    @classmethod
    def merge(self, *multiple_lists, **kwargs):
        """
        Merge multiple lists of of the same type into a single, ordered list.
         
        Note that this method will return a new list object, but will modify the contained objects (eg, object.level) in place.
        Inheriting classes may safely override this method.
        """
        # A dictionary where each key is a multiplcity and each value is a list of excited states with that mult.
        multiplicities = {}
        for excited_states in multiple_lists:
            group = excited_states.group()
            for multiplicity in group:
                if multiplicity not in multiplicities:
                    # Firt time we've seen this mult, add.
                    multiplicities[multiplicity] = group[multiplicity]
                
                else:
                    # Already got this type of excited states.
                    warnings.warn("Attempting to merge excited states of multiplicity '{}' but this multiplicity has already been provided by an earlier result, ignoring".format(multiplicity))
        
        merged = self(itertools.chain(*multiplicities.values()), **kwargs)
        merged.sort()
        merged.assign_levels()
        return merged
    
    def _get_dump_(self):
        """
        Method used to get a dictionary used to generate on-demand values for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
        """
        return {"spectrum": self.generate_spectrum}
    
    def generate_spectrum(self, digichem_options):
        """
        Abstract function that is called to generate an on-demand value for dumping.
        
        This functionality is useful for hiding expensive properties from the normal dump process, while still exposing them when specifically requested.
        
        :param key: The key being requested.
        :param digichem_options: Digichem options being used to dump.
        """
        from digilib.result.spectroscopy import Spectroscopy_graph, Absorption_emission_graph
        
        
        # Get spectrum data.
        # For excited states, we dump two spectra. One in nm, one in eV (which have different scaling).
        # TODO: It's weird that these spectra are only available in dumped format, there should be some property/function on the class that also returns them...
        spectrum_nm = Absorption_emission_graph.from_excited_states(
            self,
            fwhm = digichem_options['absorption_spectrum']['fwhm'],
            resolution = digichem_options['absorption_spectrum']['gaussian_resolution'],
            cutoff = digichem_options['absorption_spectrum']['gaussian_cutoff'],
            filter = digichem_options['absorption_spectrum']['y_filter'],
            use_jacobian = digichem_options['absorption_spectrum']['use_jacobian']
        )
        
        
        spectrum_ev = Spectroscopy_graph(
            [(excited_state.energy, excited_state.oscillator_strength) for excited_state in self],
            fwhm = digichem_options['absorption_spectrum']['fwhm'],
            resolution = digichem_options['absorption_spectrum']['gaussian_resolution'],
            cutoff = digichem_options['absorption_spectrum']['gaussian_cutoff'],
            filter = digichem_options['absorption_spectrum']['y_filter']
        )
        
        try:
            spectrum_nm_data = spectrum_nm.plot_cumulative_gaussian()
            y_units = "arb. unit" if digichem_options['absorption_spectrum']['use_jacobian'] else "oscillator_strength"
            
            spectrum_nm_data = [{"x":{"value": float(x), "units": "nm"}, "y": {"value":float(y), "units": y_units}} for x,y in spectrum_nm_data]
            spectrum_nm_peaks = [{"x":{"value": float(x), "units": "nm"}, "y": {"value":float(y), "units": y_units}} for x, y in spectrum_nm.peaks()]
        
        except Exception:
            spectrum_nm_data = []
            spectrum_nm_peaks = []
        
        try:
            spectrum_ev_data = spectrum_ev.plot_cumulative_gaussian()
            
            spectrum_ev_data = [{"x":{"value": float(x), "units": "eV"}, "y": {"value":float(y), "units": "oscillator_strength"}} for x,y in spectrum_ev_data]
            spectrum_ev_peaks = [{"x":{"value": float(x), "units": "eV"}, "y": {"value":float(y), "units": "oscillator_strength"}} for x, y in spectrum_ev.peaks()]
        
        except Exception:
            spectrum_ev_data = []
            spectrum_ev_peaks = []
            
        return {
            "nm": {
                "values": spectrum_nm_data,
                "peaks": spectrum_nm_peaks
            },
            "ev": {
                "values": spectrum_ev_data,
                "peaks": spectrum_ev_peaks
            }    
        }
    
    def _dump_(self, digichem_options, all):
        dump_dict = {
            "values": super()._dump_(digichem_options, all),
        }
        
        # Add extra properties.
        mult_pairs, mults = self.group_pairs()
        
        # Number of states of each multiplicity.
        for mult, states in mults.items():
            dump_dict['num_{}'.format(Energy_state.multiplicity_number_to_string(mult))] = len(states)
            
        # DeST and other values.
        for mult_pair in mult_pairs:
            state1 = mults[mult_pair[0]][0]
            state2 = mults[mult_pair[1]][0]
            dump_dict['dE({}{})'.format(state1.multiplicity_symbol, state2.multiplicity_symbol)] = {
                "value": float(state1.energy - state2.energy),
                "units": "eV"
            }
            
        return dump_dict
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(Excited_state.list_from_dump(data['values'], result_set, options))


class Excited_state_transition(Result_object):
    """
    Class that represents a transition that contributes to a particular excited state.
    """
    
    def __init__(self, level, starting_mo, ending_mo, coefficient):
        """
        Constructor for excited state transitions.
        
        :param level: The 'level' of this transition. The most significant (highest probability) transition has level 1.
        :param starting_mo: The Molecular_orbital object from which this transition begins.
        :param ending_mo: The Molecular_orbital object to which the transition ends.
        :param coefficient: The coefficient of this orbital. Square to get the probability.
        """
        self.level = level
        self.starting_mo = starting_mo
        self.ending_mo = ending_mo
        self.coefficient = coefficient
        
    @property
    def probability(self):
        """
        The probability of this transition.
        """
        return self.coefficient **2
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "index": self.level,
            "start": {
                "spin": self.starting_mo.spin_type,
                "index": self.starting_mo.level,
                "label": self.starting_mo.label
            },
            "end": {
                "spin": self.ending_mo.spin_type,
                "index": self.ending_mo.level,
                "label": self.ending_mo.label
            },
            "coefficient": float(self.coefficient),
            "probability": float(self.probability)
        }
        
    @classmethod
    def list_from_parser(self, parser, threshold = 1e-3):
        """
        Create a list of excited state transitions from an output file parser.
        
        :param parser: An output file parser.
        :param threshold: The probability threshold below which transitions will be discarded.
        :return: A list of Excited_state_transition objects.
        """
        try:
            # Create a tuple of our MOs (helps us later).
            MOs = (parser.results.orbitals, parser.results.beta_orbitals)
            
            transitions_list = []
            
            # etsecs is a list of list, the outer list corresponds to each excited state, the inner list contains actual transitions.
            for excited_state_transitions in parser.data.etsecs:
    
                # We'll first create an intermediate list of keyword dicts which we'll then sort.
                data_list = []

                for (starting_mo_index, starting_mo_AB), (ending_mo_index, ending_mo_AB), coefficient in excited_state_transitions:
                    try:
                        if coefficient **2 < threshold:
                            # Skip tiny contributions.
                            continue

                        data_list.append({
                            'starting_mo': MOs[starting_mo_AB][starting_mo_index],
                            'ending_mo': MOs[ending_mo_AB][ending_mo_index],
                            'coefficient': coefficient
                        })
                    
                    except IndexError:
                        # This is fairly common in Orca 6, where only a subset of virtual orbitals are printed by default.
                        digilib.log.get_logger() \
                            .warning("Unable to construct excited state transition; transition is to/from an orbital that is not available ({} and {})".format(starting_mo_index, ending_mo_index))
                
                # Sort by probability/coefficient.
                data_list.sort(key=lambda keywords: math.fabs(keywords['coefficient']), reverse=True)
                
                # Now get a list objects and append to our big list.
                transitions_list.append(
                    [self(index+1, keywords['starting_mo'], keywords['ending_mo'], keywords['coefficient']) for index, keywords in enumerate(data_list)]
                )
                
            # All done.
            return transitions_list
        
        except AttributeError:
            # No data.
            return []
        
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        trans = []
        for tran_dict in data:
            # Get our start and end MOs.
            orbitals = {}
            for key in ("start", "end"):
                if tran_dict[key]['spin'] in ("none", "alpha"):
                    orbitals[key] = result_set.orbitals[tran_dict[key]['index'] -1]
                
                else:
                    orbitals[key] = result_set.beta_orbitals[tran_dict[key]['index'] -1]
            
            trans.append(self(tran_dict['index'], orbitals['start'], orbitals['end'], tran_dict['coefficient']))
        
        return trans


class Energy_state(Result_object, Floatable_mixin):
    """
    Class for representing different energy states of the same system.
    """
    
    # Colour categories.
    colors = [
        {"max": 400, "name": "Ultraviolet"},
        {"max": 420, "name": "Violet"},
        {"max": 470, "name": "Blue"},
        {"max": 505, "name": "Cyan"},
        {"max": 555, "name": "Green"},
        {"max": 595, "name": "Yellow"},
        {"max": 625, "name": "Orange"},
        {"max": 740, "name": "Red"},
        {"max": float("inf"), "name": "Infrared"}    
    ]
    
    def __init__(self, level, multiplicity, multiplicity_level, energy):
        """
        Constructor for Energy_state objects.
        
        :param level: The ordered (by energy) index of this energy state, where 0 is the GS and 1 is the lowest excited state.
        :param multiplicity_level: the ordered (by energy) index of this energy state in terms of states that share the same multiplicity.
        :param multiplicity: The multiplicity of this state as a number (eg, 1 for singlet, 2 for doublet). Fractional multiplicities are also accepted.
        :param energy: The energy of this state in eV. Whether this value is absolute or relative to another state depends on the implementing class.
        """
        self.level = level
        self.multiplicity = round(multiplicity) if multiplicity is not None else None
        # 'True' multiplicity is unrounded (do something smarter)
        self.true_multiplicity = multiplicity
        self.multiplicity_level = multiplicity_level
        self.energy = energy
        
    @property
    def multiplicity(self):
        """
        """
        return self._multiplicity if self._multiplicity is not None else 0
    
    @multiplicity.setter
    def multiplicity(self, value):
        """
        """
        self._multiplicity = value
    
    def __float__(self):
        """
        Float of this class.
        """
        return float(self.energy)
    
    @classmethod
    def multiplicity_number_to_string(self, multiplicity):
        if multiplicity == 1:
            return "singlet"
        elif multiplicity == 2:
            return "doublet"
        elif multiplicity == 3:
            return "triplet"
        elif multiplicity == 4:
            return"quartet"
        elif multiplicity is None or multiplicity == 0:
            return "no multiplicity"
        elif multiplicity % 1 == 0:
            # Multiplicity is an integer, so return as a stringy whole number.
            return str(int(multiplicity))
        else:
            return str(multiplicity)
        
    
    @property
    def multiplicity_string(self):
        return self.multiplicity_number_to_string(self.multiplicity)
    
    @property
    def multiplicity_symbol(self):
        multiplicity = self.multiplicity
        
        # Get a shorthand symbol if we can.
        if  multiplicity == 1:
            return "S"
        elif  multiplicity == 2:
            return "D"
        elif  multiplicity == 3:
            return "T"
        elif  multiplicity == 4:
            return "Q"
        elif multiplicity == 0:
            return "?"
        elif multiplicity % 1 == 0:
            # Multiplicity is an integer, so return as a stringy whole number.
            return str(int(multiplicity))
        else:
            return str(multiplicity)
    
    @property
    def state_symbol(self):
        """
        A short hand notation to identify this excited state.
        
        If the multiplicity is well defined (singlet, doublet, triplet etc), the symbol starts with an appropriate letter (S, D, T etc), otherwise a numeric multiplicity is used. The symbol ends with an integer in brackets, indicating the excited state's level.
        eg, S(1) is the first singlet excited state, S(2) is the second and so on.
        """
        return "{}({})".format(self.multiplicity_symbol, self.multiplicity_level)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "index": self.level,
            "symbol": self.state_symbol,
            "multiplicity": self.multiplicity,
            "multiplicity_index": self.multiplicity_level,
            "energy": {
                "value": float(self.energy),
                "units": "eV"
            }
        }


class Excited_state(Energy_state):
    """
    Class for representing an excited state.
    """
    
    def __init__(self, level, multiplicity, multiplicity_level, symmetry, energy, oscillator_strength, transitions, transition_dipole_moment = None):
        """
        Constructor for excited state objects.
        
        :param level: The 'level' of this excited state (essentially an index), where the lowest state has a level of 1, increasing by 1 for each higher state.
        :param multiplicity: The multiplicity of this excited state (as an integer or float).
        :param multiplicity_level: The 'level' of this state within the excited states that have the same multiplicity.
        :param symmetry: The symmetry of this excited state; a complicated term that contains our multiplicity among other things.
        :param energy: The energy of this excited state (in eV).
        :param oscillator_strength: The oscillator strength of this transition.
        :param transitions: The singly excited transitions which make up this transition.
        :param transition_dipole_moment: Optional transition dipole moment of the excited state.
        """
        super().__init__(level, multiplicity, multiplicity_level, energy)
        self.symmetry = symmetry
        # Oscillator strength is a dimensionless float that describes the probability of the transition from the reference state (normally ground) to this excited state.
        self.oscillator_strength = oscillator_strength
        # The transitions which contribute to this state.
        self.transitions = transitions
        
        # If we were given a TDM, set its excited state to ourself.
        if transition_dipole_moment is not None:
            transition_dipole_moment.set_excited_state(self)
                
        self.transition_dipole_moment = transition_dipole_moment
        
    @classmethod
    def get_multiplicity_from_symmetry(self, symmetry_string):
        """
        Get the multiplicity of of an excited state from its symmetry.
        
        See multiplicity() for a more detailed description of the format of multiplicity strings.
        :param symmetry_string: The symmetry string.
        :return: The multiplicity as a number (possibly an int, possibly a float).
        """
        # Split the symmetry string on the dash (-) character.
        multiplicity_string = (symmetry_string.split('-', 1))[0]
        # Convert to a number if we have a string.
        if multiplicity_string == "Singlet":
            return 1
        elif multiplicity_string == "Doublet":
            return 2
        elif multiplicity_string == "Triplet":
            return 3
        elif multiplicity_string == "Quartet":
            return 4
        elif multiplicity_string == "???":
            return None
        else:
            # Try and cast to float.
            # Float multiplicities don't make sense in the real world, but are possible results from calculations.
            return float(multiplicity_string)
    
    @property
    def wavelength(self):
        """
        The wavelength that corresponds to the energy of this excited state (in nm).
        
        :raises Result_unavailable_error: If the energy of this state is 0.
        """
        try:
        # λ = (c * h) / e
            return self.energy_to_wavelength(self.energy)
        except FloatingPointError:
            # Our energy is zero.
            raise Result_unavailable_error('excited state wavelength', "excited state '{}' energy is 0 eV".format(self.state_symbol))
        
    @property
    def joules(self):
        """
        The energy of this excited state in Joules.
        """
        return self.energy * 1.602176634e-19
        
    @property
    def color(self):
        """
        The 'color' that corresponds to the energy of this excited state (as a string).
        """
        # Go through our color presets and see which one we match.
        for color in self.colors:
            if self.wavelength < color["max"]:
                # Good match.
                return color["name"]
        return "???"
    
    @property
    def CIE_XYZ(self):
        """
        The CIE XYZ tristimulus values of the 'color' that corresponds to the energy of this excited state (as a numpy array).
        """
        import colour
        try:
            return colour.wavelength_to_XYZ(self.wavelength)
        except ValueError:
            # Wavelength is out of our colour range, so we can't see it.
            return numpy.zeros(3)
        
    @property
    def CIE_xy(self):
        """
        The CIE xy chromaticity coordinates of the 'color' that corresponds to the energy of this excited state.
        """
        import colour
        try:
            return colour.XYZ_to_xy(self.CIE_XYZ)
        except Exception:
            return numpy.zeros(2)
        
                
    @property
    def rgb(self):
        """
        The RGB values of the 'color' that corresponds to the energy of this excited state (as a list of [r, g, b] from 0 -> 255).
        """
        return self.xyz_to_rgb(self.CIE_XYZ)
        
    @classmethod
    def xyz_to_rgb(self, XYZ):
        import colour
        rgb =  [numpy.clip(clr, 0, math.inf) for clr in colour.XYZ_to_sRGB(XYZ)]
        
        # Now we normalise if one of our values exceeds 1.
        if max(rgb) > 1:
            rgb = [clr / max(rgb) for clr in rgb]
            
        # Now convert to 0 -> 255 and return.
        return [int(clr * 255) for clr in rgb]    
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dict = super()._dump_(digichem_options, all)
        dump_dict.update({
            "wavelength": {
                "value": float(self.wavelength),
                "units": "nm"
            },
            "colour": self.color,
            "cie": {
                "x": float(self.CIE_xy[0]),
                "y": float(self.CIE_xy[1]),
            },
            "symmetry": self.symmetry,
            "oscillator_strength": float(self.oscillator_strength) if self.oscillator_strength is not None else None,
            "tdm": self.transition_dipole_moment.dump(digichem_options, all) if self.transition_dipole_moment is not None else None,
            "transitions": [tran.dump(digichem_options, all) for tran in self.transitions],
        })
        return dump_dict
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of Excited_state objects from an output file parser.
                
        :param parser: An output file parser.
        :return: A list of Excited_state objects in the same order as given in cclib.
        """
        # List of excited states.
        excited_states = []
        
        # Assemble cclib's various arrays into a single list.
        # If we're missing etoscs, that's ok.
        etoscsc = getattr(parser.data, "etoscs", [])
        
        try:
            excited_states_data = list(zip_longest(parser.data.etsyms, parser.data.etenergies, etoscsc, fillvalue = 0.0))
        
        except AttributeError:
            if not hasattr(parser.data, "etsyms") and not hasattr(parser.data, "etenergies"):
                return []
            
            elif not hasattr(parser.data, "etsyms"):
                digilib.log.get_logger().warning("Unable to fully parse excited states because only energies are available, missing symmetries")
                return []
            
            elif not hasattr(parser.data, "etenergies"):
                digilib.log.get_logger().warning("Unable to fully parse excited states because only symmetries are available, missing energies")
                return []
            
            else:
                raise
        
        # First get our transitions.
        transitions = Excited_state_transition.list_from_parser(parser)
        
        # Loop through our data.
        for index, (symmetry, energy, oscillator_strength) in enumerate(excited_states_data):                    
            # Relevant transition dipole moment.
            try:
                tdm = parser.results.transition_dipole_moments[index]
            except IndexError:
                tdm = None
                
            # Get and append our object.
            # We'll set level and mult level once we've got all our objects, so just use None for now.
            excited_states.append(
                self(
                    level = None,
                    multiplicity = Excited_state.get_multiplicity_from_symmetry(symmetry),
                    multiplicity_level = None,
                    symmetry = symmetry,
                    energy = self.wavenumbers_to_energy(energy),
                    oscillator_strength = oscillator_strength,
                    transitions = transitions[index] if len(transitions) != 0 else [],
                    transition_dipole_moment = tdm
                )
            )
        
        # Now assign total and multiplicity levels which we skipped earlier.
        Excited_state_list.assign_levels(excited_states)
        
        # All done, return our list.
        return excited_states
    
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        states = []
        for es_dict in data:
            # Get our tdm (if available).
            if es_dict['tdm'] is not None:
                tdm = result_set.transition_dipole_moments[es_dict['index']-1]
            
            else:
                tdm = None
            
            states.append(self(
                es_dict['index'],
                es_dict['multiplicity'],
                es_dict['multiplicity_index'],
                es_dict['symmetry'],
                es_dict['energy']['value'],
                es_dict['oscillator_strength'],
                Excited_state_transition.list_from_dump(es_dict['transitions'], result_set, options),
                tdm
            ))
        
        return states
    