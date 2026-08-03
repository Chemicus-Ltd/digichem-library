import math
from scipy.constants import epsilon_0, Planck, c

from digilib.exception.base import Digichem_exception, Result_unavailable_error
from digilib.result.excited_state import Excited_state
from digilib.result import Result_object

# Dirac constant.
h_bar = Planck / (math.pi *2)


class Emissions(Result_object):
    """
    Class that holds all emissions from a result.
    """
    
    def __init__(self, adiabatic = None, vertical = None):
        """
        """
        # For now, we have no way of knowing which state is being optimised.
        # As such, we assume it's the lowest of each mult (because this is most
        # common thanks to kasha's rule). This being the case, each emission is
        # stored in a dict, which each key is the multiplicity. Each value is not
        # a list (because we only have one emission for each mult), but the 
        # emission object itself.
        # TODO: Improve once we can detect which state is being optimised.
        self.adiabatic = adiabatic if adiabatic is not None else {}
        self.vertical = vertical if vertical is not None else {}
        
    def _dump_(self, digichem_options, all):
        # Several dumpers (JSON, DB backends based on JSON etc) don't support non-string keys...
        return {
            "adiabatic": {str(key):value.dump(digichem_options, all) for key,value in self.adiabatic.items()},
            "vertical": {str(key):value.dump(digichem_options, all) for key,value in self.vertical.items()}
        }
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        emissions = self()
        
        for transition_type in ("adiabatic", "vertical"):
            if transition_type not in data:
                continue
            
            setattr(emissions, transition_type,
                {mult: Relaxed_excited_state.from_dump(state_data, result_set, options, transition_type = transition_type) for mult, state_data in data[transition_type].items()}
            )
            
        return emissions


class Relaxed_excited_state(Excited_state):
    """
    Class for representing an emission energy from an excited state to a ground state.
    
    Note that while emission and absorption can be approximated as the reverse of each other, emission and absorption strictly have different energies.
    Excited states, as represented by the Excited_state class, are for absorption energies.
    This class is for emission energies, which are typically lower in energy (because the excited state has relaxed). 
    """
    
    def __init__(self,
            level,
            multiplicity,
            multiplicity_level,
            ground_multiplicity,
            excited_energy,
            ground_energy,
            oscillator_strength,
            transition_type,
            # TODO: Add Support.
            transitions = None,
            # TODO: Add Support.
            symmetry = None,
            transition_dipole_moment = None,
        ):
        """
        """
        self.ground_energy = ground_energy
        self.excited_energy = excited_energy
        self.ground_multiplicity = ground_multiplicity
        # The emission type (as a string), either fluorescence or phosphorescence.
        if ground_multiplicity == multiplicity:
            self.emission_type = "fluorescence"
        
        else:
            self.emission_type = "phosphorescence"
        
        # Either adiabatic or vertical.
        self.transition_type = transition_type
        
        super().__init__(level, multiplicity, multiplicity_level, symmetry, excited_energy - ground_energy, oscillator_strength, transitions if transitions is not None else [], transition_dipole_moment)
        
    
    @classmethod
    def from_results(self,
                 ground_state_result,
                 excited_state_result,
                 transition_type,
                 excited_state = None,
                 # For now we assume this is the lowest possible excited state (may change in future).
                 level = 1,
                 multiplicity_level = 1,
        ):
        """
        Constructor for Relaxed_excited_state objects.
        
        :param ground_state_result: A Result_set object representing the ground state.
        :param excited_state_result: A Result_set object representing the excited state.
        :param transition_type:  A string describing the type of transition, either 'adiabatic' (GS and ES relaxed) or 'vertical' (ES relaxed, GS @ ES geom).
        :param excited_state: An optional Excited_state object. This is required (for example) in time dependent DFT where the total energy of excited_state_result is the ground state energy (at the excited state geometry).
        :param level: The level (ordered index) of this excited state, this has no effect if excited_state is not None (in which case it is taken from the given excited state).
        :param multiplicity_level: The ordered index of this excited state out of states with the same multiplicity, this has no effect if excited_state is not None (in which case it is taken from the given excited state).
        """            
        if excited_state is not None:
            # If we have an excited state we can inherit certain properties from it.
            level = excited_state.level
            multiplicity_level = excited_state.multiplicity_level
            oscillator_strength = excited_state.oscillator_strength
            transition_dipole_moment = excited_state.transition_dipole_moment
        else:
            # We don't have a concept of oscillator strength (yet?).
            oscillator_strength = None
            transition_dipole_moment = None
        
        # The total energy of the excited state in this transition.
        # Start with our excited_state_result energy.
        excited_energy = excited_state_result.energies.final
        
        # Add the excited state energy (if we have it).
        if excited_state is not None:
            excited_energy += excited_state.energy
        
        # The total energy of the ground state in this transition.
        ground_energy = ground_state_result.energies.final
        
        # The multiplicity (as a number) of the excited state in this emission transition.
        if excited_state is not None:
            excited_multiplicity =  excited_state.multiplicity
        else:
            excited_multiplicity =  excited_state_result.ground_state.multiplicity
        
        # The multiplicity (as a number) of the ground state in this emission transition.
        ground_multiplicity = ground_state_result.ground_state.multiplicity
            
        return self(
            level = level,
            multiplicity = excited_multiplicity,
            multiplicity_level = multiplicity_level,
            ground_multiplicity = ground_multiplicity,
            # TODO: Add Support.
            symmetry = None,
            excited_energy = excited_energy,
            ground_energy = ground_energy,
            oscillator_strength = oscillator_strength,
            # TODO: Add Support.
            transitions = None,
            transition_dipole_moment = transition_dipole_moment,
            transition_type = transition_type,
        )
    
    @property
    def excited_multiplicity(self):
        """
        This is an alias of multiplicity
        """
        return self.multiplicity
        
    @classmethod
    def guess_from_results(self, *results):
        """
        Try and find emission energies from a number of calculation results.
        
        Attempts will be made to calculate both adiabatic and vertical emission energies based on the properties of the excited states in results.
        
        :param results: A number of result sets.
        :returns: A tuple of dictionaries of the form (vertical, adiabatic), where the key of each dict is the multiplicity of the corresponding emission.
        """        
        # The two types of emission we consider.
        # Each is a dictionary where the keys are multiplicities and the items are the corresponding emission.
        adiabatic = {}
        vertical = {}
        
        # We're only interested in optimisations (by definition emission needs to be from relaxed geom).
        opt_results = [result for result in results if "Optimisation" in result.metadata.calculations]
        
        for excited_state_result in reversed(opt_results):
            # See if we can use this calc type as an excited state.
            # The easiest situation is if we have excited states in the calc (and it is an opt).
            if len(excited_state_result.excited_states) > 0:
                # This will do fine for vertical.
                vertical.update(self.for_each_multiplicity(excited_state_result, excited_state_result, "vertical"))
                
                # For adiabatic, we need the earliest suitable optimisation from our list that doesn't have excited states.
                try:
                    ground_state_result = [result for result in opt_results if len(result.excited_states) == 0][0]
                    adiabatic.update(self.for_each_multiplicity(ground_state_result, excited_state_result, "adiabatic"))
                except IndexError:
                    # No suitable calcs.
                    pass                    
                    
            else:
                # We have no explicit excited states, but so long as we are an opt we can still find an emission if we have calcs of different mult.
                # First look for adiabatic.
                try:
                    # Try and find a suitable ground.
                    ground_state_result = [result for result in opt_results if result.ground_state.multiplicity != excited_state_result.ground_state.multiplicity][0]
                    
                    # Get 'emission' object.
                    emission = self.from_results(ground_state_result, excited_state_result, "adiabatic")
                    
                    # If the energy is negative we will ignore.
                    if emission.energy >= 0:
                        adiabatic[excited_state_result.ground_state.multiplicity] = emission
                except IndexError:
                    # No good.
                    pass
                
                # And now vertical (we need a single point at different geom.
                try:
                    # Try and find a suitable ground.
                    ground_state_result = [result for result in results if "Single Point" in result.metadata.calculations and result.ground_state.multiplicity != excited_state_result.ground_state.multiplicity][0]
                    
                    # Get 'emission' object.
                    emission = self.from_results(ground_state_result, excited_state_result, "vertical")
                    
                    # If the energy is negative we will ignore.
                    if emission.energy >= 0:
                        vertical[excited_state_result.ground_state.multiplicity] = emission
                except IndexError:
                    # No good.
                    pass
                
        return (vertical, adiabatic)
                
                
    @classmethod
    def for_each_multiplicity(self, ground_state_result, excited_state_result, transition_type):
        """
        Return one emission object for each multiplicity in an excited state result.
        
        This method is useful because frequently only the lowest excited state of each mult is considered for emission, so there is one emission type per multiplicity. 
        """
        return {multiplicity: self.from_results(ground_state_result, excited_state_result, transition_type, states[0]) for multiplicity, states in excited_state_result.excited_states.group().items()}
                
    # TODO: This might not be used anywhere?
    @classmethod
    def determine_from_results(self, main_result, *, ground_state_result = None, excited_state_result, transition_type = None, excited_state = None):
        """
        Try and determine some emission energies from a number of result sets.
        
        :param main_result: A Result_set object containing main calculation results.
        :param ground_state_result: A Result_set object representing the ground state.
        :param excited_state_result: A Result_set object representing the excited state.
        :param transition_type:  A string describing the type of transition, either 'adiabatic' (GS and ES relaxed) or 'vertical' (ES relaxed, GS @ ES geom).
        :param excited_state: An optional Excited_state object. This is required (for example) in time dependent DFT where the total energy of excited_state_result is the ground state energy (at the excited state geometry). If excited_state is not an Excited_state object (and isn't None), it is passed as criteria to Excited_state_list.get_state(). 
        """                
        # Decide on what to use for our ground state if not given explicitly.
        if ground_state_result is None:
            # If a transition type wasn't given, we can try and guess.
            if transition_type is None:
                # We don't know what type of emission we've been asked for.
                # If the main_result is an opt, we'll assume we're adiabatic and use that as the ground state.
                # Otherwise, we'll assume vertical and use the excited state.
                if "Optimisation" in main_result.metadata.calculations:
                    transition_type = "adiabatic"
                else:
                    transition_type = "vertical"
            
            # Distinguish between the two types of excited states
            if len(excited_state_result.excited_states) > 0:
                # Emission from TD style excited states.
                
                if transition_type == "vertical":
                    # For vertical emission, the excited state result is also the ground state result.
                    ground_state_result = excited_state_result
                    
                elif transition_type == "adiabatic":
                    # For adiabatic we need an optimised ground state geometry.
                    if "Optimisation" in main_result.metadata.calculations:
                        ground_state_result = main_result
                        
                    else:
                        # No explicit ground given and out main result is not an Opt, so it probably isn't suitable.            
                        raise Digichem_exception("Unable to determine ground state in adiabatic emission; no explicit ground state given and this calculation is not an optimisation")
                        
            else:
                # Emission from unrestricted triplet calc.
                # This method requires another calc type to compare to.
                if transition_type == "vertical":
                    # For vertical, we want a singlet calc at the same geometry as the triplet calc.
                    if 'Single Point' in main_result.metadata.calculations:
                        ground_state_result = main_result
                    else:
                        raise Digichem_exception("Unable to determine ground state in vertical emission; no explicit ground state given and this calculation is not a single point")
                    
                elif transition_type == "adiabatic":
                    # For adiabatic, we need an opt.
                    if "Optimisation" in main_result.metadata.calculations:
                        ground_state_result = main_result
                    else:
                        # No explicit ground given and out main result is not an Opt, so it probably isn't suitable.            
                        raise Digichem_exception("Unable to determine ground state in adiabatic emission; no explicit ground state given and this calculation is not an optimisation")
            
        # Now decide which excited state to use.
        excited_state = excited_state if excited_state is None or isinstance(excited_state, Excited_state) else excited_state_result.excited_states.get_state(excited_state)
        
        # No excited state given and generating emission from TD style excited states, get all possible.
        if len(excited_state_result.excited_states) > 0 and excited_state is None:
            return self.for_each_multiplicity(ground_state_result, excited_state_result, transition_type)
        
        else:
            # Get emission.
            emission = self.from_results(
                ground_state_result,
                excited_state_result,
                transition_type,
                excited_state)
            
            # Return as single dict.
            return {emission.multiplicity: emission}
    
   
        
    @property
    def emission_rate(self):
        """
        The rate of emission (k(F) or k(Phos) etc) from this excited state.
        
        Calculated according to Shizu, K., Kaji, H. Commun Chem 5, 53 (2022).
        """
        try:
            return ((4 * self.joules **3) / (3 * epsilon_0 * h_bar **4 * c **3) ) * self.transition_dipole_moment.coulomb_meters **2
        except Exception:
            if self.transition_dipole_moment is None:
                raise Result_unavailable_error("emission_rate", "there is no transition dipole moment associated with this emission") from None
            
            else:
                raise
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        
        dump_dict = super()._dump_(digichem_options, all)
        dump_dict['emission_type'] = self.emission_type
        dump_dict['ground_multiplicity'] = self.ground_multiplicity
        dump_dict['excited_energy'] = {
            "value": float(self.excited_energy),
            "units": "eV"
            }
        dump_dict['ground_energy'] = {
            "value": float(self.ground_energy),
            "units": "eV"
            }
        
        
        try:
            dump_dict["emission_rate"] = {
                "value": float(self.emission_rate),
                "units": "s^-1"
            }
        
        except Result_unavailable_error:
            pass
        
        return dump_dict
    
    @classmethod
    def from_dump(self, data, result_set, options, transition_type):
        """
        Get an instance  of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        # Get our tdm (if available).
        if 'tdm' in data and data['tdm'] is not None:
            tdm = result_set.transition_dipole_moments[data['index']-1]
        
        else:
            tdm = None
        
        return self(
            level = data['index'],
            multiplicity = data['multiplicity'],
            multiplicity_level = data['multiplicity_index'],
            ground_multiplicity = data['ground_multiplicity'],
            excited_energy = data['excited_energy']['value'],
            ground_energy = data['ground_energy']['value'],
            oscillator_strength = data['oscillator_strength'],
            transition_type = transition_type,
            transition_dipole_moment = tdm,
        )
        
        
