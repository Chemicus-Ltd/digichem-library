"""
Contains classes for representing the ground state of a system.

This file is closely related to the excited_states.py file.
These definitions may change or move.
"""

from digilib.result.excited_state import Energy_state


class Ground_state(Energy_state):
    """
    Class for representing the ground state of a system/molecule.
    """
    
    def __init__(self, charge, multiplicity, energy):
        """
        Constructor for ground state objects.
        
        :param charge: The charge of the ground states as a signed integer.
        :param multiplicity: The multiplicity as an integer of float.
        :param energy: The absolute energy of the GS in eV.
        """
        # The 'level' refers to the ordering of energy states, the ground state is always the lowest so have level 0.
        # Similarly, the multiplicity_level refers to the ordering of energy states that share the same multiplicity. Regardless of what multiplicity our GS has, we will be the lowest by definition.
        # Energy is our energy above the GS, so that will always be 0.
        super().__init__(level = 0, multiplicity = multiplicity, multiplicity_level = 0, energy = 0.0)
        self.charge = charge
        # Also save our absolute energy.
        self.absolute_energy = energy
        
    def __eq__(self, other):
        """
        Is this ground state object equal to another?
        """
        return self.charge == other.charge and self.multiplicity == other.multiplicity and self.energy == other.energy
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        parent_dict = super()._dump_(digichem_options, all)
        return {
            "index": parent_dict['index'],
            "symbol": parent_dict['symbol'],
            "charge": self.charge,
            "multiplicity": parent_dict['multiplicity'],
            "multiplicity_index": parent_dict['multiplicity_index'],
            # The 'energy' field from Energy_state is an excited state energy, relative to the ground state (so is always 0).
            # Replace with total energy.
            "energy": {
                "value": float(self.absolute_energy),
                "units": "eV"
            }
        }
    
    @classmethod
    def from_parser(self, parser):
        """
        Create a Ground_state object from an output file parser.
        
        As different energy types can be present in a single calculation, the following (arbitrary?) precedence is used:
        1) Coupled-cluster energy.
        2) The highest level Moller-Plesset energy.
        3) SCF energy.
        
        :param parser: An output file parser.
        """
        return self.from_energies(parser.results.metadata.charge, parser.results.metadata.multiplicity, parser.results.energies)
    
    @classmethod
    def from_dump(self, data, result_set, options):
        return self.from_energies(data['charge'], data['multiplicity'], result_set.energies)
    
    @classmethod
    def from_energies(self, charge, multiplicity, energies):
        """
        Create a ground state object selecting from various lists of energies.
        
        :param charge: Charge of the ground state.
        :param multiplicity: Multiplicity of the ground state.
        :param energies: List of energies.
        """
        if len(energies.cc) > 0:
            energy = energies.cc.final
        elif len(energies.mp) > 0:
            energy = energies.mp.final
        elif len(energies.scf) > 0:
            energy = energies.scf.final
        else:
            # We have no energies that we can understand. Give up.
            return None
        # Return our constructed object.
        return self(charge, multiplicity, energy)