"""Classes for total system energies."""
import itertools
import scipy.constants

from digichem.exception.base import Result_unavailable_error
from digichem.result import Result_container
from digichem.result import Unmergeable_container_mixin
from digichem.result import Result_object

class Energies(Result_object):
    """
    Class that holds all the energies calculated from a computation.
    """
    
    def __init__(self, scf = None, mp = None, cc = None):
        """
        :param scf: List of self-consistent field energies.
        :param mp: List of lists of Moller-Plesset energies.
        :param cc: List of coupled cluster energies.
        """
        self.scf =  SCF_energy_list(scf if scf is not None else [])
        self.mp_energies = mp
        self.cc = CC_energy_list(cc if cc is not None else [])
    
    @property
    def mp(self):
        try:
            return self.mp_energies[-1]
        
        except IndexError:
            return MP_energy_list([], 0)
        
    def __iter__(self):
        for energy in (self.scf, self.mp, self.cc):
            yield energy
        
    @property
    def final(self):
        """
        The total energy of this calculation.
        
        This convenience property is the energy at the highest level of theory available (CC > MP > SCF).
        
        :raises Result_unavailable_error: If no total energy is available.
        """
        # Try CC first.
        if len(self.cc) > 0:
            return self.cc.final
        elif len(self.mp) > 0:
            return self.mp.final
        else:
            return self.scf.final
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an Energies object from an output file parser.
        
        :param parser: An output file parser.
        :return: The populated Energies object.
        """
        return self(
            cc = CC_energy_list.from_parser(parser),
            mp = MP_energy_list.list_from_parser(parser),
            scf = SCF_energy_list.from_parser(parser)
        )
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        try:
            final = float(self.final)
        
        except Result_unavailable_error:
            final = None
        
        dump = {
            "final": {
                "value": final,
                "units": "eV"
            },
            "scf": self.scf.dump(digichem_options, all),
            "mp": self.mp.dump(digichem_options, all),
            "cc": self.cc.dump(digichem_options, all)
        }
        dump.update({"mp{}".format(index+2): energy.dump(digichem_options, all) for index, energy in enumerate(self.mp_energies)})
        return dump
        
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        scf = SCF_energy_list.from_dump(data['scf'], result_set, options)
        #mp = MP_energy_list.from_dump(data['mp'], result_set, options)
        cc = CC_energy_list.from_dump(data['cc'], result_set, options)
        
        # TODO: consider using a dict so list indices make sense?
        mp_energies = []
        
        # Other MP energies.
        for mp_order in itertools.count(2):
            if "mp{}".format(mp_order) in data:
                mp_energies.append(MP_energy_list.from_dump(data["mp{}".format(mp_order)], result_set, options))
            
            else:
                break
        
        return self(scf, mp_energies, cc)


class Energy_list(Result_container, Unmergeable_container_mixin):
    """
    Class that represents a list of calculated energies.
    
    Storing energies as a list is useful because optimisations will print energies at each step, and it can be useful to look back and see how the opt has progressed.
    """
    
    # The type of energy; corresponds to the name of the attribute provided by cclib.
    cclib_energy_type = ""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def __float__(self):
        """
        Our float representation, which is the last energy in our list.
        """
        return self.final
    
    def __str__(self):
        """
        Our string representation, same as the float but in text form.
        """
        return str(self.__float__())
    
    def assign_levels(self):
        """
        (Re)assign total levels of the objects in this list.
        
        The contents of this list will be modified in place.
        """
        # Energies don't have levels.
        pass
    
    @property
    def human_energy_type(self):
        """
        Human readable name of this type of energy (for example, self-consistent field).
        """
        if self.energy_type == "SCF":
            return "self-consistent field"
        elif self.energy_type == "MP":
            return "Møller–Plesset"
        elif self.energy_type == "CC":
            return "coupled-cluster"
        else:
            raise ValueError("Unrecognised energy type {}".format(self.energy_type))

    @property
    def final(self):
        """
        The 'final' energy in this list, useful for getting the final optimised energy.
        
        :raises Result_unavailable_error: If this list is empty.
        """
        try:
            return self[-1]
        except IndexError:
            raise Result_unavailable_error(self.energy_type + ' energy', 'there is no {} energy'.format(self.energy_type))
    
    @classmethod
    def eV_to_kJ(self, eV):
        """
        Convert a value in eV to kJ.
        
        :param eV: The value in electron volts.
        :return: The value in kilojoules
        """
        return eV * 1.602176634e-22
    
    @classmethod
    def kJ_to_kJmol(self, kJ):
        """
        Convert a value in kJ to kJmol-1.
        
        :param kJ: The value in kilojoules.
        :return: The value in kilojoules per mol
        """
        return kJ * scipy.constants.Avogadro
    
    @classmethod
    def eV_to_kJmol(self, eV):
        """
        Convert a value in eV to kJmol-1.
        
        :param eV: The value in electron volts.
        :return: The value in kilojoules per mol
        """
        return self.kJ_to_kJmol(self.eV_to_kJ(eV))

    @classmethod
    def from_parser(self, parser):
        """
        Get an Energy_list from an output file parser.
        
        :param parser: An output file parser.
        :return: The populated Energy_list object. The object will be empty if the energy is not available.
        """
        try:
            return self(getattr(parser.data, self.cclib_energy_type))
        except AttributeError:
            return self()
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "num_steps": len(self),
            "final": {
                "value": float(self.final) if len(self) > 0 else None,
                "units": "eV"
            },
            "values": [{"value": float(energy), "units": "eV"} for energy in self]
        }
        
    @classmethod
    def from_dump(self, data, result_set, options, **kwargs):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self([energy['value'] for energy in data['values']], **kwargs)


class SCF_energy_list(Energy_list):
    """
    List of Self-consistent field energies.
    """
    cclib_energy_type = "scfenergies"
    energy_type = "SCF"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of SCF energies that are not identical; non-equivalent energies will be ignored"


class CC_energy_list(Energy_list):
    """
    List of coupled-cluster energies.
    """
    cclib_energy_type = "ccenergies"
    energy_type = "CC"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of CC energies that are not identical; non-equivalent energies will be ignored"


class MP_energy_list(Energy_list):
    """
    List of Moller-Plesset energies.
    """
    cclib_energy_type = "mpenergies"
    energy_type = "MP"
    
    # A warning issued when attempting to merge non-equivalent objects
    MERGE_WARNING = "Attempting to merge list of MP energies that are not identical; non-equivalent energies will be ignored"
    
    def __init__(self, items, order):
        super().__init__(items)
        self.order = order
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump = super()._dump_(digichem_options, all)
        dump['order'] = self.order
        return dump
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return super().from_dump(data, result_set, options, order = data['order'])
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of MP_energy_list objects from an output file parser.
        
        Note that unlike other calculated energies, MP energies are calculated sequentially up to the requested order. So for example, MP4 first calculates the MP1 energy (which is the same as the uncorrected energy), then MP2, MP3 and finally MP4.
        
        :param parser: An output file parser.
        :return: The populated Energy_list object. The object will be empty if the energy is not available.
        """
        # First, split our list of lists (for each opt step and each MP energy).
        mp_energies = []
        
        try:
            parser_energies = getattr(parser.data, self.cclib_energy_type)
            
            for energy_step in parser_energies:
                for index, energy in enumerate(energy_step):
                    try:
                        mp_energies[index].append(energy)
                        
                    except IndexError:
                        mp_energies.append([energy])
        
        except AttributeError:
            if not hasattr(parser.data, self.cclib_energy_type):
                return []
            
            else:
                raise
        
        return [self(energies, index+2) for index, energies in enumerate(mp_energies)]
    