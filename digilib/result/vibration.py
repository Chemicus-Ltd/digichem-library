# General imports.
from itertools import zip_longest
import warnings

# Digichem imports.
from digilib.result import Result_container
from digilib.result import Result_object
from digilib.result import Floatable_mixin, Unmergeable_container_mixin

# Hidden.
#from digilib.result.spectroscopy import Spectroscopy_graph


class Vibrations_list(Result_container, Unmergeable_container_mixin):
    """
    Class for representing a group of molecular vibrations.
    """
    
    MERGE_WARNING = "Attempting to merge lists of vibrations that are not identical; non-equivalent vibrations will be ignored"
        
    @property
    def negative_frequencies(self):
        """
        Get a Vibrations_list object of the vibrations in this list that have negative frequencies (they are imaginary).
        """
        warnings.warn("negative_frequencies is deprecated, use negative instead")
        return self.negative
    
    @property
    def negative(self):
        """
        Get a Vibrations_list object of the vibrations in this list that have negative frequencies (they are imaginary).
        """
        return type(self)([vibration for vibration in self if vibration.frequency < 0])
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an Vibrations_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: A Vibrations_list object. The list will be empty if no vibration frequency data is available.
        """
        return self(Vibration.list_from_parser(parser))
    
    def _get_dump_(self):
        """
        Method used to get a dictionary used to generate on-demand values for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
        """
        return {
            "spectrum": self.generate_spectrum,
        }
    
    def generate_spectrum(self, digichem_options):
        """
        Abstract function that is called to generate an on-demand value for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        :param key: The key being requested.
        :param digichem_options: Digichem options being used to dump.
        """
        from digilib.result.spectroscopy import Spectroscopy_graph
        
        
        spectrum = Spectroscopy_graph.from_vibrations(
            self,
            fwhm = digichem_options['IR_spectrum']['fwhm'],
            resolution = digichem_options['IR_spectrum']['gaussian_resolution'],
            cutoff = digichem_options['IR_spectrum']['gaussian_cutoff'],
            filter = digichem_options['IR_spectrum']['y_filter'],
        )
        
        
        try:
            spectrum_data = spectrum.plot_cumulative_gaussian()
            
            spectrum_data = [{"x":{"value": float(x), "units": "c m^-1"}, "y": {"value":float(y), "units": "km mol^-1"}} for x,y in spectrum_data]
            spectrum_peaks = [{"x":{"value": float(x), "units": "c m^-1"}, "y": {"value":float(y), "units": "km mol^-1"}} for x, y in spectrum.peaks()]
        
        except Exception:
            spectrum_data = []
            spectrum_peaks = []
        
        return {
            "values": spectrum_data,
            "peaks": spectrum_peaks
        }
    
    def _dump_(self, digichem_options, all):
        dump_dict = {
            "num_vibrations": len(self),
            "num_negative": len(self.negative),
            "values": super()._dump_(digichem_options, all),
        }
        return dump_dict

    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(Vibration.list_from_dump(data['values'], result_set, options))


class Vibration(Result_object, Floatable_mixin):
    """
    Class for representing vibrational frequencies.
    """
    
    def __init__(self, level, frequency, intensity, symmetry = None):
        """
        Vibration object constructor.
        
        :param level: The numerical index of the vibration.
        :param frequency: Frequency of the vibration (cm-1).
        :param intensity: Intensity of the vibration in IR (km mol-1)
        :param symmetry: Symmetry term of the vibration (string).
        """
        self.level = level
        self.symmetry = symmetry
        self.frequency = frequency
        self.intensity = intensity
        
    def __float__(self):
        """
        A float representation of this object.
        """
        return float(self.frequency)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "index": self.level,
            "frequency": {
                "value": float(self.frequency),
                "units": "c m^-1",
            },
            "intensity": {
                "value": float(self.intensity),
                "units": "km mol^-1"
            },
            "symmetry": self.symmetry
        }
        
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return [self(vib_dict['index'], vib_dict['frequency']['value'], vib_dict['intensity']['value'], vib_dict['symmetry']) for vib_dict in data]
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of Vibration objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of Vibration objects. The list will be empty if no frequency data is available.
        """
        try:
            return [self(index+1, frequency, intensity, symmetry) for index, (symmetry, frequency, intensity) in enumerate(zip_longest(getattr(parser.data, 'vibsyms', []), parser.data.vibfreqs, parser.data.vibirs, fillvalue = None))]
        except AttributeError:
            return []