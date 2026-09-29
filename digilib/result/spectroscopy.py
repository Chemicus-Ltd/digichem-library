import math
import numpy
import itertools
import scipy.constants
import scipy.signal

from digilib.misc.base import transpose
from digilib.exception.base import Digichem_exception
import digilib.result.excited_state

class Spectroscopy_graph_abc():
    """
    ABC for classes used to plot Gaussians over spectroscopic results. 
    """
    
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        raise NotImplementedError("Implement in subclass")
    
    def peaks(self):
        """
        Find peaks in the cumulative graph.
        
        :return: A list of x,y coords that are peaks.
        """
        coords = self.plot_cumulative_gaussian()
        y_coords = [y for x,y in coords]
        indexes = scipy.signal.find_peaks(y_coords, height = max(y_coords) * self.cutoff)[0]
        return [coords[index] for index in indexes]
        
    def plot_gaussian(self):
        """
        Plot a gaussian distribution around a set of coordinates.
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        raise NotImplementedError("Implement in subclass")
    
    def plot_cumulative_gaussian(self):
        """
        Plot an additive gaussian distribution around a set of coordinates.
        
        :return: A single list of tuples of (x, y) coordinates plotted by the gaussian function.
        """
        # First, get our gaussian plots.
        gplots = self.plot_gaussian()
        
        # Sum our y values.
        coords = {}
        for plot in gplots:
            for x, y in plot:
                coords[x] = coords.get(x, 0) + y
                
        # Now return as list of tuples.
        return sorted(list(coords.items()), key = lambda coord: coord[0])
    
    @classmethod
    def gaussian(self, a, b, c, x):
        """
        An implementation of the gaussian function.
        
        :param a: The maximum height of the peak.
        :param b: The x position of the center of the peak.
        :param c: The width of the peak.
        :param x: The x-value to plot for.
        :return: The corresponding y value.
        """
        return a * math.exp(-( (x - b)**2 / (2 * c**2) ))
    
    @classmethod
    def gaussian_x(self, a, b, c, y):
        """
        An implementation of the gaussian function rearranged to find x.
        
        :param a: The maximum height of the peak.
        :param b: The x position of the center of the peak.
        :param c: The width of the peak.
        :param y: The y-value to plot for.
        :return: A tuple of the two corresponding x values.
        """
        # We have two solutions of the form x = ± d + b
        # Calculate d first.
        try:
            d = math.sqrt( -math.log( y/a ) * 2 * c **2  )
        except ZeroDivisionError:
            # a (max height) is zero; the intensity is zero.
            # We could instead return b?
            raise Digichem_exception("'{}' cannot calculate Gaussian function limits; a (intensity) is zero".format(self.__name__))
        
        # Now return our two solutions.
        return (-d + b, d + b)
    
    @classmethod
    def fwhm_to_c(self, fwhm):
        """
        Convert a full width at half maximum to the corresponding c-value in the Gaussian function.
        
        :param fwhm: The desired full width at half maximum (in nm).
        :return: The c-value.
        """
        return fwhm / (2 * math.sqrt(2 * math.log(2)))


class Spectroscopy_graph(Spectroscopy_graph_abc):
    """
    Top level class for graphing spectroscopy results (energy vs intensity).
    
    For generating pictures of these graphs, see digilib.image.spectroscopy
    """
    
    def __init__(self, coordinates, fwhm, resolution = 1, cutoff = 0.01, filter = 1e-6, adjust_zero = False):
        """
        Constructor for Spectroscopy_graph objects
        
        :param coordinates: A list of (energy, intensity) tuples to plot. The units of energy and intensity are irrelevant here (but should be consistent). Note that coordinates with 0 intensity will be removed.
        :param fwhm: The full-width at half-maximum to plot peaks with (in units of the x axis).
        :param resolution: The spacing (or step-size) between points to plot using the gaussian function, in units of the x-axis.
        :param cutoff: The minimum y value to plot using the gaussian function, as a fraction of the intensity.
        :param adjust_zero: If True and all y values are 0, set all y values to 1 (so that something can be plotted).
        """
        # We save our coordinates under two properties.
        # base_coordinates are untransformed.
        # coordinates (which is by default the same as base_coordinates) are transformed.        
        self.false_intensity = False
        # Set arb y value if we've been asked to.
        if adjust_zero:
            # Check if all the y axis is zero.
            if all([y == 0 for x,y in coordinates]):
                coordinates = [(x, 1) for x, y in coordinates]
                self.false_intensity = True
        
        self.base_coordinates = [(x, y) for x,y in coordinates if y != 0]
        self.fwhm = fwhm
        self.resolution = resolution
        self.cutoff = cutoff
        self.filter = filter
        
        # Caches
        self._gaussians = []
        self._cumulative_gaussians = []
        
    @classmethod
    def from_vibrations(self, vibrations, *args, **kwargs):
        """
        Alternative constructor from a Vibrations_list object.
        """
        return self([(vibration.frequency, vibration.intensity) for vibration in vibrations], *args, **kwargs)
    
    @classmethod
    def from_excited_states(self, excited_states, *args, **kwargs):
        """
        An alternative constructor that takes a list of excited states as argument.
        
        :param excited_states: An Excited_states_list object to construct from.
        :param adjust_zero: If all the intensities of the given excited states are zero, whether to arbitrarily set the y coords to 1.
        :param **kwargs: Passed to the real constructor.
        """
        coords = [(excited_state.energy, excited_state.oscillator_strength if excited_state.oscillator_strength is not None else 0) for excited_state in excited_states]
        return self(coords, *args, **kwargs)
    
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        return self.base_coordinates
        
    def plot_gaussian(self, refresh = False):
        """
        Plot a gaussian distribution around a set of coordinates.
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        if len(self._gaussians) == 0 or refresh:
            self._gaussians = self.get_gaussian()
        
        return self._gaussians
        
    def get_gaussian(self):
        # First, determine our c value.
        c = self.fwhm_to_c(self.fwhm)
        
        # Next, determine the limits in which we'll plot.
        limits = self.gaussian_limits()
        
        # Apply rounding to nearest resolution step to cover up floating-point errors.
        points = [self.resolution *round(point / self.resolution) for point in numpy.linspace(*limits)]
        
        # Plot and return.
        digilib.log.get_logger().debug(
            "Plotting gaussian peaks from {:0.2f} to {:0.2f} with a step size of {} ({} total points) for {} peaks ({} total iterations)".format(
                limits[0],
                limits[1],
                self.resolution,
                limits[2],
                len(self.base_coordinates),
                len(self.base_coordinates) * limits[2]
            )
        )
        gaussians = [
            [(x, y) for x in points if (y  := self.gaussian(a, b, c, x)) >= self.filter]
            for b, a in self.base_coordinates
        ]
        
        return gaussians
    
    def plot_cumulative_gaussian(self, refresh = False):
        """
        Plot an additive gaussian distribution around a set of coordinates.
        
        :return: A single list of tuples of (x, y) coordinates plotted by the gaussian function.
        """
        if len(self._cumulative_gaussians) == 0 or refresh:
            self._cumulative_gaussians = self.get_cumulative_gaussian()
        
        return self._cumulative_gaussians
        
    def get_cumulative_gaussian(self):
        # First, determine our c value.
        c = self.fwhm_to_c(self.fwhm)
        
        # Next, determine the limits in which we'll plot.
        limits = self.gaussian_limits()
        
        # Apply rounding to nearest resolution step cover up floating-point errors.
        points = [self.resolution *round(point / self.resolution) for point in numpy.linspace(*limits)]
        
        # Plot and return.
        digilib.log.get_logger().debug("Plotting cumulative gaussian peaks from {:0.2f} to {:0.2f} with a step size of {} ({} total points) for {} peaks ({} total iterations)".format(limits[0], limits[1], self.resolution, limits[2], len(self.base_coordinates), len(self.base_coordinates) * limits[2]))
        gaussian = [
            (x, y) for x in points if ( y := sum((self.gaussian(a, b, c, x) for b, a in self.base_coordinates))) > self.filter
        ]
        
        return gaussian
    
    def gaussian_limits(self):
        """
        Determine min and max x limits to plot a gaussian function.
        
        :param c: The width of the peak.
        :param cutoff: The minimum y value to plot using the gaussian function, as the fraction of a.
        :param resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :return: A tuple of (minlim, maxlim, num) where minlim is the most negative value, maxlim the most positive and num the integer number of points to plot between them to achieve resolution.
        """
        # Calculate limits for each set of coordinates given to us.
        all_limits = [self.gaussian_x(y, x, self.fwhm_to_c(self.fwhm), self.cutoff * y) for x, y in self.base_coordinates]
        
        try:
            limits = (min(itertools.chain.from_iterable(all_limits)), max(itertools.chain.from_iterable(all_limits)))
        
        except ValueError:
            if len(list(itertools.chain.from_iterable(all_limits))) == 0:
                # Nothing to plot.
                raise Digichem_exception("'{}' cannot plot spectrum; there are no values".format(type(self).__name__))
            else:
                raise
            
        # Align our start and end points to an integer multiple of resolution.
        # This is necessary to facilitate easy summation of gaussians (because they all align to the same grid,
        # assuming the same value of resolution)
        start = math.floor(limits[0] / self.resolution) * self.resolution
        end = math.ceil(limits[1] / self.resolution) * self.resolution
        
        # Now we need to generate the x values which we'll plot for.
        # This is a little more complicated than it needs to be because there's no simple range() for floats.
        # Calculate the number of points (and round up).
        num_points = round( math.fabs(end - start) / self.resolution)
        
        # Extend our limits so they are a clean multiple of num_points.
        limits = (start, start + num_points * self.resolution)

        return (limits[0], limits[1], num_points +1)
    
    
class Absorption_emission_graph(Spectroscopy_graph):
    """
    Class for graphing absorption/emission spectra in nm.
    
    Use the normal Spectroscopy_graph object for plotting in eV.
    """
    
    def __init__(self, *args, use_jacobian = True, **kwargs):
        """
        Constructor for Absorption_emission_graph objects.
        
        :param false_intensity: A flag to indicate that the intensity units (y-axis) are arbitrary.
        :param use_jacobian: Whether to use the jacobian transform to scale the y axis.
        """
        # Get our x,y values from our excited states.
        super().__init__(*args, **kwargs)
        self.use_jacobian = use_jacobian
        
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        return self.nm_coordinates
    
    @property
    def nm_coordinates(self):
        """
        Get the list of coordinates of this graph scaled to nm.
        """
        return [self.energy_to_wavelength(coord, self.use_jacobian) for coord in self.base_coordinates]
    
    @classmethod
    def energy_to_wavelength(self, coord, use_jacobian):
        """
        Convert a pair of x,y coordinates in ev to nm.
        
        :param coord: A tuple of (eV, f).
        :return: A tuple of (nm, i), where i is intensity.
        """
        x, y = coord
        return (digilib.result.excited_state.Excited_state.energy_to_wavelength(x), self.jacobian(x, y) if use_jacobian else y)
        
    @classmethod
    def jacobian(self, E, f_E):
        """
        An implementation of the jacobian transform that scales intensity in energy units to intensity in wavelength units.
        
        See J. Phys. Chem. Lett. 2014, 5, 20, 3497 for why this is necessary.
        
        Note that the jacobian transform will maintain the area under the curve regardless of x units (nm or x).
        Sadly, this has the consequence of mangling the intensity units (it becomes tiny; an oscillator strength of 1 at 3 eV becomes 1.163e-12).
        """
        return ((E * scipy.constants.electron_volt)**2  * f_E) / (scipy.constants.Planck * scipy.constants.c)
    
    @classmethod
    def inverse_jacobian(self, E, f_nm):
        """
        An implementation of the jacobian transform that scales intensity in wavelength units to intensity in energy units.
        
        See J. Phys. Chem. Lett. 2014, 5, 20, 3497 for why this is necessary.
        
        Note that the jacobian transform will maintain the area under the curve regardless of x units (nm or x).
        Sadly, this has the consequence of mangling the intensity units (it becomes tiny; an oscillator strength of 1 at 3 eV becomes 1.163e-12).
        """
        # TODO: Might be better to rearrange this to accept nm rather than eV?
        return (
            (f_nm * scipy.constants.Planck * scipy.constants.c) / (E * scipy.constants.electron_volt) **2
        )
    
    @classmethod
    def shift_coord(self, coord, delta_eV):
        """
        Shift a coordinate (in nm) by a given energy value.

        :param delta_eV: The energy (in eV) to shift by. A positive value will blueshift (higher energy).
        """
        old_x_nm, old_y_nm = coord
        # Convert x to energy.
        old_x_ev = digilib.result.excited_state.Excited_state.wavelength_to_energy(old_x_nm)
        # Transform y.
        old_y_ev = self.inverse_jacobian(old_x_ev, old_y_nm)

        # Shift by given amount.
        new_x_ev = old_x_ev + delta_eV

        # Convert back to nm.
        new_x_nm, new_f_nm = self.energy_to_wavelength((new_x_ev, old_y_ev), True)

        return (new_x_nm, new_f_nm)
    
    def shift(self, delta_eV):
        """
        """
        return map(lambda coord: self.shift_coord(coord, delta_eV), self.coordinates)

    
    def plot_gaussian(self):
        """
        Plot a gaussian distribution around our excited state energies.
        
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        # All we need to do over our parent is convert x values from e to wavelength.
        # And scale y values using the jacobian transform.
        return [[self.energy_to_wavelength(coord, self.use_jacobian) for coord in plot] for plot in super().plot_gaussian()]
    
    def plot_cumulative_gaussian(self):
        """
        Plot an additive gaussian distribution around a set of coordinates.
        
        :return: A single list of tuples of (x, y) coordinates plotted by the gaussian function.
        """
        return [self.energy_to_wavelength(coord, self.use_jacobian) for coord in super().plot_cumulative_gaussian()]


# TODO: Find this a better home?
def unpack_coupling(couplings, satellite_threshold = 0.02):
    """
    Unpack a nested dict of dict of NMR_group_spin_coupling objects into a single, ordered dict.
    
    The keys of the returned dict will be a tuple of (foreign_atom_group, foreign_isotope).
    """
    couplings = {
            (coupling_group, coupling_isotope): isotope_coupling for coupling_group, atom_dict in couplings.items() for coupling_isotope, isotope_coupling in atom_dict.items()
            if coupling_group.element[coupling_isotope].abundance / 100 > satellite_threshold
        }
    # Sort couplings.
    couplings = dict(sorted(couplings.items(), key = lambda coupling: abs(coupling[1].total), reverse = True))
    
    return couplings


class NMR_graph(Spectroscopy_graph):
    """
    A class for plotting a single NMR peak.
    
    A 'single' NMR peak here refers to the signal that would be observed from only one atom (or atom group).
    For plotting an entire spectrum, use a Combined_graph of NMR_graph objects.
    """
    
    def multiplicity(self, atom_group, coupling, satellite_threshold = 0.05):
        """
        Determine the multiplicity of this peak.
        
        Note that the returned multiplicity might not match exactly the observed multiplicity of the peak
        due to line-broadening and overlapping signals.
        
        :param atom_group: The atom_group that this peak corresponds to.
        :param coupling: NMR_group_spin_coupling objects between this atom_group and other groups.
        """        
        # First, determine how many peaks are visible.
        peaks = self.peaks()
        
        # Filter out peaks which are significantly smaller than the tallest peak.
        max_peak = max(transpose(peaks, 2)[1])
        peaks = [peak for peak in peaks if peak[1] > max_peak * satellite_threshold]
        
        if len(peaks) == 1:
            return [{"symbol": "s", "number": 1, "multiplicity": "singlet"}]
        
        # A list of multiplicity dicts.
        mults = []
        total_peaks = 1
        coupling_index = 0
        
        # Get couplings.
        couplings = unpack_coupling(coupling)
        couplings = list(couplings.values())
        
        # Unless coupling has been calculated for all available isotopes for each atom group (unlikely),
        # there will be one additional peak for each atom group from non-NMR active nuclei.
        # Make sure we don't count this peak multiple times.
        #
        # For atoms in which the majority of the abundance is already accounted for (1H, for example),
        # exclude the residual peak.
        residual_isotope_peaks = set(
            [a_group for a_group, isotopes in coupling.items() if sum((a_group.element[isotope].abundance for isotope in isotopes.keys())) / 100 > satellite_threshold]
        )
        
        # We will keep requesting more splitting until we are able to account for all the peaks we can see.
        while len(peaks) > total_peaks:
            # Get the next coupling available.
            try:
                group_coupling = couplings[coupling_index]
            
            except IndexError:
                # Ran out of coupling.
                break
            
            other_group = group_coupling.groups[group_coupling.other(atom_group)]
            residual = other_group not in residual_isotope_peaks
            multiplicity = group_coupling.multiplicity(atom_group)
            residual_isotope_peaks.add(other_group)
            
            mults.append(multiplicity)
            num = multiplicity['number']
            if residual:
                num += 1
            
            total_peaks *= num
            
            # Update counter.
            coupling_index += 1
        
        return mults
        
    @classmethod
    def from_nmr(self, nmr_peaks, *args, **kwargs):
        """
        An alternative constructor that takes a list of simulated NMR peaks as argument.
        
        :param nmr_peaks: Simulated NMR peaks. See result.nmr.NMR_spectrometer for how to obtain this data.
        """
        return self(nmr_peaks, *args, **kwargs)


class Combined_graph(Spectroscopy_graph_abc):
    """
    A class for plotting one graph from multiple, individual graphs.
    """
    
    def __init__(self, graphs):
        """
        Constructor for Combined_graph objects.
        
        :param graphs: A dictionary of individual Spectroscopy_graph objects to combine.
        """
        self.graphs = graphs
        
    @property
    def cutoff(self):
        """
        """
        return max((graph.cutoff for graph in self.graphs.values()))
    
    @property
    def fwhm(self):
        """
        """
        return sum((graph.fwhm for graph in self.graphs.values())) / len(self.graphs)
        
    @classmethod
    def from_nmr(self, grouped_peaks, *args, **kwargs):
        """
        Construct a Combined_graph from a dictionary of grouped NMR peaks
        
        :param grouped_peaks: A dictionary of lists of NMR peaks. Each key in the dictionary should correspond to one atom group.
        """
        return self(
            {
                peak_key: NMR_graph.from_nmr(peaks, *args, **kwargs) for peak_key, peaks in grouped_peaks.items()
            }
        )
    
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        return list(itertools.chain(*[graph.coordinates for graph in self.graphs.values()]))
        
    def plot_gaussian(self):
        """
        Plot a gaussian distribution around a set of coordinates.
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        # Because this class represents one additional layer of abstractions,
        # we'll combine the results returned by each sub graph together, to obtain a single nested list.
        return [
            graph.plot_cumulative_gaussian()
            for graph in self.graphs.values()
        ]
    