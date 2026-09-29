from matplotlib import ticker
import matplotlib.collections
import statistics
import adjustText

from digilib.exception.base import File_maker_exception
from digilib.result.spectroscopy import Spectroscopy_graph,\
    Absorption_emission_graph, unpack_coupling
from digilib.image.graph import Graph_image_maker


class Spectroscopy_graph_maker(Graph_image_maker):
    """
    A graph type for displaying spectroscopy data (UV-vis, IR etc.).
    """
    def __init__(
        self,
        output,
        graph,
        *,
        # Keyword only arguments from here:
        
        # Options that control what to show in the graph.
        plot_bars = True,
        plot_peaks = False,
        plot_cumulative_peak = True,
        
        # Options controlling the axes limits.
        x_padding = 50,
        peak_cutoff = 0.01,
        x_limits = 'auto',
        y_limits = 'auto',
        **kwargs
    ):
        """
        Constructor for spectroscopy graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        :param fwhm: The desired full width at half maximum of the peaks that we are going to plot.
        :param gaussian_cutoff: The minimum y value to plot using the gaussian function, as the fraction of the max intensity.
        :param gaussian_resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :param plot_bars: If True, vertical bars are plotted on the graph for each excited state.
        :param plot_peaks: If True, Gaussian peaks are plotted for each excited state.
        :param plot_cumulative_peak: If True, a cumulative graph is plotted from gaussian peaks from all excited states.
        :param x_padding: An amount, in units of the x-axis, to extend beyond the lowest & highest energy peak (so the graph is framed). Note that the graph will never extend beyond 0.
        :param y_cutoff: The minimum y value to include when using auto x limits (the x axis will extend
        :param x_limits_method: String controlling how the x axis limits are set. Options are 'auto' for automatic scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
        :param y_limits_method: String controlling how the y axis limits are set. Options are 'auto' for automatic scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
        
        """
        # Call our image maker parent.
        super().__init__(output, **kwargs)
        
        # Save our graph object (holds coordinates to plot).
        self.graph = graph
        
        # Options controlling what to plot.
        self.should_plot_bars = plot_bars
        self.should_plot_peaks = plot_peaks
        self.should_plot_cumulative_peak = plot_cumulative_peak
        
        self.init_image_width = 10
        self.init_image_height = 5.26
        
        # Axis scaling.
        self.x_padding = x_padding
        self.peak_cutoff = peak_cutoff
        self.x_limits_method = x_limits
        self.y_limits_method = y_limits
        
        # Line-widths.
        self.linewidth = 1
        self.cumulative_linewidth = 2
        self.column_linewidth = 1.5
        
    @property
    def fwhm(self):
        return self.graph.fwhm
        
    @classmethod
    def transpose(self, coordinates):
        """
        Transpose a list of coordinates from [(x1,y1), (x2,y2)...] to [[x1, x2...], [y1, y2...]]
        """
        return list(map(list, zip(*coordinates)))
    
    def plot(self):
        """
        Plot the contents of our graph.
        """
        # Make sure we actually have something to plot.
        if len(self.graph.coordinates) == 0:
            # We have nothing to plot; we could plot an empty graph (and maybe should?), but for now we'll stop.
            raise File_maker_exception(self, "cannot plot spectrum; there are no values with intensity above 0")
        
        # What we do next depends on what options we've been given.
        if self.should_plot_bars:
            self.plot_columns()
            # If we aren't going to plot any curves, plot some invisible points so we get automatic axis scaling.
            if not self.should_plot_peaks and not self.should_plot_cumulative_peak:
                self.plot_hidden_points()
            
        if self.should_plot_cumulative_peak:
            self.plot_cumulative_line()
                
        if self.should_plot_peaks:
            self.plot_lines()
                        
    def plot_columns(self):
        """
        Plot vertical columns for each plot on the graph we are building.
        """
        columns = [[(x, 0) , (x, y)] for x, y in self.graph.coordinates]
        self.axes.add_collection(matplotlib.collections.LineCollection(columns, linewidths = self.column_linewidth, colors = (0,0,0)))
        
    def plot_hidden_points(self):
        """
        Plot invisible points to take advantage of matplotlib auto scaling.
        """        
        # Plot. Set size (s) to zero so we don't see our points, but we still get good automatic scaling.
        self.axes.scatter(*self.transpose(self.graph.coordinates), s=0)
        
    def plot_lines(self):
        """
        Plot x and y values as individual lines on the graph we are building.
        """            
        for plot in self.graph.plot_gaussian():
            data = self.transpose(plot)
            self.axes.plot(*data, 'C0-', linewidth = self.linewidth)
            
    def plot_cumulative_line(self):
        """
        Plot x and y values as a single line on the graph we are building.
        """
        data = self.transpose(self.graph.plot_cumulative_gaussian())
        self.axes.plot(*data, 'C0-', linewidth = self.cumulative_linewidth)
    
    @property
    def peaks(self):
        """
        A list of peaks (in x units).
        """
        return [x for x,y in self.graph.peaks()]
    
    def selected_peaks(self, decimals = 0, number = None):
        """
        A unique list of sorted peaks, rounded to a given number of decimal points.
        
        :param decimals: The number of decimal points to round to.
        :param number: The number of peaks to return (the most intense peaks will be returned first).
        :returns: A list of ordered peaks as floats (if decimals > 0) or ints.
        """
        if number is not None:
            peaks = self.graph.peaks()
            peaks.sort(key = lambda coord: coord[1])
            peaks = [x for x,y in peaks[-number:]]
        
        else:
            peaks = self.peaks
        
        return sorted(list(set([round(peak, decimals) if decimals > 0 else int(peak) for peak in peaks])))
        
    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        # We need to get a list of all peaks that are above our cutoff point.
        # First determine our highest point.
        highest_point = max(self.transpose(self.graph.coordinates)[1])
        
        # Now filter by a fraction of that amount.
        visible_x_values = [x for x, y in self.graph.plot_cumulative_gaussian() if y >= (highest_point * self.peak_cutoff)]
        
        # Now we can just get our min-max, remembering to include our x-padding, and never extending beyond 0 (because negative energy has no meaning in this context(?)).
        self.axes.set_xlim(max(min(visible_x_values) - self.x_padding, 0), max(visible_x_values) + self.x_padding)
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        # Use matplotlib autoscaling.
        self.axes.autoscale(enable = True, axis = 'y')
        
        # Clamp to 0 -> pos.
        self.axes.set_ylim(0, self.axes.get_ylim()[1])
            
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """
        super().adjust_axes()
        
        # Call the appropriate scaling method depending on the value of the limits_method.
        for axis in ("x", "y"):
            method = getattr(self, axis + "_limits_method")
            if method == "auto":
                # Call the automatic function.
                getattr(self, "auto_{}_limits".format(axis))()
            
            elif isinstance(method, tuple) or isinstance(method, list):
                # Manual, explicit limits.
                getattr(self.axes, "set_{}lim".format(axis))(method[0], method[1])
            
            else:
                # Don't recognise method.
                raise ValueError("Unknown {} limits method '{}'".format(axis, method))
            
        # Use matplotlib to automatically layout our graph.
        self.figure.tight_layout()
        
        # Get a constant scale on our x axis.
        self.constant_scale(0, self.inch_per_x)
        

class Absorption_emission_graph_maker(Spectroscopy_graph_maker):
    """
    A graph for displaying simulated absorptions.
    """
    
    # The key/name of where our options are stored in the main config.
    options_name = None
    
    def __init__(self, output, graph, *args, **kwargs):
        """
        Constructor for UV-Vis type absorption graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        :param adjust_zero: If all the intensities of the given excited states are zero, whether to arbitrarily set the y coords to 1.
        :param peak_cutoff: The minimum oscillator strength a peak must have to be drawn, as a fraction of the tallest peak. Set to 0 for no cutoff.
        :param max_width: The maximum width in pixels of the graph.
        """        
        # Call our parent.
        super().__init__(output, graph, *args, **kwargs)
        
        # The amount of space (in matplotlib inches) to allocate per unit of the x axis.
        self.inch_per_x = 0.0192
        
        # The DPI to save our image with.
        self.output_dpi = 100
        
        # Axes labels.
        self.x_label = "Wavelength /nm"
        if not self.graph.false_intensity:
            self.y_label = "Intensity" if self.graph.use_jacobian else "Oscillator Strength"
        else:
            self.hide_y = True
    
    @classmethod
    def from_options(self, output, *, excited_states, options, adjust_zero = True, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            Absorption_emission_graph.from_excited_states(
                excited_states,
                options[self.options_name]['fwhm'],
                options[self.options_name]['gaussian_resolution'],
                options[self.options_name]['gaussian_cutoff'],
                use_jacobian = options[self.options_name]['use_jacobian'],
                filter = options[self.options_name]['y_filter'],
                adjust_zero = adjust_zero),
            **{key: value for key, value in options[self.options_name].items() if key not in ["gaussian_cutoff", "gaussian_resolution", "fwhm", "use_jacobian", "y_filter"]},
            **kwargs
        )
    
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """
        # Change spacing of tick markers.
        self.axes.xaxis.set_major_locator(ticker.MultipleLocator(50))
        self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))

        # Call parent.
        super().adjust_axes()
        
        # Check to see if we've exceeded our max width.
        if self.max_width is not None and (self.figure.get_size_inches()[0] * self.output_dpi) > self.max_width:
            # Adjust our tick labels.
            x_difference = abs(self.axes.get_xlim()[0] - self.axes.get_xlim()[1])
            if x_difference < 1200:
                self.axes.xaxis.set_major_locator(ticker.MultipleLocator(100))
                self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(50))
            else:
                self.axes.xaxis.set_major_locator(ticker.MultipleLocator(200))
                self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(100))
            
            # Update our layout.
            self.figure.tight_layout()
            
class Absorption_graph_maker(Absorption_emission_graph_maker):
    """
    Class for creating absorption (UV-Vis) spectra.
    """
    options_name = "absorption_spectrum"
    
class Emission_graph_maker(Absorption_emission_graph_maker):
    """
    Class for creating emission spectra.
    """
    options_name = "emission_spectrum"

class Frequency_graph_maker(Spectroscopy_graph_maker):
    """
    A graph for displaying vibrational frequencies.
    """
            
    def __init__(self, output, graph, *args, **kwargs):
        """
        Constructor for frequency graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        :param x_limits_method: String controlling how the x axis limits are set. Options are 'auto' for standard auto scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
        :param y_limits_method: String controlling how the y axis limits are set. Options are 'auto' for standard auto scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
        """
        # Call our parent.
        super().__init__(output, graph, *args, **kwargs)
        
        # Axes labels.
        self.x_label = "Frequency /cm$^{-1}$"
        self.y_label = "Intensity /km mol$^{-1}$"
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = 0.00252
        
        
    @classmethod
    def from_options(self, output, *, vibrations, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            graph = Spectroscopy_graph.from_vibrations(
                vibrations,
                options['IR_spectrum']['fwhm'],
                options['IR_spectrum']['gaussian_resolution'],
                options['IR_spectrum']['gaussian_cutoff'],
                filter = options['IR_spectrum']['y_filter']
            ),
            **{key: value for key, value in options['IR_spectrum'].items() if key not in ["gaussian_cutoff", "gaussian_resolution", "fwhm", "y_filter"]},
            **kwargs
        )

    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        # We use automatic X scaling, except we don't go past 0, we show at least 3500 cm-1 and our scale is reversed (this is typical for IR).
        self.axes.autoscale(enable = True, axis = 'x')
        self.axes.set_xlim(max(self.axes.get_xlim()[1], 3500), 0)
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        # We use automatic Y scaling, except we don't go past 0 and we are reversed (this is typical for IR).
        self.axes.autoscale(enable = True, axis = 'y')
        self.axes.set_ylim(self.axes.get_ylim()[1], 0)


class NMR_graph_maker_abc(Spectroscopy_graph_maker):
    """
    ABC for NMR spectra.
    """
            
    def __init__(self, output, graph, coupling = None, plot_labels = True, **kwargs):
        """        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        """
        # Call our parent.
        super().__init__(output, graph, **kwargs)
        
        # A dictionary of dicts. The outer key is the atom group 
        self.coupling = coupling if coupling is not None else {}
        
        # Axes labels.
        self.x_label = "Chemical Shift /ppm"
        self.y_label = "Arbitrary Units"
        #self.hide_y = True
        self.max_width = 1600
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = 15
        
        self.linewidth = 0.5
        self.cumulative_linewidth = 0.5
        self.column_linewidth = 0.5
        
        self.should_plot_labels = plot_labels
        
        self.x_padding = 0.05
        self.y_padding_percent = 1.5
        
        self.annotations = []
        self.focus_label = None
        
    def plot(self):
        """
        Plot the contents of our graph.
        """
        super().plot()
        if self.should_plot_labels:
            self.plot_labels()
            
    def get_label_for_peak(self, atom_group, graph, full = False):
        # Get the median of the x values (the middle point).
        x_coords, y_coords = self.transpose(graph.plot_cumulative_gaussian())
        
        x_coord = statistics.median(x_coords)
        y_coord = max(y_coords)
        
        # Now filter by a fraction of that amount.
        visible_x_values = [x for x, y in graph.plot_cumulative_gaussian() if y >= (y_coord * 0.95)]
        
        # Also get the height of the total curve at this point.
        total_x_coord, total_y_coords = self.transpose(((x,y) for x, y in self.graph.plot_cumulative_gaussian() if x >= min(visible_x_values) and x <= max(visible_x_values)))
        total_max_y = max(total_y_coords)
        
        # Get multiplicity of the peak.
        coupling = self.coupling.get(atom_group, {})
        mult = graph.multiplicity(atom_group, coupling)
        mult_string = "".join((multiplicity['symbol'] for multiplicity in mult))
        
        if full:
            label = r"$\mathdefault{{{{{}}}_{{{}}}}}$ ({})".format(atom_group.element, atom_group.index, mult_string)
            #label += "\n" + r"$\int$ = {}{}".format(len(atom_group.atoms), atom_group.element.symbol)
            #label += "\n{:.2f} ppm".format(x_coord) + r", $\int$ = {}{}".format(len(atom_group.atoms), atom_group.element.symbol)
            label += "\n{:.2f} ppm".format(x_coord) + r", ∫ = {}{}".format(len(atom_group.atoms), atom_group.element.symbol)
        
        else:
            label = r"$\mathdefault{{{{{}}}_{{{}}}}}$".format(atom_group.element, atom_group.index)
        
        return label, x_coord, total_max_y
        #return label, x_coord, y_coord
            
    def plot_label_for_peak(self, label, x_coord, y_coord):
        """
        Plot a label for a given sub-graph.
        """
#         return self.axes.annotate(
#             "{}".format(label),
#             (x_coord, y_coord),
#             textcoords = "offset pixels",
#             xytext = (0, 8),
#             horizontalalignment = "center",
#             fontsize = 8
#         )
        text = self.axes.text(
            # Pop the label 5% above the peak.
            x_coord, y_coord * 1.1,
            label,
            horizontalalignment = "center",
            fontsize = 10
        )
        text.set_bbox({"facecolor": 'white', "alpha": 0.75})
        return text
        
    def plot_labels(self):
        """
        Plot labels for each peak.
        """
        self.annotations = []
        for atom_group, graph in self.graph.graphs.items():
            self.annotations.append(self.plot_label_for_peak(*self.get_label_for_peak(atom_group, graph)))
            
    def make_graph(self):
        figure = super().make_graph()
        
                
        x, y = self.transpose(self.graph.peaks())
        # Decrease the value of all y (because adjusttext will repel labels away if they are too close)
        y = [val * 0.8 for val in y]
        
        #self.annotations = []
        
        if len(self.annotations) > 0:
            adjustText.adjust_text(
                self.annotations,
                x = x,
                y = y,
                objects = [self.focus_label] if self.focus_label is not None else None,
                ax = self.axes,
                avoid_self = False,
                expand_axes = False,
                expand = (1.5, 1.3),
                arrowprops = {"arrowstyle": '-', "linewidth": 0.5},
                min_arrow_len = 1,
                # Don't allow moving down on the y axis.
                #only_move = "xy+"
                only_move = "xy",
                #only_move = "y+",
                time_lim = 10
            )
        return figure


class NMR_graph_maker(NMR_graph_maker_abc):
    """
    A graph for displaying NMR spectra
    """
    
    def __init__(self, output, graph, coupling = None, **kwargs):
        """        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        """
        # Call our parent.
        super().__init__(output, graph, coupling = coupling, **kwargs)
        
        self.x_padding = None
        self.x_padding_percent = 0.05
        
        # NMR graphs have a range of different units depending on the atom and isotope being studied.
        # For example, 1H NMR typically range from ~10 -> 0
        # 13C typically range from 300-200 -> 0 etc.
        # To account for this, the 'scale' (in pixels) of the x-axis will be adjusted depending on the
        # range of values we have.
        # This is in inches.
        self.target_width = 10
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = None
        
        self.linewidth = 0.5
        self.cumulative_linewidth = 0.5
        self.column_linewidth = 0.5
    
    @classmethod
    def from_options(self, output, graph, *, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """
        return self(
            output,
            graph,
            enable_rendering = options['nmr']['enable_rendering'],
            plot_bars = options['nmr']['plot_bars'],
            **kwargs
        )
        
    def get_zoomed(self, output, focus, options, **kwargs):
        """
        Return a zoomed spectrum (a NMR_graph_zoom_maker object) focusing on a single peak.
        
        :param focus: The name of a peak.
        """
        return NMR_graph_zoom_maker.from_options(output, self.graph, coupling = self.coupling, focus = focus, options = options, **kwargs)
    
    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        # We need to get a list of all peaks that are above our cutoff point.
        # First determine our highest point.
        highest_point = max(self.transpose(self.graph.coordinates)[1])
        
        # Now filter by a fraction of that amount.
        visible_x_values = [x for x, y in self.graph.plot_cumulative_gaussian() if y >= (highest_point * self.peak_cutoff)]
        
        # Now we can just get our min-max.
        # NMR is typically shown on a relative scale, meaning that negative values are absolutely possible.
        # We always want to make sure that zero is shown however.
        #
        # NMR is also typically shown on a reversed scale.

        # First, work out how wide our graph will be.
        min_x = min(min(visible_x_values), 0)
        max_x = max(max(visible_x_values), 0)
              
        x_padding = (
            max_x - min_x
        ) * self.x_padding_percent
        
        # If we have no negative shifts, set zero as the end of one scale.
        if min(visible_x_values) < 0:
            minimum = min(visible_x_values) - x_padding
        
        else:
            # Don't extend beyond zero if there are no negative values.
            minimum = min(visible_x_values) - x_padding
            #minimum = 0
            
        # If we have no positive shifts, set zero as the end of one scale.
        if max(visible_x_values) > 0:
            maximum = max(visible_x_values) + x_padding
        
        else:
            maximum = 0
            
        # Space to allocate per unit of the y axis.
        self.inch_per_x = self.target_width / (maximum - minimum) 
        
        return self.axes.set_xlim(maximum, minimum)
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        super().auto_y_limits()
        self.axes.set_ylim(0, self.axes.get_ylim()[1] * self.y_padding_percent)


class NMR_graph_zoom_maker(NMR_graph_maker_abc):
    """
    A class for generating a zoomed NMR spectrum of a single (split) peak.
    """
    
    def __init__(self, output, graph, focus, coupling = None, plot_background_peaks = False, **kwargs):
        """
        Constructor for frequency graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a digilib.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        :param focus: An atom group to focus on.
        """
        # Call our parent.
        super().__init__(
            output,
            graph,
            coupling = coupling,
            plot_peaks = True,
            plot_cumulative_peak = True,
            x_limits = 'auto',
            peak_cutoff = 0.02,
            **kwargs)
        
        self.plot_background_peaks = plot_background_peaks
        self.focus = focus
        
        self.x_padding = None
        self.x_padding_percent = 0.5
        
        self.target_width = 3.5
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = None
        
        self.cumulative_linewidth = 1
        
        if self.focus not in self.graph.graphs:
            raise ValueError("The sub-graph '{}' is not recognised".format(self.focus))
        
    @classmethod
    def from_options(self, output, graph, *, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """
        return self(
            output,
            graph,
            enable_rendering = options['nmr']['enable_rendering'],
            plot_bars = options['nmr']['plot_zoom_bars'],
            **kwargs
        )
        
    def visible_window(self):
        """
        Determine the x-axis window that is to be displayed.
        """
        # Get the graph object corresponding to the peak we are interested in.
        graph = self.graph.graphs[self.focus]
        
        # We need to get a list of all peaks that are above our cutoff point.
        # First determine our highest point.
        coords = graph.plot_cumulative_gaussian()
        highest_point = max(self.transpose(coords)[1])
        
        # Now filter by a fraction of that amount.
        visible_x_values = [x for x, y in coords if y >= (highest_point * self.peak_cutoff)]
        
        x_padding = (
            max(visible_x_values) -min(visible_x_values)
        ) * self.x_padding_percent
        
        maximum = max(visible_x_values) + x_padding
        minimum = min(visible_x_values) - x_padding
        
        return (minimum, maximum)
        
    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        minimum, maximum = self.visible_window()
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = self.target_width / (maximum - minimum) 
        
        # Now we can just get our min-max.
        self.axes.set_xlim(
            maximum,
            minimum
        )
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        # Get the graph object corresponding to the peak we are interested in.
        graph = self.graph.graphs[self.focus]
        focus_spectrum = self.transpose(graph.plot_cumulative_gaussian())
        
        # Get all peaks.
        minimum, maximum = self.visible_window()
        spectrum = self.transpose([(x,y) for x, y in self.graph.plot_cumulative_gaussian() if x >= minimum and x <= maximum])
        
#         # Cut the full spectrum within the limits of our focus graph.
#         # NOTE: This works well, but as the plotted graph is normally wider than
#         # just the values of the peak of interest (because of x_padding), there
#         # may still be peaks in view which are cut off...
#         cut_start = spectrum[0].index(focus_spectrum[0][0])
#         cut_end = spectrum[0].index(focus_spectrum[0][-1])
#         
#         cut_spectrum = (spectrum[0][cut_start:cut_end], spectrum[1][cut_start:cut_end])
        
        
        
        highest_point = max(spectrum[1])

        height = highest_point * 1.3
        
        # Clamp to 0 -> pos.
        self.axes.set_ylim(0 - height * 0.025, height)
        
    def plot_lines(self):
        """
        Plot x and y values as individual lines on the graph we are building.
        """
        if self.plot_background_peaks:
            for graph_name, graph in self.graph.graphs.items():
                if graph_name != self.focus:
                    plot = graph.plot_cumulative_gaussian()
                    data = self.transpose(plot)
                    self.axes.plot(*data, 'C0--', linewidth = self.linewidth)
                
        # Plot our focus last (so it appears on top).
        plot = self.graph.graphs[self.focus].plot_cumulative_gaussian()
        data = self.transpose(plot)
        self.axes.plot(*data, 'r-', linewidth = self.cumulative_linewidth)
        
    def plot_columns(self):
        """
        Plot vertical columns for each plot on the graph we are building.
        """
        for atom_group, graph in self.graph.graphs.items():
            if atom_group != self.focus:
                columns = [[(x, 0) , (x, y)] for x, y in graph.coordinates]
                colors = (0,0,0)
                self.axes.add_collection(matplotlib.collections.LineCollection(columns, linewidths = self.column_linewidth, colors = colors))
                
        graph = self.graph.graphs[self.focus]
        columns = [[(x, 0) , (x, y)] for x, y in graph.coordinates]
        colors = (1,0,0)
        self.axes.add_collection(matplotlib.collections.LineCollection(columns, linewidths = self.column_linewidth, colors = colors))
        
    def plot_labels(self):
        """
        Plot labels for each peak.
        """
        self.annotations = []
        minimum, maximum = self.visible_window()
        
        for atom_group, graph in self.graph.graphs.items():
            if atom_group != self.focus:
                # Only plot the label if the peak is within our visible window.
                x_coords, y_coords = self.transpose(graph.plot_cumulative_gaussian())
                x_coord = statistics.median(x_coords)
                if x_coord > minimum and x_coord < maximum:
                    self.annotations.append(
                        self.plot_label_for_peak(*self.get_label_for_peak(atom_group, graph))
                    )
                
        # Plot a more complete label for our focus peak.
        label, x_coord, y_coord = self.get_label_for_peak(self.focus, self.graph.graphs[self.focus], full = True)
        
        # Get couplings for our main atom group.
        coupling = self.coupling.get(self.focus, {})
        
        # Get the multiplicities of our peak.
        mult = self.graph.graphs[self.focus].multiplicity(self.focus, coupling)
        couplings = unpack_coupling(coupling)
        
        # Add coupling info if we have it.
#         couplings = self.coupling.get(self.focus, {})
#         couplings = {
#             (coupling_group, coupling_isotope): isotope_coupling for coupling_group, atom_dict in couplings.items() for coupling_isotope, isotope_coupling in atom_dict.items()
#             if isotope_coupling.groups[isotope_coupling.other(atom_group)].element[isotope_coupling.isotopes[isotope_coupling.other(atom_group)]].abundance / 100 > satellite_threshold
#         }
#         # Sort couplings.
#         couplings = dict(sorted(couplings.items(), key = lambda coupling: abs(coupling[1].total), reverse = True))
        
        # Add coupling info.
        if len(couplings) > 0 and mult[0]["number"] != 1:
            # Only show couplings for peaks we can actually distinguish.
            for (coupling_group, coupling_isotope), coupling in list(couplings.items())[:len(mult)]:
                label += "\n" + r"$\mathdefault{{^{{{}}}}}$J = {:.2f} Hz ($\mathdefault{{^{{{}}}{}_{{{}}}}}$, {}{})".format(
                    coupling.distance if coupling.distance is not None else "",
                    coupling.total,
                    coupling_isotope,
                    coupling_group.element,
                    coupling_group.index,
                    #len(coupling_group.atoms),
                    coupling.num_coupled_atoms(self.focus),
                    coupling_group.element
                )
        #self.focus_label = self.plot_label_for_peak(label, x_coord, y_coord)
        self.annotations.append(self.plot_label_for_peak(label, x_coord, y_coord))
