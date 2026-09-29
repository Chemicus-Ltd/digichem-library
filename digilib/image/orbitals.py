from matplotlib import ticker
import statistics
import math

from digilib.exception.base import Result_unavailable_error
from digilib.image.graph import OneD_graph_image_maker


class Orbital_diagram_maker(OneD_graph_image_maker):
    """
    Class for making orbital energy diagrams.
    """
        
    def __init__(self, output, orbitals, full_axis_lines = False, y_limits = "auto", limits_fallback = True, **kwargs):
        """
        Constructor for graph makers.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg: png, jpeg).
        :param orbitals: A list of orbitals to plot.
        :param full_axis_lines: If False, a boundary line is drawn for the x axis but nowhere else. If True, boundary lines are drawn on all 4 sides of the plotting area.
        :param y_limits: String controlling how the y axis limits are set. Options are 'all' for showing all orbitals, 'auto' for standard auto scaling showing the HOMO, LUMO and 0 point, or 'center' for placing the HOMO-LUMO gap in the center of the y axis. Alternatively, limits_method can be a tuple of (y_min, y_max) which will be used directly as axis limits.
        :param limits_fallback: If True, use simple_limits() if the chosen limits_method fails. If false, an exception is raised.
        """
        super().__init__(output, **kwargs)
        self.orbitals = orbitals
        
        # Set our output width and height.
        self.init_image_width = 4.1
        # We use a wider graph if we are using unrestricted orbitals
        if self.orbitals.spin_type != "none":
            self.init_image_width = 4.9
        self.init_image_height = 5.6
        self.output_dpi = 130
        
        # Axes label (no x label here).
        self.y_label = 'Energy /eV'
        # Y-axis scaling.
        self.inch_per_y = 0.62
        
        # Save options.
        self.full_axis_lines = full_axis_lines
        self.limits_method = y_limits
        self.limits_fallback = limits_fallback
        
    @classmethod
    def from_options(self, output, *, orbitals, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            orbitals = orbitals,
            full_axis_lines = options['orbital_diagram']['full_axis_lines'],
            y_limits = options['orbital_diagram']['y_limits'],
            enable_rendering = options['orbital_diagram']['enable_rendering'],
            **kwargs
        )
        
    def plot_lines(self):
        """
        Plot the lines that make up the body of the graph.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        # We'll plot occupied and virtual orbitals separately, with different line styles.
        occupied = [mo.energy for mo in self.orbitals if mo.HOMO_difference <= 0]
        virtual = [mo.energy for mo in self.orbitals if mo.HOMO_difference > 0]
        
        # We can use an event plot which is good for 1D graphs like this.
        self.axes.eventplot(occupied, linestyles = "solid", **self.plot_options)
        self.axes.eventplot(virtual, linestyles = "dashed", **self.plot_options)
        
    def plot_labels(self):
        """
        Plot the labels on our graph that identify the HOMO, LUMO and energy gap.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        # First work our what our x coord is (this is the same for all points).
        x_pos = self.plot_options["lineoffsets"]
        start_x_pos = self.plot_options["lineoffsets"] - self.plot_options["linelengths"] /2
        
        # Try just plotting labels for HOMO and LUMO (because we're pretty sure there should be space).
        try:
            HOMO = self.orbitals.get_orbital(HOMO_difference = 0)
        except Exception:
            # Couldn't get HOMO.
            HOMO = None
            
        try:
            LUMO = self.orbitals.get_orbital(HOMO_difference = 1)
        except Exception:
            # Couldn't get LUMO.
            LUMO = None
        
        if HOMO is not None:
            # HOMO text goes above the line.
            self.axes.annotate(
                "{}: {:0.2f} eV".format(HOMO.label, HOMO.energy),
                (x_pos, HOMO.energy),
                textcoords = "offset pixels",
                xytext = (0, 7),
                horizontalalignment = "center"
            )

        if LUMO is not None:    
            # LUMO text goes below the line.
            self.axes.annotate(
                "{}: {:0.2f} eV".format(LUMO.label, LUMO.energy),
                (x_pos, LUMO.energy),
                textcoords = "offset pixels",
                xytext = (0, -7),
                verticalalignment = "top",
                horizontalalignment = "center"
            )
        
        if HOMO is not None and LUMO is not None:
            # Now plot our dE HOMO-LUMO arrow.
            self.axes.annotate(
                "",
                xy = (start_x_pos, HOMO.energy),
                xytext = (start_x_pos, LUMO.energy),
                arrowprops=dict(arrowstyle="<->")
            )
            
            # And then our dE HOMO-LUMO energy.
            self.axes.annotate(
                "ΔE: {:0.2f} eV".format(LUMO.energy - HOMO.energy),
                (x_pos, statistics.median([HOMO.energy, LUMO.energy])),
                verticalalignment = "center",
                horizontalalignment = "center"
            )
            
    def all_limits(self):
        """
        Limit the Y axis so all orbitals are visible.
        """
        self.axes.autoscale(enable = True, axis = 'y')
    
    def center_limits(self):
        """
        Limit the Y axis so the HOMO and LUMO are in the center of the diagram, and the 0 energy point is visible.
        
        :raises Result_unavailable_error: If either the HOMO or LUMO is not available.
        """
        # Get our FMOs.
        HOMO = self.orbitals.get_orbital(HOMO_difference = 0)
        LUMO = self.orbitals.get_orbital(HOMO_difference = 1)
        
        # Calculate the midpoint between the FMOs.
        midpoint = HOMO.energy - ((HOMO.energy - LUMO.energy) /2)
        
        # Work out y_max, which is set depending on the position of the LUMO.
        y_max = max(0, math.ceil(LUMO.energy +1.1))
        
        # Y min is then set so that midpoint is indeed in the midpoint.
        y_min = midpoint - (y_max - midpoint)
        
        #
        y_min = round(y_min)
        
        # Set our axis limits.
        self.axes.set_ylim(y_min, y_max)
        
    def standard_limits(self):
        """
        Limit the Y axis so that the HOMO, LUMO and 0 energy point are all visible.
        
        :raises Result_unavailable_error: If either the HOMO or LUMO is not available.
        """
        # Get our FMOs.
        HOMO = self.orbitals.get_orbital(HOMO_difference = 0)
        LUMO = self.orbitals.get_orbital(HOMO_difference = 1)
        
        # Add on a bit to both HOMO and LUMO so we expand if the orbital is close to the end of the axis.
        y_min = math.floor(HOMO.energy -1.1) if HOMO is not None else -8
        y_max = max(0, math.ceil(LUMO.energy +1.1)) if LUMO is not None else 0
        self.axes.set_ylim(y_min, y_max)
    
    def simple_limits(self, y_min = -8, y_max = 0):
        """
        Limit the Y axis between between two points.
        
        :param y_min: The start of the y axis.
        :param y_max: The end of the y axis.
        """
        self.axes.set_ylim(y_min, y_max)
    
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        super().adjust_axes()
        
        # Orbital diagrams only have the left axis spine to enable easy stacking.
        if not self.full_axis_lines:
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.spines['bottom'].set_visible(False)
                
        # Adjust axis
        # Our x axis is always the same.
        self.axes.set_xlim(0, 1.5)
        
        # Our Y axis changes with our given data. If we a HOMO and/or LUMO, we'll stretch appropriately to fit them in.
        # If We're missing one of the FMO (unlikely, but possible), we'll resort to semi-random defaults (for now).
        try:
            if self.limits_method == "all":
                self.all_limits()
            elif self.limits_method == "center":
                self.center_limits()
            elif self.limits_method == "auto":
                self.standard_limits()
            elif isinstance(self.limits_method, tuple) or isinstance(self.limits_method, list):
                self.simple_limits(self.limits_method[0], self.limits_method[1])
            else:
                raise ValueError("Unknown limits method '{}'".format(self.limits_method))
        except Result_unavailable_error:
            # HOMO-LUMO not available.
            if self.limits_fallback:
                self.simple_limits()
            else:
                raise
        
        # Change spacing of tick markers.
        # No tick marks for x axis.
        self.axes.set_xticks([])
        self.axes.yaxis.set_major_locator(ticker.MultipleLocator(1))
                
        # Our y axis scales with data.
        self.figure.tight_layout()
        self.constant_scale(1, self.inch_per_y)
        
        