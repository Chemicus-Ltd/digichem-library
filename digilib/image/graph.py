from matplotlib import ticker
from matplotlib.transforms import Bbox
import matplotlib.figure
import math

from digilib.image import Image_maker

class Graph_image_maker(Image_maker):
    """
    A class used for creating graph images.
    
    Most of the heavy lifting here is achieved by the matplotlib library.
    """
    
    def __init__(self, *args, max_width = None, max_height = None, **kwargs):
        """
        Constructor for Graph_image_maker objects.
        
        :param max_width: The maximum image width in pixels, None for no limit.
        :param max_height: The maximum image height in pixels, None for no limit.
        """
        super().__init__(*args, **kwargs)
        
        # Set our output dpi (higher values = bigger image/quality).
        self.output_dpi = 130
        
        # The width and height (in inch) that we initially create our figure with. Note that the true figure size may change from this value.
        self.init_image_width = 5
        self.init_image_height = 5
        
        # If true, the label and tick labels will be hidden for that axis.
        self.hide_x = False
        self.hide_y = False
        
        # The label to use for the x and y axes.
        self.x_label = ""
        self.y_label = ""
        
        # The matplotlib figure object that represents our graph.
        self.figure = None
        
        # The plotting area within our figure.
        self.axes = None
        
        # Image clamping.
        self.max_width = max_width
        self.max_height = max_height
        
    def make_files(self):
        """
        Make the image described by this object.
        """
        # Change font and other defaults.
        #pyplot.rcParams.update({'font.size': 14, 'font.family': 'DejaVu Sans', 'axes.linewidth': 2})
        matplotlib.rcParams.update({'font.size': 14, 'font.family': 'DejaVu Sans', 'axes.linewidth': 2})
                
        # Make our (empty) plot.
        # pyplot uses a nasty interface, will switch to OO at some point.
        #pyplot.figure(figsize = (self.init_image_width, self.init_image_height), tight_layout = True)
        #pyplot.figure(figsize = (self.init_image_width, self.init_image_height), tight_layout = False)
        
        # Plot our graph.
        self.make_graph()
        
        # And save to file.
        self.figure.savefig(self.output, dpi = self.output_dpi)
        #pyplot.savefig(self.output, dpi = self.output_dpi)
        
    def make_graph(self):
        """
        Make the graph data described by this object.
        
        The matplotlib figure is available at self.figure (and is also returned by this method for convenience). The graph is not saved to file.
        
        :return: A reference to the created matplotlib Figure object.
        """
        # First make our figure object, which is the top-level matplotlib container.
        self.figure = matplotlib.figure.Figure(figsize = (self.init_image_width, self.init_image_height))
        
        # Also get our main axes (all the graph types we support so far contain a single plot/axes).
        self.axes = self.figure.add_subplot(111)
        
        # Call plot(), which is defined by sub-classes to actually draw the graph.
        self.plot()
        
        # Call adjust_axes(), which does what it suggests (also adds axis labels etc.)
        self.adjust_axes()
        
        # Finally, clamp our image size if we're going to be too big.
        if self.max_width is not None:
            self.clamp_dimension(0, self.max_width)
        if self.max_height is not None:
            self.clamp_dimension(1, self.max_height)
        
        # Return a reference to the figure (for convenience).
        return self.figure
        
    def plot(self):
        """
        Main workhorse of the make_files method. Inheriting classes should write their own implementation.
        
        This is called automatically by make_graph().
        """
        pass
    
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """ 
        if self.hide_x:
            self.x_label = "Arbitrary units"
            self.axes.xaxis.set_ticks([])
        if self.hide_y:
            self.y_label = "Arbitrary units"
            self.axes.yaxis.set_ticks([])
        
        # Add our axes labels.
        if not self.x_label == "": 
            self.axes.set_xlabel(self.x_label, labelpad = 8, fontsize = 16, weight = 'bold')
        if not self.y_label == "":
            self.axes.set_ylabel(self.y_label, labelpad = 8, fontsize = 16, weight = 'bold')
    
    def clamp_dimension(self, dim, maximum):
        """
        Ensure the image does not exceed a certain number of pixels in a given dimension.
        
        :param dim: The dimension to clamp; 0 for x, 1 for y.
        :param maximum: The maximum allowed pixels.
        """
        # Check to see if we've exceeded our max size.
        fig_size = self.figure.get_size_inches()
        if (fig_size[dim] * self.output_dpi) > maximum:
            # We're too big, reduce to max.
            fig_size[dim] = maximum / self.output_dpi
            self.figure.set_size_inches(fig_size)
                        
            # Update layout. Calling this twice is necessary for some reason loool
            self.figure.tight_layout()
            self.figure.tight_layout()
    
    def constant_scale(self, axis, inch_per_axis):
        """
        Adjust the size of our figure so a given axis has a constant scale.
        
        You should adjust the corresponding axis' limits before calling this method.
        
        :param axis: The axis to adjust, 0 for the x axis, 1 for the y axis.
        :param inch_per_axis: The number of inches per unit of the given axis to adjust to. 
        """
        
        # First determine how much space we need for our margins.
        axes_position = self.axes.get_position(False).get_points()
        figure_size = self.figure.get_size_inches()
        # Axes_position is in fractions, so multiply by our current figsize.
        start_margin = axes_position[0][axis] * figure_size[axis]
        end_margin = (1 - axes_position[1][axis]) * figure_size[axis]
        
        # Next determine how much space we need for our actual data.
        axes_lim = (self.axes.get_xlim(), self.axes.get_ylim())
        axis_min = axes_lim[axis][0]
        axis_max = axes_lim[axis][1]
        #axis_min = pyplot.axis()[axis *2]
        #axis_max = pyplot.axis()[axis *2 +1]
        
        # We take the absolute because min can actually be bigger than max if our axes are inverted.
        axis_length = math.fabs(axis_max - axis_min)
        axis_size = axis_length * inch_per_axis
        
        # Update size.
        figure_size[axis] = start_margin + end_margin + axis_size
        self.figure.set_size_inches(figure_size)
        # And update the position of our axes.
        axes_position[0][axis] = start_margin / figure_size[axis]
        axes_position[1][axis] = 1 - (end_margin / figure_size[axis])
        self.axes.set_position(Bbox(axes_position))
        
    
class OneD_graph_image_maker(Graph_image_maker):
    """
    Classes for creating '1D' graphs using matplotlib's event plot.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Options that we'll pass to matplotlib's event plot function.
        self.plot_options = {
            "orientation": "vertical",
            "lineoffsets": 0.75,
            "linelengths": 1,
            "linewidths": None,
            "colors": None
        }
        
        # Set our output width and height.
        self.init_image_width = 4.1
        self.init_image_height = 5.6
        
    def plot(self):
        """
        Make the image described by this object.
        """
        # Plot our 1D data.
        self.plot_lines()
        
        # Annotate with data labels.
        self.plot_labels()
            
    def plot_lines(self):
        """
        Plot the lines that make up our graph.
        
        This default implementation does nothing, inheriting classes should write their own.
        """
        pass
    
    def plot_labels(self):
        """
        Plot annotation labels that explain the data in out graph.
        
        This default implementation does nothing, inheriting classes should write their own.
        """
        pass
        
class Convergence_graph_maker(Graph_image_maker):
    """
    A class for creating graphs that show the convergence of energy in optimisation calculations.
    """
    
    def __init__(self, output, energies, **kwargs):
        """
        Constructor for graph makers.
        
        :param energies: A list of energies to plot (probably an Energy_list object).
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        """
        super().__init__(output, **kwargs)
        self.energies = energies
        
        # Our output width and height (in inch).
        self.init_image_width = 7
        self.init_image_height = 5.265
        
        # Axis titles.
        self.x_label = "Step number"
        self.y_label = "Energy /eV"
        
    @classmethod
    def from_options(self, output, *, energies, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            energies = energies,
            **kwargs
        )
    
    def plot(self):
        """
        Plot the contents of our graph.
        """
        # Plot a simple 'scatter' style graph.
        self.axes.plot(range(1, len(self.energies) +1), self.energies)    
            
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """
        super().adjust_axes()
        
        # Axes range, the format is [xmin, xmax, ymin, ymax].
        # We'll add a bit on to the y axis so we can see better.
        y_diff = max(self.energies) - min(self.energies)
        y_fudge = y_diff * 0.05
        # Also the x.
        x_diff = len(self.energies) - 1
        x_fudge = x_diff * 0.05
        
        self.axes.set_xlim(1 - x_fudge, len(self.energies) + x_fudge)
        self.axes.set_ylim(min(self.energies) - y_fudge, max(self.energies) + y_fudge)
        
        # Change spacing.
        self.axes.xaxis.set_major_locator(ticker.MaxNLocator(nbins = '12', steps = [1, 2, 2.5, 5, 10], integer = True))
        
        # Use matplotlib to automatically layout our graph.
        self.figure.tight_layout()
        
        