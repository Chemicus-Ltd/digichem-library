from matplotlib import ticker
from matplotlib.ticker import FuncFormatter
from itertools import chain
import math

from digilib.image.graph import OneD_graph_image_maker
#from digilib.result.excited_state import Excited_state, Excited_state_list
import digilib.result.excited_state

class Excited_states_diagram_maker(OneD_graph_image_maker):
    """
    Class for making orbital energy diagrams.
    """
    
    def __init__(self, output, excited_states, ground_state, y_limits = "auto", show_dest = True, **kwargs):
        super().__init__(output, **kwargs)
        self.excited_states = excited_states
        self.ground_state = ground_state
                        
        # Each multiplicity will be plotted as a separate column, so first we need to organise our data.
        states = chain(self.excited_states, [self.ground_state]) if self.ground_state is not None else self.excited_states
        self.grouped_states = (type(self.excited_states)(sorted(states, key = lambda state: state.energy))).group()
        # We'll also keep hold of just our excited states (so no ground).
        self.grouped_excited_states = self.excited_states.group()
        
        # Whether to display our dEst label, this is forced off if we don't have only singlets and triplets.
        self.show_dest = show_dest if 1 in self.grouped_excited_states and 3 in self.grouped_excited_states and len(self.grouped_excited_states) == 2 else False
        self.left_padding = 0
        self.right_padding = 0    
            
        # Set our output dpi (higher values = bigger image/quality).
        self.output_dpi = 140
        
        # Initial image dimensions (we'll change these to suit).
        self.init_image_width = 5.6
        self.init_image_height = 5.6
        
        # Amount of space to allocate per unit of each axis.
        self.inch_per_x = 1.66
        #self.inch_per_y = 0.94
        self.inch_per_y = 1.1
        
        # Y axis label (we don't have an x axis label).
        self.y_label = "Energy /eV"
        
        self.plot_options = {
            "orientation": "vertical",
            "lineoffsets": 1,
            "linelengths": 0.9,
            "linewidths": None,
            "colors": None
        }
        
        # Save our scaling method.
        self.limits_method = y_limits
    
    @classmethod
    def from_options(self, output, *, excited_states, ground_state, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            excited_states = excited_states,
            ground_state = ground_state,
            show_dest = options['excited_states_diagram']['show_dest'],
            y_limits = options['excited_states_diagram']['y_limits'],
            enable_rendering = options['excited_states_diagram']['enable_rendering'],
            **kwargs
        )
            
    def plot_lines(self):
        """
        Plot the lines that make up the body of the graph.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        # Turn our dictionary of grouped excited states into an ordered list for matplotlib.
        data = [
            [excited_state.energy for excited_state in self.grouped_states[grouped_key]]
            #self.grouped_states[grouped_key] 
            for grouped_key 
            in sorted(self.grouped_states)
        ]
            
        # We can use an event plot which is good for 1D graphs like this.
        self.axes.eventplot(data, linestyles = "solid", **self.plot_options)
    
    def get_x_pos_from_multiplicity(self, multiplicity):
        """
        Get the X coordinate of an energy state.
        
        :raises ValueError: If the given multiplicity is not plotted. 
        :param multiplicity: The multiplicity of the energy state, which determines its x position.
        :return: The X coord.
        """
        # Get all the multiplicities that we're using.
        mult = sorted(self.grouped_states)
        
        # Our X position depends on the multiplicity of our GS and how many total multiplicities there are (as well as how spaced out each column on the graph is).
        if len(self.grouped_states) > 1:
            return mult.index(multiplicity) * self.plot_options['lineoffsets']
        else:
            # It appears matplotlib uses a different numbering if we only have 1 column, where x == 1 (normally the first column is x == 0). This makes me sad. Might be bug...
            return 1
    
    def annotate_ground_state(self):
        """
        Add a label to our ground state.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        # First work out where our GS is on the graph.
        x_pos = self.get_x_pos_from_multiplicity(self.ground_state.multiplicity)
        
        # Now plot our text.
        self.axes.annotate(
             r"$\bf{{{}}}_{{{}}}$".format(self.ground_state.multiplicity_symbol, self.ground_state.multiplicity_level),
             (x_pos, self.ground_state.energy),
             textcoords = "offset pixels",
             xytext = (0, 7),
             horizontalalignment = "center",
             fontsize = 12
         )
        
    def annotate_excited_state(self, state):
        """
        Add a label to an excited state
        
        This method is called automatically as part of the make() method.
        :param multiplicity_symbol: The symbol (a string) to annotate (eg, S1, T1).
        :return: Nothing.
        """
        # First work out where our state is on the graph.
        #state = self.excited_states.get_state(multiplicity_symbol)
        x_pos = self.get_x_pos_from_multiplicity(state.multiplicity)
        
        # Assemble our label.
        label_string =     r"$\bf{{{}}}_{{{}}}$: ".format(state.multiplicity_symbol, state.multiplicity_level) + "{:0.2f} eV".format(state.energy)
                        
        # Add oscillator strength (if we have it).
        if state.oscillator_strength is not None:
            label_string +=    "\n" + "f: {:0.2f}".format(state.oscillator_strength)
        
        # Now plot our text.
        self.axes.annotate(
             label_string,
             (x_pos, state.energy),
             textcoords = "offset pixels",
             xytext = (0, -7),
             verticalalignment = "top",
             horizontalalignment = "center",
             fontsize = 12
         )
    
    def annotate_dest(self, x, y):
        """
        """
        # First the text label.
        self.axes.annotate(
                r"$\bf{{ΔE_{{ST}}}}$: {:0.2f} eV".format(self.excited_states.singlet_triplet_energy),
                (x, y),
                textcoords = "offset pixels",
                xytext = (0,0),
                verticalalignment = "center",
                horizontalalignment = "left",
                fontsize = 12
            )
    
    def plot_labels(self):
        """
        Plot the labels on our graph that identify various excited states and the singlet/triplet splitting energy (dEst) if available.
        
        This method is called automatically as part of the make() method.
        :return: Nothing.
        """
        # First add a label for our ground state, unless there isn't room (this hack again).
        lowest_state = self.excited_states[0]
        # We need about 0.6 eV space.
        if (lowest_state.energy - 0.6) >= 0:
            self.annotate_ground_state()
        
        # Now we'll add a label for the lowest state of each multiplicity.
        for multiplicity in self.grouped_excited_states:
            self.annotate_excited_state(self.grouped_excited_states[multiplicity][0])
            
            
        # We'll also try and add a dEst label (obviously this is only possible if there is an S1 and T1 available).
        if self.show_dest:
            # First get our two states.
            S1 = self.excited_states.get_state("S(1)")
            T1 = self.excited_states.get_state("T(1)")
            
            # Find the midpoint between them.
            y_pos = min([S1.energy, T1.energy]) +  (abs(S1.energy - T1.energy) /2)
            
            # The label will be to the right of S1/T1.
            x_pos = max(self.get_x_pos_from_multiplicity(S1.multiplicity), self.get_x_pos_from_multiplicity(T1.multiplicity)) + (self.plot_options['linelengths'] /2) + 0.2
            # Add the label
            self.annotate_dest(x_pos, y_pos)
            
            # Now draw our vertical arrow.
            self.axes.annotate(
                "",
                xy = (x_pos-0.05, S1.energy),
                xytext = (x_pos-0.05, T1.energy),
                arrowprops=dict(arrowstyle="<->")
            )
            
            # Now extend a line from S1 and T1 to meet the vertical arrow.
            for state in (S1, T1):
                self.axes.annotate(
                    "",
                    xy = (self.get_x_pos_from_multiplicity(state.multiplicity) + (self.plot_options['linelengths'] /2), state.energy),
                    xytext = (x_pos-0.05, state.energy),
                    arrowprops=dict(arrowstyle="-", linestyle="--")
                )
            
            # Pad to the right so our dest label is chopped off.
            self.right_padding = 0.7            

        
    def x_axis_formatter(self, x, pos):
        """
        This method tells matplotlib how to label our X axis (with multiplicity labels).
        """
        # The multiplicities that we're using.
        mult = sorted(self.grouped_states)
        
        if len(self.grouped_states) > 1:
            # Check to see if our x is in our list of multiplicities.
            if x % self.plot_options['lineoffsets'] == 0 and x >= 0:
                try:
                    return digilib.result.excited_state.Excited_state.multiplicity_number_to_string((mult[int(x / self.plot_options['lineoffsets'])]))
                    
                except IndexError:
                    # Out of range.
                    return ""
        else:
            # Matplotlib uses different counting when there's only 1 column :(.
            if x == 1:
                return digilib.result.excited_state.Excited_state.multiplicity_number_to_string(mult[0])
            else:
                return ""
            
    def all_limits(self):
        """
        Limit the Y axis so all states are visible.
        """
        # Use auto scaling.
        self.axes.autoscale(enable = True, axis = 'y')
        
        # We'll go a small amount below 0 so our S0 is off the bottom of the graph.
        min_y = -0.1
        # Use a sneaky hack so that if the max energy is close to an integer, we'll use the next highest integer.
        #max_y = math.ceil(self.excited_states[-1].energy +0.1)
        max_y = self.excited_states[-1].energy + 0.1
        self.axes.set_ylim(min_y, max_y)
        
        # Adjusty hack.
        self._adjust_limits_for_zero()
        
    def auto_limits(self):
        """
        Limit the Y axis so S1, T1 ... N1 etc are all visible.
        """
        self.all_limits()
        
        # We'll show the lowest excited state of each mult.
        highest_state = max([self.grouped_excited_states[multiplicity][0].energy for multiplicity in self.grouped_excited_states])
        # Use a sneaky hack so that if the max energy is close to an integer, we'll use the next highest integer.
        max_y = math.ceil(highest_state +0.1)
        self.axes.set_ylim(None, max_y)
        
        
    def _adjust_limits_for_zero(self):
        """
        This method is a hack which lowers the y axis limit in cases where S1/T1 is very close to S0
        
        This is necessary because otherwise the state annotations will overflow the bottom of the diagram.
        """
        # First get our lowest E excited state.
        lowest_state = self.excited_states[0]
        # We need about 0.7 eV space.
        if (lowest_state.energy - 0.5) < 0:
            # Add on a bit extra.
            self.axes.set_ylim(bottom = lowest_state.energy - (0.5))
        
        
    def simple_limits(self, y_min = -0.1, y_max = 5):
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
                
        # Adjust axis
        # Add a bit o' padding (often we actually reduce padding from the default).
        xlim = self.axes.get_xlim()
        self.axes.set_xlim(xlim[0] - (self.left_padding -0.25), xlim[1] + (self.right_padding -0.25))
        #self.axes.set_xlim(xlim[0] +0.25, xlim[1] -0.25)
        
        # Now decide how big our y axis will be.
        if self.limits_method == "all":
            self.all_limits()
        elif self.limits_method == "auto":
            self.auto_limits()
        elif isinstance(self.limits_method, tuple) or isinstance(self.limits_method, list):
            self.simple_limits(self.limits_method[0], self.limits_method[1])
        else:
            raise ValueError("Unknown limits method '{}'".format(self.limits_method))
        
        # Set custom text on our x axis.
        self.axes.xaxis.set_major_formatter(FuncFormatter(self.x_axis_formatter))
        #self.axes.xaxis.set_major_locator(ticker.MultipleLocator(1))
        self.axes.xaxis.set_major_locator(ticker.FixedLocator([self.get_x_pos_from_multiplicity(multiplicity) for multiplicity in self.grouped_excited_states]))
        
        # Both our dimensions can scale to fit more data.
        self.figure.tight_layout()
        self.constant_scale(0, self.inch_per_x)
        self.constant_scale(1, self.inch_per_y)

