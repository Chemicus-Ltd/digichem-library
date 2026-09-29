# Methodologies for determining molecule linearity.
import math
import numpy
from statistics import  mean
try:
    from bayes_opt import BayesianOptimization
    from bayes_opt.acquisition import UpperConfidenceBound
    from bayes_opt.exception import NotUniqueError

except ImportError:
    BayesianOptimization = None
    UpperConfidenceBound = None
    NotUniqueError = None

from digichem.result.alignment import Alignment, Axis_swapper_mix
from digichem.exception import Digichem_exception


class Grid_search(Alignment, Axis_swapper_mix):
    """
    The brute force method for determining molecule alignment. 
    """

    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Grid Search", "GRID"]

    def __init__(self, atoms, steps = 90, *args, charge = None, **kwargs):
        """
        Constructor for Brute_force.

        :param atoms: List of atom objects to orientate.
        :param charge: Optional overall electronic charge of the molecule.
        :param steps: The number of steps to search over for each axis. The total number of iterations will be steps*steps*steps.
        """
        self.steps = steps
        super().__init__(atoms, *args, charge=charge, **kwargs)

    @property
    def human_method_type(self):
        return "{} ({})".format(self.CLASS_HANDLE[0], self.steps)
    
    @property
    def method_type(self):
        return "{}{}".format(self.CLASS_HANDLE[1], self.steps)
        
    def align_axes(self):
        """
        Realign the axes of our coordinate system.
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # Begin by translating to centre of coordinates.
        coords = self.get_coordinate_list()
        self.translate((-mean(coords[0]) , -mean(coords[1]), -mean(coords[2])))

        # Now find the optimal rotation angles.
        x_theta, y_theta, z_theta = self.grid_search(self.steps)

        # Rotate by the requested amount.
        self.rotate_YZ(x_theta)
        self.rotate_XZ(y_theta)
        self.rotate_XY(z_theta)

        # Swap the axes 
        self.reassign_axes()

    def kernel(self, x_angle, y_angle, z_angle, points):
        new_points = [
            self.rotate_coords_XY(
                self.rotate_coords_XZ(
                    self.rotate_coords_YZ(coord, x_angle), y_angle), z_angle)
            for coord in points
        ]
        #x, y, z = self.get_coordinate_list(new_points)
        x, y, z = list(map(list, zip(* [(coord[0], coord[1], coord[2]) for coord in new_points])))

        X = max(x) - min(x)
        Y = max(y) - min(y)
        Z = max(z) - min(z)

        volume = X * Y * Z
        return volume

    def grid_search(
        self,
        steps,
        start = 0,
        stop = 0.5*math.pi,
        *args,
        x_start = None,
        x_stop = None,
        y_start = None,
        y_stop = None,
        z_start = None,
        z_stop = None,
        memory = None,
    ):
        """
        Perform a grid search over the three rotation axes.
        """
        # Start by creating a copy of our coordinates to experiment with.
        points = [(atom.coords[0], atom.coords[1], atom.coords[2]) for atom in self]

        smallest = math.inf
        x_theta = 0
        y_theta = 0
        z_theta = 0

        if not x_start:
            x_start = start
        
        if not y_start:
            y_start = start
        
        if not z_start:
            z_start = start
        
        if not x_stop:
            x_stop = stop
        
        if not y_stop:
            y_stop = stop
        
        if not z_stop:
            z_stop = stop

        for x_angle in numpy.linspace(x_start, x_stop, steps):
            for y_angle in numpy.linspace(y_start, y_stop, steps):
                for z_angle in numpy.linspace(z_start, z_stop, steps):
                    volume = self.kernel(x_angle, y_angle, z_angle, points)
                    # print("{}, {}, {}, {}".format(x_angle, y_angle, z_angle, volume))
                    if memory is not None:
                        memory.append(({"x_angle": x_angle, "y_angle": y_angle, "z_angle": z_angle}, volume))

                    if volume < smallest:
                        smallest = volume
                        x_theta = x_angle
                        y_theta = y_angle
                        z_theta = z_angle

        # print("Smallest: {}".format(smallest))
        return x_theta, y_theta, z_theta

    
class Nested_grid(Grid_search, Axis_swapper_mix):
    """
    The nested grid search algorithm for determining molecule alignment. 
    """

    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Nested grid", "NEST", "GRID+"]

    def __init__(self, atoms, *grids, charge = None, **kwargs):
        """
        Constructor for Nested_grid.

        :param atoms: List of atom objects to orientate.
        :param grids: The number of grids to optimise over (and the number of steps in each grid).
        :param charge: Optional overall electronic charge of the molecule.
        """
        if len(grids) == 0:
            grids = (35, 10, 10, 10, 10)
        self.grids = grids
        super().__init__(atoms, charge=charge, steps = 0, **kwargs)

    @property
    def human_method_type(self):
        return "{} ({})".format(self.CLASS_HANDLE[0], ":".join((str(ele) for ele in self.grids)))
    
    @property
    def method_type(self):
        return "{}{}".format(self.CLASS_HANDLE[1], ":".join((str(ele) for ele in self.grids)))

    def align_axes(self):
        """
        Realign the axes of our coordinate system.
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # Begin by translating to centre of coordinates.
        coords = self.get_coordinate_list()
        self.translate((-mean(coords[0]) , -mean(coords[1]), -mean(coords[2])))

        # Perform the nested grid search.
        x_start = 0
        x_stop = 0.5*math.pi
        y_start = 0
        y_stop = 0.5*math.pi
        z_start = 0
        z_stop = 0.5*math.pi

        for steps in self.grids:
            x_theta, y_theta, z_theta = self.grid_search(
                steps,
                x_start = x_start,
                y_start = y_start,
                z_start = z_start,
                x_stop = x_stop,
                y_stop = y_stop,
                z_stop = z_stop
            )

            # Decide on parameters for next search.
            x_step = (x_stop - x_start) / (steps - 1)
            x_start = x_theta - x_step
            x_stop = x_theta + x_step

            y_step = (y_stop - y_start) / (steps - 1)
            y_start = y_theta - y_step
            y_stop = y_theta + y_step

            z_step = (z_stop - z_start) / (steps - 1)
            z_start = z_theta - z_step
            z_stop = z_theta + z_step

        # Rotate by the requested amount.
        self.rotate_YZ(x_theta)
        self.rotate_XZ(y_theta)
        self.rotate_XY(z_theta)

        # Swap the axes 
        self.reassign_axes()


class Bayesian_grid(Grid_search, Axis_swapper_mix):
    """
    The bayesian grid search algorithm for determining molecule alignment. 
    """

    def __init__(self, atoms, initial_steps = 1, opt_steps = 50, *args, charge = None, **kwargs):
        """
        Constructor for Bayesian_grid.

        :param atoms: List of atom objects to orientate.
        :param charge: Optional overall electronic charge of the molecule.
        :param initial_steps: The number of steps to search over for each axis before starting the optimiser.
        """
        self.opt_steps = opt_steps
        super().__init__(atoms, *args, charge=charge, steps = initial_steps, **kwargs)

    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Bayesian grid", "BAY"]

    def align_axes(self):
        """
        Realign the axes of our coordinate system.
        
        :return: Nothing. The atoms are rearranged in place.
        """
        if BayesianOptimization is None:
            raise Digichem_exception("Cannot perform Bayesian optimisation; the bayes_opt package is not available")

        # Begin by translating to centre of coordinates.
        coords = self.get_coordinate_list()
        self.translate((-mean(coords[0]) , -mean(coords[1]), -mean(coords[2])))

        # Then create a copy of our coordinates to experiment with.
        points = [(atom.coords[0], atom.coords[1], atom.coords[2]) for atom in self]

        acq = UpperConfidenceBound(
                            kappa=2.576,
                            exploration_decay=1,
                            exploration_decay_delay=0)

        optimizer = BayesianOptimization(
            f = None,
            acquisition_function=acq,
            pbounds = {
                "x_angle": [0, 0.5*math.pi], 
                "y_angle": [0, 0.5*math.pi],
                "z_angle": [0, 0.5*math.pi]
            },
            verbose = 0,
            random_state = 2
        )

        # First do a coarse grid search.
        if self.steps > 0:
            memory = []
            self.grid_search(self.steps, memory = memory)
            # Update the optimizer with the known values.
            for angles, volume in memory:
                optimizer.register(params = angles, target = -volume)
            del memory
        
        try:
            for i in range(self.opt_steps):
                # Get optimizer to suggest new parameter values to try using the
                # specified acquisition function
                next_point = optimizer.suggest()

                # Evaluate the output of the black_box_function using 
                # the new parameter values.
                target = self.kernel(**next_point, points = points)
                
                # Update the optimizer with the evaluation results. 
                optimizer.register(params = next_point, target = -target)
        except NotUniqueError:
            # We're already done.
            pass

        # Rotate by the requested amount.
        self.rotate_YZ(optimizer.max["params"]['x_angle'])
        self.rotate_XZ(optimizer.max["params"]['y_angle'])
        self.rotate_XY(optimizer.max["params"]['z_angle'])

        # Swap the axes 
        self.reassign_axes()