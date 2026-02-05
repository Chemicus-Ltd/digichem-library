import yaml
from deepmerge import Merger
from pathlib import Path

from configurables import Configurable, Options, Option
from configurables.option import Nested_dict_type

from digichem.config.locations import user_config_location
from digichem.misc.io import atomic_write
from digichem.translate import Cube_grid_points


class Auto_type():
    """
    A class that when called, converts a string into a more appropriate type automatically.
    
    The rules for deciding which type to convert to are the same as for parsing yaml using pyyaml, as that is the module that is relied on for conversion.
    """
    
    def __new__(cls, value):
        return yaml.safe_load(value)


class Digichem_options(Configurable):
    """
    Class for holding main digichem options from various sources.
    """
    

    alignment = Option(help = "The default alignment method to use, FIX: fixed, MIN: minimal, FAP: furthest atom pair, AA: average angle, AAA: advanced average angle.", choices = ["FIX", "MIN", "FAP", "AA", "AAA", "GRID", "NEST", "BAY"], default = "NEST")
    angle_units = Option(help = "The default angle units to use, deg: degrees, rad: radians.", choices = ["deg", "rad"], default = "deg")
    
    external = Options(
        help = "Options specifying paths to various external programs that digichem may use. If no path is given, then these programs will simply be executed by name (so relying on OS path resolution to find the necessary executables, which is normally fine.)",
        formchk = Option(help = "This option is deprecated, and does nothing", no_edit = True),
        cubegen = Option(help = "This option is deprecated, and does nothing", no_edit = True),
        cubegen_parallel = Option(help = "What type of parallelism to use with cubegen, multithreaded runs a single instance of cubegen across multiple CPUs, pool runs multiple instances of cubegen", choices = [None, "multithreaded", "pool"], default = "pool")
    )

    parse = Options(
        help = "Options for controlling parsing of config files",
        profiling_rows = Option(help = "The maximum number of rows to parse from the calculation profile file (if available); if more rows than this are available then the data will be downsampled to at most this number of data points", type = int, default = 1000)
    )
    
    skeletal_image = Options(
        help = "Options for controlling the rendering of 2D skeletal images.",
        resolution = Options(help = "The resolution for rendered images", 
            absolute = Option(help = "The output resolution (x and y) in pixels", type = int, default = None),
            relative = Option(help = "The output resolution (x and y) in multiples of the length of the molecule", type = int, default = 100)
        )
    )
    
    logging = Options(
        help = "Options relating to output of error messages. Note that the final logging level is determined by combining both 'log_level' and 'verbose', so a 'log_level' of 'OFF' and 'verbose' of '2' is equal to 'ERROR'.",
        log_level = Option(help = "The level of messages to output, one of OFF (no logging at all), CRITICAL (fewest messages), ERROR, WARNING, INFO or DEBUG (most messages)", choices = ["OFF", "CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], default = "INFO"),
        verbose = Option(help = "Increase the verbosity of the program by this amount. Each integer increase of verbosity will increase 'log_level' by 1 degree.", type = int, default = 0),
        render_logging = Option(help = "Whether to print output from render engines.", type = bool, default = False)
    )
    
    render = Options(
        help = "Options for controlling the appearance of 3D molecule images.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        engine = Option(help = "The rendering engine to use", choices = ["vmd", "batoms"], default = "batoms"),
        vmd = Options(help = "VMD specific options (only applies if engine == 'vmd'",
            executable = Option(help = "Path to the VMD (Visual Molecular Dynamics) executable", default = "vmd"),
            tachyon = Option(help = "The tachyon ray-tracing library, performs the actual rendering. Tachyon is typically packaged with VMD, but often isn't added to the path automatically", default = "tachyon"),
            rendering_style = Option(help =\
"""The render/display mode, changes the appearance of rendered molecules/orbitals.
Possible options are:
 pastel: The default style, uses light, pastel-ish colours for orbitals with low transparency. Normal atom colours.
 light-pastel: Similar to pastel, but with lighter orbital colours.
 dark-pastel: Similar to pastel, but with darker orbital colours, reminiscant of the 'sharp' style.
 sharp: Orbitals are darkly coloured around their edge but have high transparency in their center. Normal atom colours.
 gaussian: A style attempting to mimic that of GaussView.
 vesta: A style attempting to mimic that of VESTA.""",
                choices = ["pastel", "light-pastel", "dark-pastel", "sharp", "gaussian", "vesta"], default = "pastel"
            ),
        ),
        batoms = Options(help = "Beautiful Atoms/Blender specific options (only applies if engine == 'batoms'",
            blender = Option(help = "Path to the blender executable, in which beautiful atoms should be installed", default = "batoms-blender"),
            cpus = Option(help = "The number of CPUs/threads to use. This option is overridden if running in a calculation environment (where it uses the same number of CPUs as the calculation did)", type = int, default = 1),
            render_samples = Option(help = "The number of render samples (or passes) to use. Higher values result in higher image quality and greater render times", type = int, default = 32),
            perspective = Option(help = "The perspective mode", choices = ["orthographic", "perspective"], default = "perspective"),
            stacking = Option(help = "The number of image copies to composite together to avoid transparency artifacts", type = int, default = 10)
            # TODO: Colour options.
        ),
        safe_cubes = Option(help = "Whether to sanitize cubes so older software can parse them (VMD < 1.9.2 etc)", type = bool, default = False),
        use_existing = Option(help =\
"""If True, previously created files will be reused. If False, new images will rendered, replacing the old.
This is on by default for 3D rendered images because they are expensive (time-consuming) to render.""", type = bool, default = True
        ),
        
        auto_crop = Option(help = "Whether to enable  automatic cropping of excess whitespace around the border of generated images.  If False, overall image rendering is likely to take less time, but molecules may only occupy a small portion of the true image.", type = bool, default = True),
        resolution = Option(help = "The target resolution for rendered images. Higher values will increase image quality, at the cost of increased render time and file size.", type = int, default = 512),
        orbital = Options(help = "Specific options for orbital density plots.",
            cube_grid_size = Option(help =\
"""The size of grid used to generate cube files.
Densities are plotted on a 3D grid of points, this option controls how many points there are in the grid (per dimension).
In addition to an integer number of points (~100 is often sufficient), any of the following keywords can also be specified:
Tiny, Small, Medium, Large, Huge or Default.""", type = Cube_grid_points, default = Cube_grid_points("Default")
            ),
            isovalue = Option(help = "The isovalue to use for rendering orbital density.", type = float, default = 0.02)
        ),
        spin = Options(help = "Specific options for spin density plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data. As cubes of this type are rendered with a smaller isovalue, it is often necessary to use a larger grid size than normal to maintain quality.", type = Cube_grid_points, default = Cube_grid_points("Large")),
            isovalue = Option(help = "The isovalue to use for plotting spin density.", type = float, default = 0.0004)
        ),
        density = Options(help = "Specific options for total density plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data.", type = Cube_grid_points, default = Cube_grid_points("Large")),
            isovalue = Option(help = "The isovalue to use for plotting total density.", type = float, default = 0.02)
        ),
        difference_density = Options(help = "Specific options for excited states difference density plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data.", type = Cube_grid_points, default = Cube_grid_points("Default")),
            isovalue = Option(help = "The isovalue to use for difference density plots.", type = float, default = 0.001)
        ),
        natural_transition_orbital = Options(help = "Specific options for NTO plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data.", type = Cube_grid_points, default = Cube_grid_points("Default")),
            isovalue = Option(help = "The isovalue to use for NTO plots.", type = float, default = 0.02)
            
        ),
        dipole_moment = Options(help = "Specific options for permanent dipole moment plots.",
            scaling = Option(help = "The value (x) to scale the TDM by, where 1 D = x Å.", type = float, default = 1.0)
            
        ),
        transition_dipole_moment = Options(help = "Specific options for transition dipole moment plots.",
            electric_scaling = Option(help = "The value (x) to scale the TEDM by, where 1 D = x Å.", type = float, default = 5.0),
            magnetic_scaling = Option(help = "The value (x) to scale the TMDM by, where 1 au = x Å.", type = float, default = 10.0),
        ),
            
    )
    
    orbital_diagram = Options(help = "Options that control the appearance of orbital energy diagrams.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'all': Limits are set so all orbitals are visible.
  - 'auto': Limits are set automatically so the 0 point (ie, 'X' axis), HOMO and LUMO are clearly visible.
  - 'center': Like 'auto' except the lower limit is set more negative so the HOMO-LUMO are closer to the middle of the diagram.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in eV).""", type = Auto_type, default = "auto" 
        ),
        full_axis_lines = Option(help = "If True, black lines are drawn around the boarder of the diagram. If False, a line is drawn for the Y axis but not for the other 3 sides of the diagram.", type = bool, default = False)
    )
    
    excited_states_diagram = Options(help = "Options that control the appearance of excited states energy diagrams.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'all': Limits are set so all excited states are visible.
  - 'auto': Limits are set automatically so the lowest excited state of each multiplicity is clearly visible (S1, D1, T1, Q1... N1 etc).
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in eV).""", type = Auto_type, default = "all"
        ),
        show_dest = Option(help = "Whether or not to show the dE(ST) label.", type = bool, default = True)
    )
    
    absorption_spectrum = Options(help = "Options for controlling the appearance of simulated UV-Vis like absorption spectra.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in oscillator strength).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in nm).""", type = Auto_type, default = "auto"
        ),
        max_width = Option(help =\
"""The maximum image width in pixels.
Absorption graphs will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1200
        ),
        peak_cutoff = Option(help =\
"""The minimum oscillator strength that a peak must have to be shown in the graph, as a fraction ofthe highest peak.
Set to 0 for no cutoff (all peaks shown), which may results in the graph being extended well beyond the drawn peaks (because many peaks are too small to see).
This option has no effect when using manual x limits.""", type = float, default = 0.01
        ),
        y_filter = Option(help = "The minimum y value to simulate using the Gaussian function (y values below this are discarded)", type = float, default = 1e-6),
        x_padding = Option(help = "The amount (in nm) to extend the x axis past the highest/lowest energy peak.", type = int, default = 40),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are eV.", type = float, default = 0.4),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 0.01),
        use_jacobian = Option(help = "Whether or not to use the jacobian transformation to correctly scale the y-axis (see J. Phys. Chem. Lett. 2014, 5, 20, 3497)", type = bool, default = True),
        plot_bars = Option(help = "Whether to plot vertical bars for each excited state.", type = bool, default = True),
        plot_peaks = Option(help = "Whether to plot individual Gaussian functions for each excited state.", type = bool, default = False),
        plot_cumulative_peak = Option(help = "Whether to plot the sum of all Gaussian functions (most closely simulates a real spectrum).", type = bool, default = True)
    )
    
    emission_spectrum = Options(help = "Options for controlling the appearance of simulated emission spectra. 'emission_spectrum' and 'absorption_spectrum 'take the same options.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in oscillator strength).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in nm).""", type = Auto_type, default = "auto"
        ),
        max_width = Option(help =\
"""The maximum image width in pixels.
Absorption graphs will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1200
        ),
        peak_cutoff = Option(help =\
"""The minimum oscillator strength that a peak must have to be shown in the graph, as a fraction of the highest peak.
Set to 0 for no cutoff (all peaks shown), which may results in the graph being extended well beyond the drawn peaks (because many peaks are too small to see).
This option has no effect when using manual x limits.""", type = float, default = 0.01
        ),
        y_filter = Option(help = "The minimum y value to simulate using the Gaussian function (y values below this are discarded)", type = float, default = 1e-6),
        x_padding = Option(help = "The amount (in nm) to extend the x axis past the highest/lowest energy peak.", type = int, default = 40),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are eV.", type = float, default = 0.4),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 0.01),
        use_jacobian = Option(help = "Whether or not to use the jacobian transformation to correctly scale the y-axis (see J. Phys. Chem. Lett. 2014, 5, 20, 3497)", type = bool, default = True),
        plot_bars = Option(help = "Whether to plot vertical bars for each excited state.", type = bool, default = True),
        plot_peaks = Option(help = "Whether to plot individual Gaussian functions for each excited state.", type = bool, default = False),
        plot_cumulative_peak = Option(help = "Whether to plot the sum of all Gaussian functions (most closely simulates a real spectrum).", type = bool, default = True)
    )
    
    IR_spectrum = Options(help = "Options for controlling the appearance of simulated IR like vibrational frequency spectra.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in km/mol).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in cm-1).""", type = Auto_type, default = "auto"
        ),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are cm-1.", type = float, default = 80),
        max_width = Option(help =\
"""The maximum image width in pixels.
IR spectra will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1500),
        y_filter = Option(help = "The minimum y value to simulate using the Gaussian function (y values below this are discarded)", type = float, default = 1e-6),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 1.0)
    )
    
    nmr = Options(help = "Options for controlling simulated NMR spectra",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        coupling_filter = Option(help = "Discard J coupling that is below this threshold (in Hz)", type = float, default = 0.1),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are ppm.", type = float, default = 0.01),
        y_filter = Option(help = "The minimum y value to simulate using the Gaussian function (y values below this are discarded)", type = float, default = 1e-6),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in ppm. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 0.001),
        frequency = Option(help = "The frequency to run the simulated spectrometer at. Larger values will result in narrower coupling. Units are MHz", type = float, default = 100),
        pre_merge = Option(help = "The threshold within which similar peaks will be merged together, before they are split by coupling. This option typically results in faster execution, but more error, compared to post_merge. Units in ppm.", type = float, default = 0.01),
        post_merge = Option(help = "The threshold within which similar peaks will be merged together. Units in ppm.", type = float, default = None),
        standards = Option(help = "The chemical shift of a standard reference peak (in ppm) to use to adjust the spectrum. One for each element.", type = Nested_dict_type, default = Nested_dict_type({
            "H": 31.68766666666667,
            "B": 111.27199999999999,
            "C": 197.90316666666664,
            "F": 183.2779999999999,
            "P": 290.8630000000001,
        })),
        plot_bars = Option(help = "Whether to plot vertical bars for each NMR shift", type = bool, default = False),
        plot_zoom_bars = Option(help = "Whether to plot vertical bars for each NMR shift in zoomed spectra", type = bool, default = True),
        
        # TODO: Validation on these sorts of options is poor and needs looking at.
        isotopes = Option(help = "Isotope specific options. Each key should consist of a tuple of (proton_number, isotope).", type = Nested_dict_type, default = Nested_dict_type({
            # Resonance frequencies calculated at 9.3947 T.
            # 1H, increase fidelity to see more detail.
            "1H": {"frequency": 400, "fwhm": 0.0015, "gaussian_resolution": 0.0001, "coupling_filter": 0.001, "pre_merge": 0.0005},
            # 11B.
            "11B": {"frequency": 128.3},
            # 13C.
            "13C": {"frequency": 100.6},
            # 19F
            "19F": {"frequency": 376.5},
            # 31P
            "31P": {"frequency": 162.0}
        }))
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, allow_unrecognised_options = True, **kwargs)
            
    @classmethod
    def _from_reduce(self, kwargs):
        """
        Function to help un-pickling this class.
        """
        return self(**kwargs)
    
    def __reduce__(self):
        """
        Magic function to help pickling this class.
        """
        return (self._from_reduce, (self.dump(),))
        
    # These methods allow the main digichem options object to be accessed as a dict.
    # This is largely to provide compatibility with legacy code, where the digichem options object was, in fact, a dict.
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
        
    def __delitem__(self, key):
        # Deleting a configurable option doesn't actually delete it (because that wouldn't make much sense),
        # Instead it reverts that option to its default.
        delattr(self, key)
    
    def save(self, path = user_config_location):
        """
        Save the current value of these options to file, so that they will be reloaded on next program start.
        
        :param path: Where to save to (the user's config file by default).
        """
        data = yaml.dump(self.dump())

        path = Path(path)
        
        try:
            path.parent.mkdir(exist_ok = True, parents = True)
            atomic_write(path, data)
        
        except FileNotFoundError as e:
            # We lost the race, give up.
            raise Exception("Failed to write settings to file '{}'; one of the parent directories does not exist".format(path)) from e
    
