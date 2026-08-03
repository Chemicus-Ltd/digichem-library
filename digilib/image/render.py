"""
Classes for rendering 3D images of molecules and densities, primarily with blender and beautiful atoms.

Also see vmd.py for an older render engine.
"""
from pathlib import Path
import os
import copy
import subprocess
import yaml
from PIL import Image
import math
import numpy
import warnings
import json

import digilib.log
from digilib.exception.base import File_maker_exception
from digilib.file.base import File_converter
from digilib.image.base import Cropable_mixin
from digilib.datas import get_resource


class Render_maker(File_converter, Cropable_mixin):
    """
    ABC for classes that make 3D renders from cube files.
    """
        
    # Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
    input_file_type = "cube"
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "render"
    
    def __init__(
            self,
            *args,
            cube_file = None,
            rotations = None,
            auto_crop = True,
            resolution = 1024,
            also_make_png = True,
            isovalue = 0.2,
            num_cpu = 1,
            **kwargs):
        """
        Constructor for Image_maker objects.
        
        :param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imgaes will use it as a hint for their output format).
        :param cube_file: The path to a cube_file to use to render images.
        :param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
        :param auto_crop: If False, images will not have excess white space cropped.
        :param resolution: The max width or height of the rendered images in pixels.
        :param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats.
        :param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
        :param num_cpu: The number of CPUs for multithreading.
        :param blender_executable:
        """
        super().__init__(*args, input_file = cube_file, **kwargs)
        # Save our translations list.
        #self.translations = translations if translations is not None else (0,0,0)
        # This is not used, why not?
        self.translations = (0,0,0)
        # Save our rotations list.
        self._rotations = rotations if rotations is not None else []
        
        # Some options that control how we function.
        self.auto_crop = auto_crop
        self.target_resolution = resolution
        self.also_make_png = also_make_png
        self.isovalue = isovalue
        # TODO: Should we not guess the number of CPUs required?
        self.num_cpu = num_cpu if num_cpu is not None else num_cpu
        
        # TODO: These.
        self.primary_colour = "blue"
        self.secondary_colour = "red"
                
        # These 4 attributes are file paths to the four images we create.
        # We'll keep the same file extension as the was given to us.
        self.file_path = {
            'x0y0z0': self.output.with_suffix(".x0y0z0" + self.output.suffix),
            'x90y0z0': self.output.with_suffix(".x90y0z0" + self.output.suffix),
            'x0y90z0': self.output.with_suffix(".x0y90z0" + self.output.suffix),
            'x45y45z45': self.output.with_suffix(".x45y45z45"  + self.output.suffix),
            # TODO: This should respect also_make_png, but works fine for now (these are actually unused).
            # There are higher quality PNG versions.
             'x0y0z0_big': self.output.with_suffix(".x0y0z0.png"),
             'x90y0z0_big': self.output.with_suffix(".x90y0z0.png"),
             'x0y90z0_big': self.output.with_suffix(".x0y90z0.png"),
             'x45y45z45_big': self.output.with_suffix(".x45y45z45.png")
            }
        
    def get_image(self, name = 'file'):
        """
        Get the path to one of the images that this class represents, rendering the image to file first if necessary.
        
        The functioning of this method is controlled by the dont_modify & use_existing flags.
        
        You can also use the normal python attribute mechanism (either through getattr() or dot notation) to get these paths.
        
        :raises KeyError: If name is not the name of one of the images this class represents.
        :param name: The name of an image to get. Depends on the images in self.file_path.
        :return: A pathlib Path object pointing to the image represented by name, or None if no image could be created.
        """
        return self.safe_get_file(name)
    
    @property
    def rotations(self):
        # Digichem rotates the wrong way round for some reason, reverse for our rendering engine.
        return [(axis, theta) for axis, theta in self._rotations]


class Batoms_renderer(Render_maker):
    """
    Class for rendering images with Beautiful Atoms and Blender.
    """
        
    # Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
    input_file_type = "cube"
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "render"
    
    # TODO: Do we need this? Can probably remove if/when VMD is dropped.
    # Name of the section where we get some specific configs.
    options_name = "orbital"
    
    test_resolution = 300
    test_samples = 5
    
    @property
    def batoms_script_path(self):
        """
        Get the file path to the VMD script to be used by this class to render images (as a pathlib Path).
        """
        return get_resource('data/batoms/batoms-renderer.py')
    
    def __init__(
            self,
            *args,
            cube_file = None,
            rotations = None,
            auto_crop = True,
            resolution = 1024,
            render_samples = 32,
            stack = 3,
            also_make_png = True,
            isovalue = 0.02,
            blender_executable = None,
            cpus = 1,
            num_cpu = 1,
            perspective = "perspective",
            logging = False,
            **kwargs):
        """
        Constructor for Image_maker objects.
        
        :param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imgaes will use it as a hint for their output format).
        :param cube_file: The path to a cube_file to use to render images.
        :param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
        :param auto_crop: If False, images will not have excess white space cropped.
        :param resolution: The max width or height of the rendered images in pixels.
        :param render_samples: The number of render samples, more results in longer render times but higher quality image.
        :param stack: The number of copies of the image to composite together to avoid transparency artifacts.
        :param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats.
        :param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
        :param blender_executable: Bath to the blender executable (can be None to use a default).
        :param cpus: DEPREACTED: Number of parallel threads to render with (use num_cpu instead)
        :param num_cpu: Number of parallel threads to render with.
        :param perspective: Perspective mode (orthographic or perspective)
        """
        if cpus != 1:
            warnings.warn("cpus is deprecated, use num_cpu instead", DeprecationWarning)
        
        if num_cpu == 1 and cpus != 1:
            num_cpu = cpus
        
        super().__init__(
            *args,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = auto_crop,
            resolution = resolution,
            also_make_png = also_make_png,
            isovalue = math.fabs(isovalue),
            num_cpu = num_cpu,
            **kwargs
        )
        
        # Blender specific options.
        self.render_samples = render_samples
        self.perspective = perspective
        self.stack = stack

        self.logging = logging
        
        # Use explicit blender location if given.
        if blender_executable is not None:
            self.blender_executable = blender_executable
        
        else:
            # Otherwise, use a default location.
            self.blender_executable = get_resource('data/batoms/blender/blender')
            
    @classmethod
    def from_options(self, output, *, cube_file = None, rotations = None, cpus = None, num_cpu = 1, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = options['render']['auto_crop'],
            resolution = options['render']['resolution'],
            render_samples = options['render']['batoms']['render_samples'],
            stack = options['render']['batoms']['stacking'],
            isovalue = options['render'][self.options_name]['isovalue'],
            use_existing = options['render']['use_existing'],
            dont_modify = not options['render']['enable_rendering'],
            blender_executable = options['render']['batoms']['blender'],
            # Deprecated...
            cpus = cpus if cpus is not None else options['render']['batoms']['cpus'],
            num_cpu = num_cpu,
            perspective = options['render']['batoms']['perspective'],
            logging = options['logging']['render_logging'],
            **kwargs
        )
        
    def blender_signature(self, *targets, padding = 1.0):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        args = [
            # Blender args.
            f"{self.blender_executable}",
            # Run in background.
            "-b",
            # Our script.
            "-P", f"{self.batoms_script_path}",
            "--",
            # Script specific args.
            f"{self.input_file}",
            # f"{output}",
            # Keywords.
            "--cpus", f"{self.num_cpu}",
            # "--orientation", "{}".format(orientation[0]), "{}".format(orientation[1]), "{}".format(orientation[2]),
            # "--resolution", f"{resolution}",
            # "--render-samples", f"{samples}",
            "--perspective", f"{self.perspective}",
            "--padding", f"{padding}",
        ]
        for orientation, resolution, samples, mini_file_name in targets:
            args.extend([
                "--multi",
                "{}".format(orientation[0]), "{}".format(orientation[1]), "{}".format(orientation[2]),
                f"{resolution}",
                f"{samples}",
                mini_file_name
            ])

        if len(self.rotations) > 0:
            args.append("--rotations")

            # Add rotations.
            for rotation in self.rotations:
                reverse = (1,)
                args.append(json.dumps(rotation if rotation[0] not in reverse else (rotation[0], -rotation[1])))
        
        return args
    
    #def run_blender(self, output, resolution, samples, orientation):

    def run_blender(self, *targets):
        """
        Render a (number of) images with blender.

        :param samples: How many render samples to use.
        :param *targets: Images to render, each is a tuple if (orientation, resolution, samples, path).
        """
        # Render with batoms.
        env = dict(os.environ)
        
        # Disabling the user dir helps prevent conflicting installs of certain packages
        env["PYTHONNOUSERSITE"] = "1"
        #env["PYTHONPATH"] = ":" + env["PYTHONPATH"]
        blender_process = None
        
        # Run Blender, which renders our image for us.
        try:
            blender_process = subprocess.run(
                self.blender_signature(*targets),
                stdin = subprocess.DEVNULL,
                stdout = subprocess.PIPE if not self.logging else None,
                stderr = subprocess.STDOUT,
                universal_newlines = True,
                check = True,
                env = env,
            )
            
            for target in targets:
                if not target[3].exists():
                    raise File_maker_exception(self, "Could not find render file '{}'".format(target[3]))

        except FileNotFoundError:
            raise File_maker_exception(self, "Could not locate blender executable '{}'".format(self.blender_executable))
        
        except subprocess.CalledProcessError as e:
            if not self.logging:
                digilib.log.get_logger().error("Blender did not exit successfully, dumping output:\n{}".format(e.stdout))

        except Exception:
            if not self.logging and blender_process is not None:
                digilib.log.get_logger().error("Blender did not exit successfully, dumping output:\n{}".format(blender_process.stdout))
            
            raise
        
    
    def make_files(self):
        """
        Make the image files referenced by this object.
        
        The new image will be written to file.
        """
        # TODO: This mechanism is clunky and inefficient if only one image is needed because its based off the old VMD renderer. With batoms we can do much better.
        angles = {
            "x0y0z0": [
                (0,0,0),
                self.test_resolution if self.auto_crop else self.target_resolution,
                self.test_samples if self.auto_crop else self.render_samples,
                self.file_path['x0y0z0'].with_suffix(".tmp.png")
            ],
            "x90y0z0": [
                (-1.5708, 0, 0),
                self.test_resolution if self.auto_crop else self.target_resolution,
                self.test_samples if self.auto_crop else self.render_samples,
                self.file_path['x90y0z0'].with_suffix(".tmp.png")
            ],
            "x0y90z0": [
                (0, 1.5708, 0),
                self.test_resolution if self.auto_crop else self.target_resolution,
                self.test_samples if self.auto_crop else self.render_samples,
                self.file_path['x0y90z0'].with_suffix(".tmp.png")
            ],
            "x45y45z45": [
                (0.785398, 0.785398, 0.785398),
                self.test_resolution if self.auto_crop else self.target_resolution,
                self.test_samples if self.auto_crop else self.render_samples,
                self.file_path['x45y45z45'].with_suffix(".tmp.png")
            ],
        }
        
        try:
            # First we'll render a test image at a lower resolution. We'll then crop it, and use the decrease in final resolution to know how much bigger we need to render in our final image to hit our target resolution.
            # Unless of course auto_crop is False, in which case we use our target resolution immediately.
            self.run_blender(*list(angles.values()))
            
            if self.auto_crop:
                # Load the test image and autocrop it.
                for angle, target in angles.items():
                    with Image.open(target[3], "r") as test_im:
                        small_test_im = self.auto_crop_image(test_im)
                        
                    # Get the cropped size. We're interested in the largest dimension, as this is what we'll output as.
                    cropped_resolution = max(small_test_im.size)
                    
                    # From this we can work out the ratio between our true resolution and the resolution we've been asked for.
                    resolution_ratio = cropped_resolution / self.test_resolution

                    # Update the target properties.
                    angles[angle] = [
                        # Orientation is the same.
                        angles[angle][0],
                        # New resolution.
                        int(self.target_resolution / resolution_ratio),
                        # New samples.
                        self.render_samples,
                        # New filename,
                        angles[angle][3]
                    ]
                    
                self.run_blender(*list(angles.values()))
            
        except Exception:
            raise File_maker_exception(self, "Error in blender rendering")
        
        # Convert to a better set of formats.
        # Open the files we just rendered.
        for angle, target in angles.items():
            with Image.open(target[3], "r") as im:
                
                # If we've been asked to autocrop, do so.
                if self.auto_crop:
                    try:
                        cropped_image = self.auto_crop_image(im)
                    except Exception:
                        raise File_maker_exception(self, "Error in post-rendering auto-crop")
                else:
                    cropped_image = im

                # Transparency stacking.
                # This 'hack' takes several copies of the same transparent image and layers them atop each other.
                # This avoids problems with isosurfaces being too transparent with respect to the background (and
                # thus loosing definition).
                stacked = cropped_image.copy()

                for _ in range(self.stack):
                    stacked.alpha_composite(cropped_image)

                cropped_image = stacked
                
                # Save as a higher quality png if we've been asked to.
                if self.also_make_png:
                    cropped_image.save(self.file_path[angle + "_big"])
                    
                # Remove transparency, which isn't supported by JPEG (which is essentially the only format we write here).
                # TODO: Check if the output format can support transparency or not.
                new_image = Image.new("RGBA", cropped_image.size, "WHITE") # Create a white rgba background
                new_image.paste(cropped_image, (0, 0), cropped_image)              # Paste the image on the background. Go to the links given below for details.
                cropped_image = new_image.convert('RGB')
                
                # Now save in our main output format.
                cropped_image.save(self.file_path[angle])
                
                # And delete the old .png.
                os.remove(target[3])


class Structure_image_maker(Batoms_renderer):
    """
    Class for creating structure images.
    """
        
        
class Orbital_image_maker(Structure_image_maker):
    """
    Class for creating orbital images.
    """
    
    def blender_signature(self, output, resolution, samples, orientation, isotype = "both", isocolor = "sign"):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        sig = super().blender_signature(output, resolution, samples, orientation)
        sig.extend(["--isovalues", f"{self.isovalue}"])
        sig.extend(["--isotype", isotype])
        sig.extend(["--isocolor", isocolor])
        
        return sig
    

class Difference_density_image_maker(Orbital_image_maker):
    
    # Name of the section where we get some specific configs.
    options_name = "difference_density"


class Spin_density_image_maker(Orbital_image_maker):
    """
    Class for creating spin density images.
    """
    
    # Name of the section where we get some specific configs.
    options_name = "spin"
    
    def __init__(self, *args, spin = "both", **kwargs):
        """
        Constructor for Spin_density_image_maker objects.
        
        See Orbital_image_maker for a complete constructor.
        :param spin: A string indicating which net-spins to render, either 'positive', 'negative' or 'both'.
        """
        super().__init__(*args, **kwargs)
        self.spin = spin
                
    def blender_signature(self, output, resolution, samples, orientation):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        sig = super().blender_signature(output, resolution, samples, orientation, isotype = self.spin)
        return sig


class Alpha_orbital_image_maker(Orbital_image_maker):
    pass


class Beta_orbital_image_maker(Orbital_image_maker):
    pass


class Density_image_maker(Orbital_image_maker):
    """
    Class for creating spin density images.
    """
    
    # Name of the section where we get some specific configs.
    options_name = "density"
    
    @property
    def type(self):
        """
        The density type.
        """
        return self.input_file.type
    
    def blender_signature(self, output, resolution, samples, orientation):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        # There is no negative density for total density.
        sig = super().blender_signature(output, resolution, samples, orientation, isotype = "positive")
        return sig


class Combined_orbital_image_maker(Orbital_image_maker):
    """
    Class for creating images with both the HOMO and LUMO shown together.
    """
    
    vmd_script ="generate_combined_orbital_images.tcl"
    
    def __init__(self, *args, cube_file = None, LUMO_cube_file = None, **kwargs):
        """
        Constructor for combined orbital image maker objects.
        
        :param output: Path to write to. See the constructor for Batoms_renderer for how this works.
        :param HOMO_cube_file: Path to the HOMO cube file.
        :param LUMO_cube_file: Path to the LUMO cube file.
        :param *args: See the constructor for Batoms_renderer for further options.
        :param **kwargs: See the constructor for VMD_image_maker for further options.
        """
        super().__init__(*args, cube_file = cube_file, **kwargs)
        self.LUMO_cube_file = LUMO_cube_file
        
    @property
    def HOMO_cube_file(self):
        return self.input_file
        
    @classmethod
    def from_options(self, output, *, HOMO_cube_file = None, LUMO_cube_file = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return super().from_options(
            output,
            cube_file = HOMO_cube_file,
            LUMO_cube_file = LUMO_cube_file,
            options = options,
            **kwargs
        ) 

    def check_can_make(self):
        """
        Check whether it is feasible to try and render the image(s) that we represent.
        
        Reasons for rendering not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
        
        This method returns nothing, but will raise a File_maker_exception exception if the rendering is not possible.
        """
        try:
            if self.HOMO_cube_file is None or self.HOMO_cube_file.safe_get_file() is None:
                raise File_maker_exception(self, "No HOMO cube file is available")
        except AttributeError:
            if not hasattr(self.HOMO_cube_file, 'safe_get_file'):
            # Input file does not have a safe_get_file() method.
                pass

            else:
                raise

        try:
            if self.LUMO_cube_file is None  or self.LUMO_cube_file.safe_get_file() is None:
                raise File_maker_exception(self, "No LUMO cube file is available")
        
        except AttributeError:
            if not hasattr(self.LUMO_cube_file, 'safe_get_file'):
            # Input file does not have a safe_get_file() method.
                pass

            else:
                raise
    
    def blender_signature(self, output, resolution, samples, orientation):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        sig = super().blender_signature(output, resolution, samples, orientation, isotype = "both", isocolor = "cube")
        sig.extend(["--second_cube", f"{self.LUMO_cube_file}"])
        return sig


class Dipole_image_maker(Structure_image_maker):
    """
    Class for creating dipole images.
    """
    
    
    def __init__(self, *args, dipole_moment = None, magnetic_dipole_moment = None, scaling = 1, magnetic_scaling = 1, **kwargs):
        """
        Constructor for Dipole_image_maker objects.
        
        :param output: Path to write to. See the constructor for VMD_image_maker for how this works.
        :param cube_file: A Gaussian cube file to use to render the new images.
        :param dipole_moment: A Dipole_moment object that will be rendered as a red arrow in the scene.
        :param magnetic_dipole_moment: A second dipole moment object that will be rendered as a blue arrow in the scene.
        :param scaling: A factor to scale the dipole moment by.
        :param magnetic_scaling: A factor to scale the magnetic dipole moment by.
        :param *args: See the constructor for VMD_image_maker for further options.
        :param **kwargs: See the constructor for VMD_image_maker for further options.
        """
        super().__init__(*args, **kwargs)
        self.dipole_moment = dipole_moment
        self.magnetic_dipole_moment = magnetic_dipole_moment
        self.scaling = scaling
        self.magnetic_scaling = magnetic_scaling
        self.electric_arrow_colour = "red"
        self.magnetic_arrow_colour = "green"
        
    def get_dipoles(self):
        """
        Get required dipole information as a nested list.
        """
        dipoles = []
        # First electric (if we have it).
        for dipole, scaling, colour in [
            (self.dipole_moment, self.scaling, [1.0, 0.0, 0.0, 1.0]),
            (self.magnetic_dipole_moment, self.magnetic_scaling, [0.0, 1.0, 0.0, 1.0]),
        ]:
            if dipole is not None:
                dipoles.append([
                    [float(coord * scaling) for coord in dipole.origin_coords],
                    [float(coord * scaling) for coord in dipole.vector_coords],
                    colour
                ])
        
        return dipoles
    
    def blender_signature(self, output, resolution, samples, orientation):
        """
        The signature passed to subprocess.run used to call Blender. Inheriting classes should write their own implementation.
        """
        sig = super().blender_signature(output, resolution, samples, orientation)
        sig.append("--dipoles")
        for dipole in self.get_dipoles():
            #sig.append(yaml.safe_dump(dipole))
            sig.append(json.dumps(dipole))
        sig.extend(["--alpha", "0.75"])
        return sig

    
class Permanent_dipole_image_maker(Dipole_image_maker):
    """
    Class for creating dipole images.
    """
        
    @classmethod
    def from_options(self, output, *, dipole_moment = None, cube_file = None, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """
        return super().from_options(
            output,
            dipole_moment = dipole_moment, 
            cube_file = cube_file,
            rotations = rotations,
            scaling = options['render']['dipole_moment']['scaling'],
            options = options,
            **kwargs
        )
    
    def check_can_make(self):
        """
        Check whether it is feasible to try and render the image(s) that we represent.
        
        Reasons for rendering not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
        
        This method returns nothing, but will raise a File_maker_exception exception if the rendering is not possible.
        """ 
        super().check_can_make()
        
        # Also make sure we have a dipole.
        if self.dipole_moment is None:
            raise File_maker_exception(self, "No dipole moment is available.")
    

class Transition_dipole_image_maker(Dipole_image_maker):
    """
    Class for creating TDM images.
    """
    
    def check_can_make(self):
        """
        Check whether it is feasible to try and render the image(s) that we represent.
        
        Reasons for rendering not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
        
        This method returns nothing, but will raise a File_maker_exception exception if the rendering is not possible.
        """ 
        super().check_can_make()
        
        # Also make sure we have a dipole.
        if self.dipole_moment is None and self.magnetic_dipole_moment is None:
            raise File_maker_exception(self, "No dipole moment is available.")
        
    @classmethod
    def from_options(self, output, *, dipole_moment = None, cube_file = None, magnetic_dipole_moment, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """
        return super().from_options(
            output,
            dipole_moment = dipole_moment,
            magnetic_dipole_moment = magnetic_dipole_moment,
            cube_file = cube_file,
            rotations = rotations,
            scaling = options['render']['dipole_moment']['scaling'],
            magnetic_scaling = options['render']['transition_dipole_moment']['magnetic_scaling'],
            options = options,
            **kwargs
        )   
