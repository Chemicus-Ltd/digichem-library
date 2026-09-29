from pathlib import Path
import subprocess
import math
from PIL import Image
from uuid import uuid4
import os
from math import fabs
import shutil

from digilib.exception.base import File_maker_exception
from digilib.image.render import Render_maker
import digilib.log
from digilib.datas import get_resource


class VMD_image_maker(Render_maker):
    """
    Class for generating image files from Gaussian outputs using VMD.
    
    Nearly all of the work here is done by other programs, primarily VMD (https://www.ks.uiuc.edu/Research/vmd/). These classes are mostly just wrappers.
    """
    
    # The name of the vmd/tcl script (relative to the vmd script folder) used to render images. Inheriting classes should set this to an appropriate script file.
    vmd_script = ""
    # The filename extension to use for molecule scene files produced by vmd.
    scene_file_extension = ".scene"
    # 'Path' to the vmd executable.
    #vmd_executable = "vmd"
    # 'Path' to the tachyon executable.
    #tachyon_executable = "tachyon"
    # The initial resolution at which an image is rendered. This test image will then be discarded once relevant info has been extracted.
    test_resolution = 300
    
    # Name of the section where we get some specific configs.
    options_name = "orbital"
    
    def __init__(self, *args, cube_file = None, rotations = None, auto_crop = True, rendering_style = "pastel", resolution = 1024, also_make_png = True, isovalue = 0.2, 
                 vmd_executable = "vmd", tachyon_executable = "tachyon", vmd_logging = False, num_cpu = 1,
                 **kwargs):
        """
        Constructor for Image_maker objects.
        
        :param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imges will use it as a hint for their output format.
        :param cube_file: The path to a cube_file to use to render images.
        :param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
        :param auto_crop: If False, images will not have excess white space cropped.
        :param rendering_style: A string describing the rendering style to use, either 'digichem' or 'gaussian'.
        :param resolution: The max width or height of the rendered images in pixels.
        :param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats. If 'output' is a .png file, then it is wise to set this option to False (otherwise two png files will be rendered, which is a waste).
        :param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
        :param num_cpu: Number of CPUs to use for multithreading.
        :param vmd_executable: 'Path' to the vmd executable to use for image rendering. Defaults to relying on the command 'vmd'.
        :param tachyon_executable: 'Path' to the tachyon executable to use for image rendering. Defaults to relying on the command 'tachyon'.
        :param vmd_logging: Whether to print output from vmd.
        """
        super().__init__(
            *args,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = auto_crop,
            resolution = resolution,
            also_make_png = also_make_png,
            isovalue = isovalue,
            num_cpu = num_cpu,
            **kwargs
        )
        
        # Some options that control how we function.
        self.rendering_style = rendering_style
        
        # Save executable paths.
        self.vmd_executable = vmd_executable
        self.tachyon_executable = tachyon_executable
        
        self.vmd_logging = vmd_logging
        
        
        if self.rendering_style is None:
            # Panic time.
            raise ValueError("rendering_style must not be None")
        
        # The colours we will be using for rendering (depends on rendering_stype and is only for reference; the actual definition of each style is in VMD).
        if self.rendering_style == "gaussian":
            self.primary_colour = "green"
            self.secondary_colour = "red"
        elif self.rendering_style == "vesta":
            self.primary_colour = "blue"
            self.secondary_colour = "yellow"
        else:
            self.primary_colour = "blue"
            self.secondary_colour = "red"
        
    @property
    def rotations(self):
        # Digichem rotates the wrong way round for some reason, reverse for our rendering engine.
        # VMD also likes degrees not radians.
        return [(axis, math.degrees(-theta)) for axis, theta in self._rotations]
    
    @classmethod
    def from_options(self, output, *, cube_file = None, rotations = None, options, num_cpu = 1, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = options['render']['auto_crop'],
            rendering_style = options['render']['vmd']['rendering_style'],
            resolution = options['render']['resolution'],
            isovalue = options['render'][self.options_name]['isovalue'],
            use_existing = options['render']['use_existing'],
            dont_modify = not options['render']['enable_rendering'],
            vmd_executable = options['render']['vmd']['executable'],
            tachyon_executable = options['render']['vmd']['tachyon'],
            vmd_logging = options['logging']['render_logging'],
            num_cpu = num_cpu,
            **kwargs
        )
    
    @classmethod
    def safe_name(self, file_name):
        """
        Get a filename free of 'unusual' characters.
        :return: The safe path name.
        """
        safe_chars = "._"
        return "".join([char if char.isalnum() or char in safe_chars else "_" for char in file_name])
    
    @property
    def vmd_script_path(self):
        """
        Get the file path to the VMD script to be used by this class to render images (as a pathlib Path).
        """
        return get_resource('data/vmd/{}'.format(self.vmd_script))
    
    @property
    def tcl_common_path(self):
        """
        Get the file path to the Tcl script that is common to all of our other Tcl scripts (as a pathlib Path).
        """
        return get_resource('data/vmd/common.tcl')
    
    @property
    def creation_message(self):
        """
        A short message that may (depending on log-level) be printed to the user before make_files() is called.
        """
        return "Rendering '{}' to file(s)".format(self.output if self.full_path_names else self.output.name)
        
    def make_files(self):
        """
        Make the image files referenced by this object.
        
        The new image will be written to file.
        """
        # Run VMD, which writes a plain text description of our scene to file.
        self.run_VMD_script()
        
        # Next run tachyon for each of our 4 scene files, which creates tga files.
        for image_name in ['x0y0z0', 'x90y0z0', 'x0y90z0', 'x45y45z45']:
            image_path = self.file_path[image_name]
            try:
                # First we'll render a test image at a lower resolution. We'll then crop it, and use the decrease in final resolution to know how much bigger we need to render in our final image to hit our target resolution.
                # Unless of course auto_crop is False, in which case we use our target resolution immediately.
                resolution = self.test_resolution if self.auto_crop else self.target_resolution
                self.run_tachyon_renderer(
                    image_path.with_name(self.safe_name(image_path.with_suffix(self.scene_file_extension).name)),
                    image_path.with_suffix(".tga"),
                    resolution
                )
                
                if self.auto_crop:
                    # Load the test image and autocrop it.
                    with Image.open(image_path.with_suffix(".tga"), "r") as test_im:
                        small_test_im = self.auto_crop_image(test_im)
                        
                    # Get the cropped size. We're interested in the largest dimension, as this is what we'll output as.
                    cropped_resolution = max(small_test_im.size)
                    
                    # From this we can work out the ratio between our true resolution and the resolution we've been asked for.
                    resolution_ratio = cropped_resolution / self.test_resolution
                    
                    # Now we can re-render asking for a larger resolution so that when we crop again, we should get our target resolution.
                    # Note that resolution ratio is a fraction, so self.target_resolution / resolution_ratio results in an increase.
                    self.run_tachyon_renderer(
                        image_path.with_name(self.safe_name(image_path.with_suffix(self.scene_file_extension).name)),
                        image_path.with_suffix(".tga"),
                        self.target_resolution / resolution_ratio
                    )
                
                # Delete the now unneeded scene file.
                os.remove(image_path.with_name(self.safe_name(image_path.with_suffix(self.scene_file_extension).name)))
            except Exception as e:
                raise File_maker_exception(self, "Error in tachyon rendering") from e
            
            # Convert to a better set of formats.
            # Open the file we just rendered.
            with Image.open(image_path.with_suffix(".tga"), "r") as im:
                
                # If we've been asked to autocrop, do so.
                if self.auto_crop:
                    try:
                        cropped_image = self.auto_crop_image(im)
                    except Exception:
                        raise File_maker_exception(self, "Error in post-rendering auto-crop")
                else:
                    cropped_image = im
                
                # Now save in our main output format.
                cropped_image.save(image_path)
                
                # And also as a higher quality png if we've been asked to.
                if self.also_make_png:
                    cropped_image.save(self.file_path[image_name + "_big"])
                
                # And delete the old .tga (unless we've just rendered to .tga!).
                if image_path.suffix != ".tga":
                    os.remove(image_path.with_suffix(".tga"))
        
    @property
    def prepared_translations(self):
        """
        Our list of translations in a from ready for VMD/Tcl.
        """
        return "{},{},{}".format(self.translations[0], self.translations[1], self.translations[2])
    
    @property
    def prepared_rotations(self):
        """
        Our list of rotations in a form ready for VMD/Tcl.
        """
        # This is a stringified form of our list of tuples that we'll pass to VMD (or rather the TCL interpreter).
        if len(self.rotations) > 0:
            rot_string = ":".join(["{},{}".format(axis, angle) for axis, angle in self.rotations])
        else:
            # An empty string will be misinterpreted on the command line.
            rot_string = "0,0"
        return rot_string
    
    @property
    def VMD_signature(self):
        """
        The signature pass to subprocess.run used to call VMD. Inheriting classes should write their own implementation.
        """
        return ""
    
    @property
    def inputs(self):
        """
        A dictionary of all the input files required by VMD for this image.
        
        inputs is a dictionary where each key is the path to the locally accessible file,
        and each value is the true location.
        """
        working_directory = self.output.parent
        return {
            Path(working_directory, self.vmd_script_path.name): self.vmd_script_path,
            Path(working_directory, self.tcl_common_path.name): self.tcl_common_path,
            Path(working_directory, self.safe_name(Path(str(self.input_file)).name)): Path(str(self.input_file)).absolute()
        }
    
    def run_VMD_script(self):
        """
        Called as part of the make() method, inheriting classes can implement this method if they have a different VMD call signature.
        
        Also see VMD_signature, which may do what you need.
        
        :return: The CompletedProcess object. 
        """
        # Setting these variables greatly improves VMD startup time.
        os.environ["VMDNOOPTIX"] = "1"
        os.environ["VMDNOOSPRAY"] = "1"
        
        # VMD has two run scripts, one in C-shell (which is obviously garbage) and one in bash.
        # Sadly, C-Shell is often the default, and it's doesn't support filenames properly.
        # We'll use the usual workaround of changing the working directory and using relative file names.
        working_directory = self.output.parent
        
        try:
            # Create local links to our input files.
            for input_dst, input_src in self.inputs.items():
                # Don't symlink if the source is already in the output dir.
                if input_dst.absolute() != input_src.absolute():
                    try:
                        os.symlink(input_src, input_dst)
                    
                    except Exception:
                        # Couldn't symlink for some reason, print a warning and copy instead.
                        digilib.log.get_logger().warning("Failed to create symlink for VMD input file, falling back to copy (this may take a while)")
                        shutil.copy(input_src, input_dst)
            
            digilib.log.get_logger().debug(str(self.VMD_signature))
            # Run VMD, which renders our image for us.
            try:
                return subprocess.run(
                    self.VMD_signature,
                    # We normally capture and discard stdout (because VMD is VERY verbose), but if we're at a suitable log level, we'll print it.
                    # Nothing useful appears to be printed to stderr, so we'll treat it the same as stdout.,
                    stdin = subprocess.DEVNULL,
                    stdout = subprocess.DEVNULL if not self.vmd_logging else None,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    cwd = working_directory,
                    # VMD has a tendency to sigsegv when closing with VMDNOOPTIX set to on (even tho everything is fine) so we can't check retval sadly.
                    #check = True
                )
            except FileNotFoundError:
                raise File_maker_exception(self, "Could not locate vmd executable '{}'".format(self.vmd_executable))
        
        finally:
            # Clean up.
            # Remove the copied inputs
            for input_dst, input_src in self.inputs.items():
                # Don't symlink if the source is already in the output dir.
                if input_dst.absolute() != input_src.absolute():
                    try:
                        input_dst.unlink()
                    
                    except Exception:
                        # Warnings are useful here, if we can't delete the files we probably failed to copy them in the first place.
                        digilib.log.get_logger().warning("Failed to delete VMD input file '{}'".format(input_dst), exc_info = True)
            
        
    def run_tachyon_renderer(self, scene_file, tga_file, resolution):
        """
        Called as part of the make() method.
        
        Take a molecule 'scene' file produced by VMD and render to an image with tachyon.
        
        :param scene_file: A VMD scene file to read in as input.
        :param tga_file: The filename to write to.
        :param resolution: The resolution to render with.
        """
        # Sadly, tachyon is old and decrepit and doesn't parse the output filename correctly (but weirdly the input file is fine).
        # To get around this, we'll use a semi-random output name that we know is safe, and change-directory to the output directory immediately before calling tachyon.
        # This way we can ensure there will be no nasty file-names for tachyon to get stuck on.
        working_directory = tga_file.parent
        tmpfile_name = ".tachyon_output_" + uuid4().hex + ".tga"
        tmpfile_full_path = Path(tga_file.parent, tmpfile_name)
        tachyon_process = None
        try:
            try:
                # Now we can run tachyon.
                tachyon_process = subprocess.run(
                    [
                        "{}".format(self.tachyon_executable),
                        scene_file.relative_to(working_directory),
                        "-aasamples", "12",
                        # Note: this can get capped in a SLURM context...
                        "-numthreads", "{}".format(self.num_cpu),
                        "-res", "{}".format(resolution), "{}".format(resolution),
                        "-o", tmpfile_name
                    ],
                    stdout = subprocess.PIPE if not self.vmd_logging else None,
                    stderr = subprocess.STDOUT,
                    universal_newlines = True,
                    check = True,
                    cwd = working_directory
                )
            except FileNotFoundError:
                # THINK this can only occur if the command itself cannot be found, but I might be wrong.
                raise File_maker_exception(self, "Could not locate tachyon executable '{}'".format(self.tachyon_executable))
            
            # Render complete (sadly we don't appear to have any ability to check for errors).
            # Now rename our tmpfile.
            #os.rename(tmpfile_full_path, tga_file)
            tmpfile_full_path.rename(tga_file)
        
        except Exception:
            # Something went wrong, remove our tmpfile.
            try:
                os.remove(tmpfile_full_path)
            except Exception:
                pass
            
            # If we didn't already show output, dump it now.
            if not self.vmd_logging and tachyon_process is not None:
                digilib.log.get_logger().error("Tachyon did not exit successfully, dumping output:\n{}".format(tachyon_process.stdout))
                
            raise
        

class Structure_image_maker(VMD_image_maker):
    """
    Class for creating structure images.
    """
    
    vmd_script ="generate_structure_images.tcl"
    
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
                "{}".format(self.vmd_executable),
                "-dispdev", "none",
                "-e", "{}".format(self.vmd_script_path.name),
                "-args",
                "{}".format(self.safe_name(Path(str(self.input_file)).name)),
                "{}".format(self.tcl_common_path.name),
                "{}".format(self.rendering_style),
                "{}".format(self.prepared_translations),
                "{}".format(self.prepared_rotations),
                "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
            ]
        
        
class Orbital_image_maker(Structure_image_maker):
    """
    Class for creating orbital images.
    """
    
    vmd_script ="generate_orbital_images.tcl"
    
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
                "{}".format(self.vmd_executable),
                "-dispdev", "none",
                "-e", "{}".format(self.vmd_script_path.name),
                "-args",
                "{}".format(self.safe_name(Path(str(self.input_file)).name)),
                "{}".format(self.tcl_common_path.name),
                "{}".format(self.rendering_style),
                "{}".format(self.isovalue),
                "{}".format(self.prepared_translations),
                "{}".format(self.prepared_rotations),
                "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
            ]
    

class Difference_density_image_maker(Orbital_image_maker):
    
    # Name of the section where we get some specific configs.
    options_name = "difference_density"


class Spin_density_image_maker(Orbital_image_maker):
    """
    Class for creating spin density images.
    """
        
    vmd_script ="generate_spin_images.tcl"
    
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
                
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
                "{}".format(self.vmd_executable),
                "-dispdev", "none",
                "-e", "{}".format(self.vmd_script_path.name),
                "-args",
                "{}".format(self.safe_name(Path(str(self.input_file)).name)),
                "{}".format(self.tcl_common_path.name),
                "{}".format(self.rendering_style),
                "{}".format(fabs(self.isovalue)),
                "{}".format(self.spin),
                "{}".format(self.prepared_translations),
                "{}".format(self.prepared_rotations),
                "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
            ]
    
class Alpha_orbital_image_maker(Orbital_image_maker):
    pass
    #cubegen_type = "AMO"
    
    
class Beta_orbital_image_maker(Orbital_image_maker):
    pass
    #cubegen_type = "BMO"
    
class Density_image_maker(Orbital_image_maker):
    """
    Class for creating spin density images.
    """
        
    vmd_script ="generate_density_images.tcl"
    
    # Name of the section where we get some specific configs.
    options_name = "density"
    
    def __init__(self, *args, **kwargs):
        """
        Constructor for Spin_density_image_maker objects.
        
        See Orbital_image_maker for a complete constructor.
        """
        super().__init__(*args, **kwargs)
    
    @property
    def type(self):
        """
        The density type.
        """
        return self.input_file.type
                
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
                "{}".format(self.vmd_executable),
                "-dispdev", "none",
                "-e", "{}".format(self.vmd_script_path.name),
                "-args",
                "{}".format(self.safe_name(Path(str(self.input_file)).name)),
                "{}".format(self.tcl_common_path.name),
                "{}".format(self.rendering_style),
                "{}".format(fabs(self.isovalue)),
                "{}".format(self.prepared_translations),
                "{}".format(self.prepared_rotations),
                "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
                "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
            ]
    
class Combined_orbital_image_maker(VMD_image_maker):
    """
    Class for creating images with both the HOMO and LUMO shown together.
    """
    
    vmd_script ="generate_combined_orbital_images.tcl"
    
    def __init__(self, *args, HOMO_cube_file = None, LUMO_cube_file = None, **kwargs):
        """
        Constructor for combined orbital image maker objects.
        
        :param output: Path to write to. See the constructor for VMD_image_maker for how this works.
        :param HOMO_cube_file: Path to the HOMO cube file.
        :param LUMO_cube_file: Path to the LUMO cube file.
        :param *args: See the constructor for VMD_image_maker for further options.
        :param **kwargs: See the constructor for VMD_image_maker for further options.
        """
        super().__init__(*args, cube_file = None, **kwargs)
        self.HOMO_cube_file = HOMO_cube_file
        self.LUMO_cube_file = LUMO_cube_file
        
    @classmethod
    def from_options(self, output, *, HOMO_cube_file = None, LUMO_cube_file = None, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            HOMO_cube_file = HOMO_cube_file,
            LUMO_cube_file = LUMO_cube_file,
            rotations = rotations,
            auto_crop = options['render']['auto_crop'],
            rendering_style = options['render']['vmd']['rendering_style'],
            resolution = options['render']['resolution'],
            isovalue = options['render'][self.options_name]['isovalue'],
            use_existing = options['render']['use_existing'],
            dont_modify = not options['render']['enable_rendering'],
            vmd_executable = options['render']['vmd']['executable'],
            tachyon_executable = options['render']['vmd']['tachyon'],
            vmd_logging = options['logging']['render_logging'],
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
    
    @property
    def inputs(self):
        """
        A dictionary of all the input files required by VMD for this image.
        
        inputs is a dictionary where each key is the path to the locally accessible file,
        and each value is the true location.
        """
        working_directory = self.output.parent
        return {
            Path(working_directory, self.vmd_script_path.name): self.vmd_script_path,
            Path(working_directory, self.tcl_common_path.name): self.tcl_common_path,
            Path(working_directory, self.safe_name(Path(str(self.HOMO_cube_file)).name)): Path(str(self.HOMO_cube_file)).absolute(),
            Path(working_directory, self.safe_name(Path(str(self.LUMO_cube_file)).name)): Path(str(self.LUMO_cube_file)).absolute()
        }
    
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
            "{}".format(self.vmd_executable),
            "-dispdev", "none",
            "-e", "{}".format(self.vmd_script_path.name),
            "-args",
            "{}".format(self.safe_name(Path(str(self.HOMO_cube_file)).name)),
            "{}".format(self.safe_name(Path(str(self.LUMO_cube_file)).name)),
            "{}".format(self.tcl_common_path.name),
            "{}".format(self.rendering_style),
            "{}".format(fabs(self.isovalue)),
            "{}".format(self.prepared_translations),
            "{}".format(self.prepared_rotations),
            "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
        ]


    
class Dipole_image_maker(Structure_image_maker):
    """
    Class for creating dipole images.
    """
    
    vmd_script ="generate_dipole_images.tcl"
    
    
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

    def get_coords(self, dipole, scaling):
        if dipole is None:
            return ( (0.0,0.0,0.0), (0.0,0.0,0.0))
        
        else:
            return (
                tuple([coord * scaling for coord in dipole._origin_coords]),
                tuple([coord * scaling for coord in dipole._vector_coords])
            ) 
        
    @property
    def VMD_signature(self):
        """
        The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
        """
        return [
            "{}".format(self.vmd_executable),
            "-dispdev", "none",
            "-e", "{}".format(self.vmd_script_path.name),
            "-args",
            "{}".format(self.safe_name(Path(str(self.input_file)).name)),
            "{}".format(self.tcl_common_path.name),
            "{}".format(self.rendering_style),
            "{}".format(self.prepared_translations),
            "{}".format(self.prepared_rotations),
            # We don't use the normal origin_coords/vector_coords because these are already rotated, while we want/need to do this rotation with the camera in VMD.
            # Hence use _origin_coords/_vector_coods, which aren't rotated.
            # Dipole 1 (electric).
            "{}:{}:{}".format(*self.get_coords(self.dipole_moment, self.scaling)[0]),
            "{}:{}:{}".format(*self.get_coords(self.dipole_moment, self.scaling)[1]),
            # Dipole 2 (magnetic).
            "{}:{}:{}".format(*self.get_coords(self.magnetic_dipole_moment, self.magnetic_scaling)[0]),
            "{}:{}:{}".format(*self.get_coords(self.magnetic_dipole_moment, self.magnetic_scaling)[1]),
            "{}".format(self.safe_name(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension).name)),
            "{}".format(self.safe_name(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension).name))
        ]

    
class Permanent_dipole_image_maker(Dipole_image_maker):
    """
    Class for creating dipole images.
    """
    
    def __init__(self, *args, dipole_moment = None, scaling = 1, **kwargs):
        """
        Constructor for Dipole_image_maker objects.
        
        :param output: Path to write to. See the constructor for VMD_image_maker for how this works.
        :param cube_file: A Gaussian cube file to use to render the new images.
        :param dipole_moment: A Dipole_moment object that will be rendered as a red arrow in the scene.
        :param scaling: A factor to scale the dipole moment by.
        :param *args: See the constructor for VMD_image_maker for further options.
        :param **kwargs: See the constructor for VMD_image_maker for further options.
        """
        super().__init__(*args, dipole_moment = dipole_moment, scaling = scaling, **kwargs)
        
    @classmethod
    def from_options(self, output, *, dipole_moment = None, cube_file = None, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            cube_file = cube_file,
            rotations = rotations,
            dipole_moment = dipole_moment,
            auto_crop = options['render']['auto_crop'],
            rendering_style = options['render']['vmd']['rendering_style'],
            resolution = options['render']['resolution'],
            isovalue = options['render'][self.options_name]['isovalue'],
            use_existing = options['render']['use_existing'],
            dont_modify = not options['render']['enable_rendering'],
            vmd_executable = options['render']['vmd']['executable'],
            tachyon_executable = options['render']['vmd']['tachyon'],
            vmd_logging = options['logging']['render_logging'],
            scaling = options['render']['dipole_moment']['scaling'],
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
    def from_options(self, output, *, dipole_moment = None, magnetic_dipole_moment, cube_file = None, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            cube_file = cube_file,
            rotations = rotations,
            dipole_moment = dipole_moment,
            magnetic_dipole_moment = magnetic_dipole_moment,
            auto_crop = options['render']['auto_crop'],
            rendering_style = options['render']['vmd']['rendering_style'],
            resolution = options['render']['resolution'],
            isovalue = options['render'][self.options_name]['isovalue'],
            use_existing = options['render']['use_existing'],
            dont_modify = not options['render']['enable_rendering'],
            vmd_executable = options['render']['vmd']['executable'],
            tachyon_executable = options['render']['vmd']['tachyon'],
            vmd_logging = options['logging']['render_logging'],
            scaling = options['render']['transition_dipole_moment']['electric_scaling'],
            magnetic_scaling = options['render']['transition_dipole_moment']['magnetic_scaling'],
            **kwargs
        )
    