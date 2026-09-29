from PIL import Image
import numpy as np

from digilib.file import File_maker
from digilib.exception import File_maker_exception

class Cropable_mixin():
    """Mixin class for those that can crop their images."""
    
    @classmethod
    def auto_crop_image(self, im):
        """
        Automatically remove excess white pixels around the outside of our image.
        
        :param im: A PIL Image object.
        :returns: A new PIL Image object with white-space removed.
        """
        # (try to) load our image from file.
        #im = Image.open(file_path, "r")
        im.load()
        
        # This solution was inspired by https://stackoverflow.com/questions/14211340/automatically-cropping-an-image-with-python-pil
        image_data = np.asarray(im)
        if im.mode == "RGBA":
            image_data_bw = image_data[:,:,3]

            # Now we just check which rows and columns are not pure white.
            non_empty_columns = np.where(image_data_bw.max(axis = 0) != 0)[0]
            non_empty_rows = np.where(image_data_bw.max(axis = 1) != 0)[0]
        
        else:
            # Axis 2 is the tuple/array of RGB values for each pixel. We want to know which ones are fully white (255,255,255), so we'll take the min of the 3 values and see if this is equal to 255. If it is, then all the pixels are 255        
            image_data_bw = image_data.min(axis = 2)
        
            # Now we just check which rows and columns are not pure white.
            non_empty_columns = np.where(image_data_bw.min(axis = 0) != 255)[0]
            non_empty_rows = np.where(image_data_bw.min(axis = 1) != 255)[0]

        if len(non_empty_rows) == 0 or len(non_empty_columns) == 0:
            raise ValueError("The image is empty!")

        # Get the bounding box of non-white stuff.
        cropBox = (min(non_empty_rows), max(non_empty_rows), min(non_empty_columns), max(non_empty_columns))
        
        # Copy the image data we want.
        image_data_new = image_data[cropBox[0]:cropBox[1]+1, cropBox[2]:cropBox[3]+1 , :]
        
        # And save over our old file.
        new_image = Image.fromarray(image_data_new)
        return new_image


class Image_maker(File_maker, Cropable_mixin):
    """
    Top level class for image maker objects.
    """
    
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "image"
    
    def __init__(self, *args, enable_rendering = True, **kwargs):
        """
        General constructor for Image_maker objects.
        
        These object are used to represent and make images. Note that a single Image_maker object can represent several image_files
        
        :param output: A path to an output file to write to. How exactly this operates depends on the inheriting class, it is often only used as a base file name. See the class you are using.
        :param dont_modify: Flag that modifies how image creation works. If True, no new images will be written to file.
        :param use_existing: Flag that modifies how image creation works. If True, existing files will be preferentially used if available (set to False to force overwriting existing files).
        """
        super().__init__(*args, existing_file = None, dont_modify = not enable_rendering, **kwargs)
        
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
    def creation_message(self):
        """
        A short message that may (depending on log-level) be printed to the user before make_files() is called.
        """
        return "Rendering '{}' to file(s)".format(self.output if self.full_path_names else self.output.name)
    
        
    def get_constrained_size(self, max_width, max_height):
        """
        Get the maximum possible dimensions of this image that retain the same aspect ratio that don't exceed a set of dimensions.
        
        :param max_width: The maximum width to resize to; set to math.inf for no max.
        :param max_heigh: The maximum height to resize to; set to math.inf for no max.
        """
        # Weasyprint hack for broken max-width.
                
        # Now we open our diagram with pillow.
        im = Image.open(self.get_image())
        
        width, height = im.size
        
        scale_factor = min(max_width/width, max_height/height)
        
        return (width *scale_factor, height * scale_factor)
    
        