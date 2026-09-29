import PIL.Image

from .base import Image_maker

# Remove PIL's built in max image size protection.
PIL.Image.MAX_IMAGE_PIXELS = None