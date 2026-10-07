# This script expects 9 arguments, which are:
# cube_file    A .cube file to load molecular data from.
set cube_file [lindex $argv 0]
# tcl_common    Path to our common library to source.
set tcl_common [lindex $argv 1]
# A string describing the rendering style to use; 'digichem' or 'gaussian'.
set rendering_style [lindex $argv 2]
# A string of x,y,z translations to perform.
set translations [lindex $argv 3]
# rotations        A string list of rotations to perform.
set rotations [lindex $argv 4]
# x0y0z0        A file name to write one of the output images to (in png format).
set x0y0z0 [lindex $argv 5]
# x90y0z0        A file name to write one of the output images to (in png format).
set x90y0z0 [lindex $argv 6]
# x0y90z0        A file name to write one of the output images to (in png format).
set x0y90z0 [lindex $argv 7]
# x45y45z45        A file name to write one of the output images to (in png format).
set x45y45z45 [lindex $argv 8]

# Load our common library.
source $tcl_common

# Set our visual style
use_style $rendering_style

# Load our molecule and keep track of its numerical handle.
set mol_handle [molecule new $cube_file]

# Rotate as we've been told.
#rotate_molecule $mol_handle $translations $rotations

# Use standard display settings.
standard_molecule_style 0 $mol_handle

# And save our pictures.
render_images $rotations $x0y0z0 $x90y0z0 $x0y90z0 $x45y45z45 

# All done.
exit
