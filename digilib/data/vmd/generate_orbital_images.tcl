# This script expects 8 arguments, which are:
# cube_file    A .cube file to load molecular data from.
set cube_file [lindex $argv 0]
# tcl_common    Path to our common library to source.
set tcl_common [lindex $argv 1]
# A string describing the rendering style to use; 'digichem' or 'gaussian'.
set rendering_style [lindex $argv 2]
# The isovalue to use.
set isovalue [lindex $argv 3]
# A string of x,y,z translations to perform.
set translations [lindex $argv 4]
# rotations        A string list of rotations to perform.
set rotations [lindex $argv 5]
# x0y0z0        A file name to write one of the output images to (in png format).
set x0y0z0 [lindex $argv 6]
# x90y0z0        A file name to write one of the output images to (in png format).
set x90y0z0 [lindex $argv 7]
# x0y90z0        A file name to write one of the output images to (in png format).
set x0y90z0 [lindex $argv 8]
# x45y45z45        A file name to write one of the output images to (in png format).
set x45y45z45 [lindex $argv 9]

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

# Display our obital (both positive and negative phases).
set pos_orbital 1
molecule addrep $mol_handle
standard_orbital_style $pos_orbital $mol_handle 0 $isovalue
#molecule modstyle $pos_orbital $mol_handle isosurface 0.02 0 0 0 1 1
#molecule modcolor $pos_orbital $mol_handle ColorID 1
#molecule modmaterial $pos_orbital $mol_handle Translucent

set neg_orbital 2 
molecule addrep $mol_handle
standard_orbital_style $neg_orbital $mol_handle 1 -$isovalue
#molecule modstyle $neg_orbital $mol_handle isosurface -0.02 0 0 0 1 1
#molecule modcolor $neg_orbital $mol_handle ColorID 0
#molecule modmaterial $neg_orbital $mol_handle Translucent

# And save our pictures.
render_images $rotations $x0y0z0 $x90y0z0 $x0y90z0 $x45y45z45 

# All done.
exit
