# This script expects 9 arguments, which are:
# HOMO_cube_file    A .cube file to load molecular data from.
set HOMO_cube_file [lindex $argv 0]
# LUMO_cube_file    A .cube file to load molecular data from.
set LUMO_cube_file [lindex $argv 1]
# tcl_common    Path to our common library to source.
set tcl_common [lindex $argv 2]
# A string describing the rendering style to use; 'digichem' or 'gaussian'.
set rendering_style [lindex $argv 3]
# The isovalue to use.
set isovalue [lindex $argv 4]
# A string of x,y,z translations to perform.
set translations [lindex $argv 5]
# rotations        A string list of rotations to perform.
set rotations [lindex $argv 6]
# x0y0z0        A file name to write one of the output images to (in png format).
set x0y0z0 [lindex $argv 7]
# x90y0z0        A file name to write one of the output images to (in png format).
set x90y0z0 [lindex $argv 8]
# x0y90z0        A file name to write one of the output images to (in png format).
set x0y90z0 [lindex $argv 9]
# x45y45z45        A file name to write one of the output images to (in png format).
set x45y45z45 [lindex $argv 10]

# Load our common library.
source $tcl_common

# Set our visual style
use_style $rendering_style

# Load our molecule structure and keep track of its numerical handle.
set HOMO_mol_handle [molecule new $HOMO_cube_file]

# Rotate it as we've been told.
#rotate_molecule $HOMO_mol_handle $translations $rotations

# Use standard display settings.
standard_molecule_style 0 $HOMO_mol_handle

# Display our obital (both positive and negative phases).
set HOMO_pos_orbital 1
molecule addrep $HOMO_mol_handle
standard_orbital_style $HOMO_pos_orbital $HOMO_mol_handle 0 $isovalue

set HOMO_neg_orbital 2
molecule addrep $HOMO_mol_handle
standard_orbital_style $HOMO_neg_orbital $HOMO_mol_handle 0 -$isovalue


# Also draw our LUMO.
# Load our molecule structure and keep track of its numerical handle.
set LUMO_mol_handle [molecule new $LUMO_cube_file]

# Rotate this one too
#rotate_molecule $LUMO_mol_handle $translations $rotations

# We don't need to show our skeleton again, so we can just reuse the default representation.
# Display our obital (both positive and negative phases).
set LUMO_pos_orbital 0
standard_orbital_style $LUMO_pos_orbital $LUMO_mol_handle 1 $isovalue

set LUMO_neg_orbital 1
molecule addrep $LUMO_mol_handle
standard_orbital_style $LUMO_neg_orbital $LUMO_mol_handle 1 -$isovalue

# And save our pictures.
render_images $rotations $x0y0z0 $x90y0z0 $x0y90z0 $x45y45z45 

# All done.
exit
