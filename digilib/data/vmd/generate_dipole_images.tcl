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
# dipole_start    {x y z} coordinates of the start of the dipole.
set dipole1_start [lindex $argv 5]
# dipole_end    {x y z} coordinates of the end of the dipole.
set dipole1_end [lindex $argv 6]
# dipole_start    {x y z} coordinates of the start of the dipole.
set dipole2_start [lindex $argv 7]
# dipole_end    {x y z} coordinates of the end of the dipole.
set dipole2_end [lindex $argv 8]
# x0y0z0        A file name to write one of the output images to (in png format).
set x0y0z0 [lindex $argv 9]
# x90y0z0        A file name to write one of the output images to (in png format).
set x90y0z0 [lindex $argv 10]
# x0y90z0        A file name to write one of the output images to (in png format).
set x0y90z0 [lindex $argv 11]
# x45y45z45        A file name to write one of the output images to (in png format).
set x45y45z45 [lindex $argv 12]

# Load our common library.
source $tcl_common

# Set our visual style
use_style $rendering_style

# Load our molecule and keep track of its numerical handle.
set mol_handle [molecule new $cube_file]

# Rotate as we've been told.
#rotate_molecule $mol_handle $translations $rotations

# Make our molecule transparent so we can better see our dipole arrow.
material change Opacity Transparent 0.5
molecule modmaterial 0 $mol_handle Transparent

# Use standard display settings.
standard_molecule_style 0 $mol_handle

set dipole1_start [split $dipole1_start ":"]
set dipole1_end [split $dipole1_end ":"]
set dipole2_start [split $dipole2_start ":"]
set dipole2_end [split $dipole2_end ":"]

# Draw our dipole moment.
# Only add them if they are not zero.
if {[lindex $dipole1_end 0] != 0 || [lindex $dipole1_end 1] != 0 || [lindex $dipole1_end 2] != 0} {
    draw color red
    draw arrow $dipole1_start $dipole1_end
}
if {[lindex $dipole2_end 0] != 0 || [lindex $dipole2_end 1] != 0 || [lindex $dipole2_end 2] != 0} {
    draw color green
    draw arrow $dipole2_start $dipole2_end
}

# And save our pictures.
render_images $rotations $x0y0z0 $x90y0z0 $x0y90z0 $x45y45z45 

# All done.
exit
