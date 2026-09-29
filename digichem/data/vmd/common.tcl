# ----------------------------------------------------------------
# Common Functions.
# ----------------------------------------------------------------

# Reset the rotation of the camera to the aligned X, Y and Z axes.
proc reset_rotation {rotations} {
    # Reset all axes to 0.
    rotate x to 0

    # rotations is a string, we want a list to iterate through.
    set rot_list [split "$rotations" ":"]

    # Iterate through.
    foreach rotation_str $rot_list {
        # Split again to get axes, angle.
        set rotation [split "$rotation_str" ","]
        set axis [lindex $rotation 0]
        set angle [lindex $rotation 1]

        # What we rotate depends on our axes.
        if {$axis == 0} {
            rotate x by $angle
        } elseif {$axis == 1} {
            rotate y by [expr -$angle]
        } elseif {$axis == 2} {
            rotate z by $angle
        } else {
            puts "Unknown axis $axis"
        }
    }
}

# Rotate a molecule by a list of rotations.
proc rotate_molecule {molecule translations rotations} {
    # Split our translations string into a list.
    set trans_list [split $translations ","]

    # Now we get an atom selection encompassing our molecule.
    set sel [atomselect $molecule all]

    # Now move our atoms.
    $sel moveby $trans_list

    # rotations is a string, we want a list to iterate through.
    set rot_list [split $rotations ":"]

    # Iterate through.
    foreach rotation_str $rot_list {
        # Split again to get axes, angle.
        set rotation [split $rotation_str ","]
        set axis [lindex $rotation 0]
        set angle [lindex $rotation 1]

        # What we rotate depends on our axes.
        if {$axis == 0} {
            $sel move [transaxis x $angle]
        } elseif {$axis == 1} {
            $sel move [transaxis y [expr -$angle]]
        } elseif {$axis == 2} {
            $sel move [transaxis z $angle]
        } else {
            puts "Unknown axis $axis"
        }
    }

    # Now we move our molecule back so it's still centered in the scene.
    $sel moveby "[expr -[lindex $trans_list 0]] [expr -[lindex $trans_list 1]] [expr -[lindex $trans_list 2]]"
}

# Taken from the VMD docs.
proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    #set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    set middle [vecadd $start [vecscale 0.85 [vecsub $end $start]]]
    #graphics $mol cylinder $start $middle radius 0.15 filled yes
    graphics $mol cylinder $start $middle radius 0.20 filled yes
    #graphics $mol cone $middle $end radius 0.25
    graphics $mol cone $middle $end radius 0.40 resolution 12
}

# Also taken from the VMD docs http://www.ks.uiuc.edu/Research/vmd/vmd-1.8.4/ug/node181.html#14132
# Find the geometric center of our atoms.
proc geom_center {selection} {
        # Init our total to 0 0 0.
        set total_coords [veczero]

        # Iterate through all our coordinates and add to our total.
        foreach coord [$selection get {x y z}] {
            set total_coords [vecadd $total_coords $coord]
        }

        # And now get our average (we do total_coords * (1 / num_coords) which is the same as total_coords / num_coords).
        return [vecscale [expr 1.0 /[$selection num]] $total_coords]
}

# Move our atoms to be centerd around the geometric center.
proc reset_center {selection} {
    # First get our center.
    set center_coords [geom_center $selection]

    # Now move.
    $selection moveby {-[lindex center_coords 0] -[lindex center_coords 1] -[lindex center_coords 2]}
}

# Render the current scene. Rotations is our alignment rotations string. The four arguments are filenames to render each of the 4 angles to.
proc render_images {rotations x0y0z0 x90y0z0 x0y90z0 x45y45z45} {

    reset_rotation $rotations
    #rotate x to 0
    render Tachyon $x0y0z0

    reset_rotation $rotations
    #rotate x to 0
    rotate x by -90
    render Tachyon $x90y0z0

    reset_rotation $rotations
    #rotate x to 0
    rotate y by 90
    render Tachyon $x0y90z0


    reset_rotation $rotations
    #rotate x to 0
    rotate x by 45
    rotate y by 45
    rotate z by 45
    render Tachyon $x45y45z45
}



# Global flags that control our molecule style.
# The colors of orbitals.
set ORBITAL_PRIMARY_COLOR 0
set ORBITAL_SECONDARY_COLOR 1
set ORBITAL_MATERIAL Translucent

# The radius of atoms in CPK
set SPHERE_SCALE 1.0
# The thickness of bonds in CPK.
set BOND_THICKNESS 0.3

# Switch to a visual style based on a string.
proc use_style {style} {
    if {$style == "pastel"} {
        use_pastel_style
    } elseif {$style == "light-pastel"} {
        use_light_pastel_style
    } elseif {$style == "dark-pastel"} {
        use_dark_pastel_style 
    } elseif {$style == "sharp"} {
        use_sharp_style
    } elseif {$style == "gaussian"} {
        use_gaussian_style
    } elseif {$style == "vesta"} {
        use_vesta_style
    } else {
        error "Unknown rendering style $style"
    }
}

# Set default atom colourings, based on https://sciencenotes.org/molecule-atom-colors-cpk-colors/ (but some have been modified).
proc use_default_atom_colours {} {
    # Set carbon colour to grey.
    color Element C 2
        
    # B to pink
    color Element B 9
    color change rgb 9 1.00 0.709803922 0.709803922
    
    # Mg to green-yellow
    color Element Mg 18
    color change rgb 18 0.5412 1.00 0.00
    
    # Al to metallic-pink
    color Element Al 27
    color change rgb 27 0.749 0.651 0.651
    
    # Si to light-yellow
    color Element Si 17
    color change rgb 17 0.941176471 0.784313725 0.62745098
    
    # P to orange
    color Element P 31
    color change rgb 31 1.00 0.501960784 0.00
    
    # S to yellow
    color Element S 4
    color change rgb 4 1.00 1.00 0.188235294
    
    # Various metals to silver.
    # Color code 14 is the default.
    #color Element Ti 14
    color change rgb 14 0.749019608 0.760784314 0.780392157
    
    # O to red
    color Element O 1
    color change rgb 1 1.00 0.050980392 0.050980392
    
    # N to blue
    color Element N 23
    #color change rgb 23 0.04 0.21 1.00
    color change rgb 23 0.188235294 0.31372549 0.97254902
    
    # F to green.
    color Element F 12
    color change rgb 12 0.564705882 0.878431373 0.31372549
    
    # Cl to green.
    color Element Cl 7
    color change rgb 7 0.121568627 0.941176471 0.121568627
    
    # Br to red.
    color Element Br 29
    color change rgb 29 0.650980392 0.160784314 0.160784314
    
    # I to purple.
    color Element I 11
    color change rgb 11 0.580392157 0.0 0.580392157
}

# Set styles used by all pastel styles.
proc use_common_pastel_style {} {
    # Orbitals are transparent
    global ORBITAL_MATERIAL 
    set ORBITAL_MATERIAL Translucent

    # Set our custom 'Translucent' material.
    material change Ambient Translucent 0.10
    material change Diffuse Translucent 0.64
    material change Specular Translucent 0.00
    material change Shininess Translucent 0.00
    material change Mirror Translucent 0.00
    material change Opacity Translucent 0.62
    material change Outline Translucent 0.00
    material change OutlineWidth Translucent 0.00
    material change TransMode Translucent 1

    # Use default bond thicknes.
    global BOND_THICKNESS SPHERE_SCALE
    set BOND_THICKNESS 0.3
    set SPHERE_SCALE 1.0
}

# Use the 'pastel' visual style.
proc use_pastel_style {} {
    # Set common styles.
    use_default_atom_colours
    use_common_pastel_style
    
    # Orbitals should use red and blue.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 24
    set ORBITAL_SECONDARY_COLOR 30
    
    # Make red and blue lighter.
    #color change rgb 24 0.00 0.72 1.00
    color change rgb 24 0.29 0.29 1.00
    color change rgb 30 1.00 0.30 0.30
}

# Use the light 'pastel' visual style.
proc use_light_pastel_style {} {
    # Set common styles.
    use_default_atom_colours
    use_common_pastel_style
    
    # Orbitals should use red and blue.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 24
    set ORBITAL_SECONDARY_COLOR 30
    
    # Make red and blue lighter.
    #color change rgb 24 0.29 0.29 1.00
    color change rgb 24 0.00 0.72 1.00
    color change rgb 30 1.00 0.30 0.30
}

# Use the dark 'pastel' visual style.
proc use_dark_pastel_style {} {
    # Set common styles.
    use_default_atom_colours
    use_common_pastel_style
    
    # Orbitals should use red and blue.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 24
    set ORBITAL_SECONDARY_COLOR 30
    
    # Make red and blue lighter.
    color change rgb 24 0.05 0.05 1.00
    color change rgb 30 1.00 0.05 0.05
}

# Use the 'sharp' style
proc use_sharp_style {} {
    # Set common styles.
    use_default_atom_colours

    # Orbitals should use red and blue.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 0
    set ORBITAL_SECONDARY_COLOR 1

    # Orbitals are transparent
    global ORBITAL_MATERIAL 
    set ORBITAL_MATERIAL Translucent

    # Set our custom 'Translucent' material.
    material change Ambient Translucent 0.04
    material change Diffuse Translucent 0.70
    material change Specular Translucent 1.00
    material change Shininess Translucent 1.00
    material change Mirror Translucent 0.00
    material change Opacity Translucent 0.30
    material change Outline Translucent 0.00
    material change OutlineWidth Translucent 0.00
    material change TransMode Translucent 1

    # Use default bond thicknes.
    global BOND_THICKNESS SPHERE_SCALE
    set BOND_THICKNESS 0.3
    set SPHERE_SCALE 1.0
}

# Set the visual style to Gaussian-like.
proc use_gaussian_style {} {
    # Set common styles.
    use_default_atom_colours

    # Make grey lighter (for carbon).
    color change rgb 2 0.630000 0.630000 0.630000

    # Orbitals should use red and green.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 19
    set ORBITAL_SECONDARY_COLOR 30

    # Orbitals are opaque
    global ORBITAL_MATERIAL
    set ORBITAL_MATERIAL Opaque

    # Make red and green darker.
    color change rgb 19 0.000000 0.500000 0.000000
    color change rgb 30 0.580000 0.000000 0.000000

    # Use thinner bonds.
    global BOND_THICKNESS SPHERE_SCALE
    set BOND_THICKNESS 0.2
    set SPHERE_SCALE 1.0
}

# Set the visual style to VESTA-like.
proc use_vesta_style {} {
    # Set carbon colour to ochre.
    # We use red2 rather than ochre because ochre is VMD's default for most atoms. This way we only change C's colour rather than all other elements.
    color Element C 29
    color change rgb 29 0.490196 0.286275 0.160784

    # Modify the real ochre so other atoms are more distinct from carbon.
    color change rgb 14 0.35 0.35 0.35

    # N to iceblue
    color Element N 15
    color change rgb 15 0.690196 0.725490 0.901961

    # H to pink
    color Element H 9
    color change rgb 9 1.000000 0.800000 0.800000

    # O to red
    color Element O 30
    color change rgb 30 0.99997 0.01328 0.00000

    # P to lilac
    color Element P 13
    color change rgb 13 0.75557 0.61256 0.76425

    # S to yellow
    color Element S 4

    # The reflectance is different for atoms in vesta and is difficult to approximate, but we give it a go.
    material change Ambient Opaque 0.10
    material change Diffuse Opaque 0.65
    #material change Specular Opaque 0.26
    material change Specular Opaque 0.52
    material change Shininess Opaque 0.51
    material change Mirror Opaque 0.00
    material change Opacity Opaque 1.00
    material change Outline Opaque 0.00
    material change OutlineWidth Opaque 0.00
    material change TransMode Opaque 0

    # Set our orbital colours to yellow and blue (they are opposites of each other?).    
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set ORBITAL_PRIMARY_COLOR 21
    set ORBITAL_SECONDARY_COLOR 4
    color change rgb 21 0.000000 1.000000 1.000000

    # Orbitals are transparent
    global ORBITAL_MATERIAL 
    set ORBITAL_MATERIAL Translucent

    # Set our custom 'Translucent' material.
    material change Ambient Translucent 0.10
    #material change Diffuse Translucent 0.70
    material change Diffuse Translucent 0.64
    material change Specular Translucent 0.00
    material change Shininess Translucent 0.00
    material change Mirror Translucent 0.00
    material change Opacity Translucent 0.62
    material change Outline Translucent 0.00
    material change OutlineWidth Translucent 0.00
    material change TransMode Translucent 1

    # Set the size of our atoms.
    global BOND_THICKNESS SPHERE_SCALE
    set BOND_THICKNESS 0.4
    set SPHERE_SCALE 0.75

}

# Set the standard display options for the body of a molecule.
proc standard_molecule_style {representation molecule} {
    # Change how our molecule looks from sticks (the default) to ball and stick.
    molecule modstyle $representation $molecule CPK

    # Use the more detailed Element colouring method rather than default name.
    mol modcolor $representation $molecule Element

    # Set connecting line thickness
    global BOND_THICKNESS SPHERE_SCALE
    mol modstyle $representation $molecule CPK $SPHERE_SCALE $BOND_THICKNESS 30.000000 30.000000
}

proc standard_orbital_style {representation molecule primary_secondary isovalue} {
    # First decide which colour we're going to use.
    global ORBITAL_PRIMARY_COLOR ORBITAL_SECONDARY_COLOR
    set orbital_color 0
    if {$primary_secondary == 0} {
        set orbital_color $ORBITAL_PRIMARY_COLOR
    } else {
        set orbital_color $ORBITAL_SECONDARY_COLOR
    }

    # Switch the representation to the isosurface mode (which shows orbitals).
    molecule modstyle $representation $molecule isosurface $isovalue 0 0 0 1 1
    
    # Change to our given color.
    molecule modcolor $representation $molecule ColorID $orbital_color

    # Set to target material
    global ORBITAL_MATERIAL
    molecule modmaterial $representation $molecule $ORBITAL_MATERIAL
}

# ----------------------------------------------------------------
# VMD Options.
# ----------------------------------------------------------------
# General options that control the look and feel of our images.

# Set background colour to white.
color Display Background 8

#display projection orthographic
axes location off

