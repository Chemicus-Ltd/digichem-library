# Render something with Beautiful Atoms.
# To use this script, run something like:
#
# # To disabled possibly conflicting packages outside of the conda env.
# PYTHONNOUSERSITE=1
# blender -b -P batoms-renderer.py
#
# Where 'blender' is the path to the Beautiful Atoms Blender executable.


# import debugpy
# debugpy.listen(5678)
# debugpy.wait_for_client()

import addon_utils
def handle_error(exception):
    raise exception
addon_utils.enable("batoms", handle_error = handle_error, default_set=True)

import sys
import argparse
import itertools
import bpy
import yaml
import math
from pathlib import Path
import logging

import ase.io
from batoms import Batoms
from batoms.utils.butils import object_mode
from batoms.render import Render

class Digichem_render(Render):
    def set_viewport_distance_center(self, center=None, padding=0, canvas=None):
        """
        Calculate canvas and direction
        """
        batoms = self.batoms
        if padding is None:
            padding = max(batoms.size) + 0.5
        if center is None:
            center = batoms.get_center_of_geometry()
        self.center = center
        if canvas is None:
            width, height, depth = batoms.get_canvas_box(
                direction=self.viewport, padding=padding
            )
        else:
            width = canvas[0]
            height = canvas[1]
            depth = canvas[2]
        if self.distance < 0:
            self.distance = max(10, depth)

        self.update_camera(width, height, depth / 2)

        # To auto centre the camera, we need to select the molecule as well as all isosurfaces that might be present.
        # Select the molecule.
        self.batoms.obj.select_set(True)

        # Isosurfaces.
        for obj in self.batoms.coll.all_objects:
            if obj.batoms.type == "ISOSURFACE":
                obj.select_set(True)

        # Set camera as active.
        bpy.context.scene.camera = self.camera.obj

        # Manually set a focal point.
        self.camera.lens = 50
        # Auto centre.
        bpy.ops.view3d.camera_to_view_selected()

        self.update_light()

def add_molecule(
        cube_file,
        name,
        rotations = None,
        visible = True,
        isovalues = None,
        isotype = "both",
        primary_color = [0.1, 0.1, 0.9, 0.55],
        secondary_color = [1, 0.058, 0.0, 0.55],
        style = "default"
        ):
    """
    """
    rotations = [] if rotations is None else rotations
    isovalues = [] if isovalues is None else isovalues
    
    surface_settings = []
    for isovalue in isovalues:
        if isotype in ["both", "positive"]:
            surface_settings.append({'level': isovalue, 'color': primary_color, 'material_style': style})
        
        if isotype in ["both", "negative"]:
            surface_settings.append({'level': -isovalue, 'color': secondary_color, 'material_style': style})
    
    
    # Load the input data.
    cube = ase.io.read(cube_file, format="cube", read_data=True, full_output=True)
    
    # The centre of the cube is often offset, fix that by shifting the atoms.
    cube["atoms"].translate(-cube["origin"][0:3])
    
    # Get the mol object.
    mol = Batoms(name, from_ase = cube["atoms"])
    
    # Set some look and feel options.    
    
    # Hide cell boundaries.
    mol.cell.hide = True
    
    # Colour tuning.
    # Carbon to 'black'.
    new_colors = {
        "C": (0.095, 0.095, 0.095, 1),
        "B": (1.0, 0.396, 0.468, 1)
    }

    for atom, color in new_colors.items():
        try:
            mol[atom].color = color
        except AttributeError:
            pass

    # And bonds.
    for key, value in mol.bond.settings.items():
        for atom, color in new_colors.items():
            if key[0] == atom:
                value['color1'] = color

            if key[1] == "C":
                value['color2'] = color
    
    # Change molecule style.
    mol.model_style = 1

    # Slightly increase volume of all atoms.
    for atom in mol.species.keys():
        mol[atom].scale *= 1.25
    
    # Increase volume of H atoms
    mol['H'].scale = 0.75
    
    # Add volumes.
    if len(surface_settings) != 0:
        mol.volumetric_data['surface'] = cube['data']
        
        for index, settings in enumerate(surface_settings):
            mol.isosurface.settings[index+1] = settings
        
        mol.isosurface.draw()
    
    # Now move the entire molecule (isosurface and all) back to the origin.
    mol.obj.select_set(True)
    bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY", center='MEDIAN')
    bpy.context.object.location = [0,0,0]
    
    # If we have any rotations, apply those.
    for axis, angle in rotations:
        # Convert to degree.
        # Not used.
        degree = angle * (180/math.pi)
        
        if axis == 0:
            axis = "x"
        
        elif axis == 1:
            axis = "y"
        
        elif axis == 2:
            axis = "z"
        
        else:
            raise ValueError("Unknown rotation axis '{}'".format(axis))
        
        #mol.rotate(degree, axis)
        
        # Taken from batoms source (collection.py: 51)
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        mol.obj.select_set(True)
        bpy.context.view_layer.objects.active = mol.obj
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(),
                                 center_override = (0,0,0))
    
    if not visible:
        mol.hide = True
    
    return mol

# Adapted from https://blender.stackexchange.com/questions/5898/how-can-i-create-a-cylinder-linking-two-points-with-python?newreg=f372ba9448694f5b97879a6dab963cee
def draw_primitive(start, end, radius, mesh_type, color, collection = None):
    """
    """
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]  
    dist = math.sqrt(dx**2 + dy**2 + dz**2)
    
    # First draw the shape.
    if mesh_type == "cylinder":
        bpy.ops.mesh.primitive_cylinder_add(
            radius = radius, 
            depth = dist,
            vertices = 32,
            location = (dx/2 + start[0], dy/2 + start[1], dz/2 + start[2]) 
        )
    elif mesh_type == "cone":
        bpy.ops.mesh.primitive_cone_add(
            radius1 = radius,
            depth = dist,
            vertices = 32,
            location = (dx/2 + start[0], dy/2 + start[1], dz/2 + start[2])
        )
    else:
        raise ValueError("Unknown mesh type '{}'".format(mesh_type))
        
    # Get a reference to the object we just made.
    obj = bpy.context.active_object
    
    phi = math.atan2(dy, dx)
    theta = math.acos(dz/dist)
    
    bpy.context.object.rotation_euler[1] = theta 
    bpy.context.object.rotation_euler[2] = phi
    
    # Get material
    mat = bpy.data.materials.new(name="Material")
    
    # assign to 1st material slot
    obj.active_material = mat
        
    mat.use_nodes = True
    tree = mat.node_tree
    nodes = tree.nodes
    bsdf = nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = color
    bsdf.inputs["Metallic"].default_value = 0.1
    try:
        # Blener 3.x
        bsdf.inputs["Specular"].default_value = 0.2
    except KeyError:
        # Blender 4.x
        bsdf.inputs["Specular IOR Level"].default_value = 0.2
    
    bsdf.inputs["Roughness"].default_value = 0.2
    mat.diffuse_color = color
    
    # If we've been asked to, asign our new object to a given collection.
    # First unlink from any old collections
    for coll in obj.users_collection:
        # Unlink the object
        coll.objects.unlink(obj)

    # Link each object to the target collection
    collection.objects.link(obj)

    return obj

    
def draw_arrow(start, end, radius, color, split = 0.9, collection = None):
    # Decide what proportion of the total vector to dedicate to the arrow stem and head.
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]
    
    join = (dx * split + start[0], dy * split + start[1], dz * split + start[2])
    cylinder = draw_primitive(start, join, radius, "cylinder", color, collection = collection)
    cone = draw_primitive(join, end, radius*2, "cone", color, collection = collection)

    # Select the two objects and join them together.
    bpy.ops.object.select_all(action='DESELECT')
    cylinder.select_set(True)
    cone.select_set(True)
    bpy.ops.object.join()

    arrow = cone

    # Set the origin of the new combined object to the origin of the arrow.
    bpy.context.scene.cursor.location = start
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')

    return arrow



def main():
    parser = argparse.ArgumentParser(
        prog='Beautiful Atoms Renderer',
        description='Render images with BAtoms')
    
    parser.add_argument("cube_file", help = "Path to the cube file to read")
    parser.add_argument("output", help = "File to write to", nargs="?", default = None)
    parser.add_argument("--second_cube", help = "Optional second cube file to read additional isosurface data from", default = None)
    parser.add_argument("--isovalues", help = "List of isovalues to render", nargs = "*", type = float, default = [])
    parser.add_argument("--isotype", help = "Whether to render positive, negative or both isosurfaces for each isovalue", choices = ["positive", "negative", "both"], default = "both")
    parser.add_argument("--isocolor", help = "The colouring method to use for isosurfaces", choices = ["sign", "cube"], default = "sign")
    parser.add_argument("--primary-color", help = "RGBA for one of the colors to use for isosurfaces", type = float, nargs = 4, default = [0.1, 0.1, 0.9, 0.5])
    parser.add_argument("--secondary-color", help = "RGBA for the other color to use for isosurfaces", type = float, nargs = 4, default = [1, 0.058, 0.0, 0.5])
    parser.add_argument("--style", help = "Material style for isosurfaces", choices = ('default', 'metallic', 'plastic', 'ceramic', 'mirror'), default = "ceramic")
    parser.add_argument("--cpus", help = "Number of parallel CPUs to use for rendering", type = int, default = 1)
    parser.add_argument("--use-gpu", help = "Whether to enable GPU rendering", action = "store_true")
    parser.add_argument("--engine", help = "The render engine to use", choices = ("cycles", "eevee", "workbench"), default = "cycles")
    parser.add_argument("--orientation", help = "The orientation to render from, as x, y, z values", nargs = 3, type = float, default = [0, 0, 0])
    parser.add_argument("--resolution", help = "The output resolution in px", type = int, default = 1024)
    parser.add_argument("--render-samples", help = "The maximum number of render samples, more generally results in higher quality but longer render times", type = int, default = 64)# default = 256)
    parser.add_argument("--rotations", help = "A list of rotations (in JSON) to rotate the molecule to a given alignment. The first item in each list item is the axis to rotate about (0=x, 1=y, 2=z), the second is the angle to rotate by (in radians)", nargs = "*", default = [])
    parser.add_argument("--dipoles", help = "Draw dipoles from a list of the following data (in JSON): 0) start coord, 1) end coord, 2) RGBA color information", nargs = "*", default = [])
    parser.add_argument("--alpha", help = "Override the opacity value for all molecule objects (but not dipoles) to  1- this value, useful for showing dipole arrows more clearly", default = None, type = float)
    parser.add_argument("--perspective", help = "The perspective mode, either orthographic or perspective", default = "perspective", choices = ["perspective", "orthographic"])
    parser.add_argument("--padding", help = "Padding", type = float, default = 1.0)
    parser.add_argument("--multi", help = "Render multiple images, each of a different angle of the scene. Each argument should consist of 6 parts, the x y z position, the resolution, the samples, and the filename (which is appended to 'output')", nargs = 6, default =[], action="append")
    # Both blender and python share the same command line arguments.
    # They are separated by double dash ('--'), everything before is for blender,
    # everything afterwards is for python (except for the first argument, wich is
    # the program name, which is for both).
    if "--" in sys.argv:
        python_argv = sys.argv[sys.argv.index("--") +1:]
    else:
        python_argv = []
    
    args = parser.parse_args(python_argv)
    
    # Batoms or blender will silently set the extension to png if it's not already.
    # This is surprising, so stop now before that happens.
    if args.output is not None and Path(args.output).suffix.lower() != ".png":
        raise ValueError("Output location must have a .png extension")
    
    if args.rotations is not None:
        rotations = [yaml.safe_load(rotation) for rotation in args.rotations]

    if args.multi != []:
        if args.orientation != [0, 0, 0]:
            raise ValueError("You cannot set both --orientation and --multi!")

        if args.resolution != 1024:
            raise ValueError("You cannot set both --resolution and --multi!")
        
        if args.output is not None:
            raise ValueError("You cannot set both 'output' and --multi!")
    
    # Remove the starting cube object.
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

    # Load the input data.
    mol = add_molecule(
        args.cube_file,
        name = "first_mol",
        visible = True,
        rotations = rotations,
        isovalues = args.isovalues,
        isotype = args.isotype,
        primary_color = args.primary_color,
        secondary_color = args.secondary_color if args.isocolor == "sign" else args.primary_color,
        style = args.style
    )
    
    # Uncomment to show atom labels.
    # Needs some tweaking to appear in render (viewport only by default).
    #mol.show_label = 'species'
    mol2 = None
    
    # If we have a second cube, add that too.
    if args.second_cube is not None:
        mol2 = add_molecule(
            args.second_cube,
            name = "second_mol",
            visible = False,
            rotations = rotations,
            isovalues = args.isovalues,
            isotype = args.isotype,
            primary_color = args.primary_color if args.isocolor == "sign" else args.secondary_color,
            secondary_color = args.secondary_color,
            style = args.style
        )
    
    if args.alpha:
        # Set all materials transparent.
        for material in bpy.data.materials:
            try:
                material.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = (1 - args.alpha)
            except Exception as e:
                pass
        
    
    # Draw any dipoles.
    arrows = []
    if args.dipoles is not None:
        
        dipoles = [yaml.safe_load(dipole) for dipole in args.dipoles]
        for start_coord, end_coord, rgba in dipoles:
            try:
                arrows.append(draw_arrow(start_coord, end_coord, 0.1, rgba, collection = mol.coll))
            
            except ZeroDivisionError:
                # Length was zero, move on.
                pass
        
    
    # Setup rendering settings.
    mol.render.engine = args.engine

    # Set up cycles for good quality rendering.
    # Prevents early end to rendering (forces us to use the actual number of samples).
    bpy.context.scene.cycles.use_adaptive_sampling = False
    # Post-processing to remove noise, works well for coloured backgrounds, useless for transparency.
    bpy.context.scene.cycles.use_denoising = True
    # Ray-tracing options
    bpy.context.scene.cycles.max_bounces = 48
    bpy.context.scene.cycles.transparent_max_bounces = 24
    if args.use_gpu:
        bpy.context.scene.cycles.device = "GPU"

    elif args.engine != "cycles":
        logging.warning("Cannot disable GPU rendering with {}".format(args.engine))
    
    # Use maximum compression.
    bpy.context.scene.render.image_settings.compression = 1000
    

    # Change light intensity.
    mol.render.lights["Default"].direction = [0.1, 0.1, 1]
    mol.render.lights["Default"].obj.data.node_tree.nodes["Emission"].inputs[1].default_value = 1
    mol.render.lights["Default"].obj.data.energy = 1
    mol.render.lights["Default"].obj.data.angle = 0.174533

    # Add a second light for depth.
    mol.render.lights.add("Accent1", direction = [1,0.5,0.75])
    mol.render.lights.add("Accent2", direction = [0.5,1,0.75])

    mol.render.lights["Accent1"].obj.data.angle = 0.0872665
    mol.render.lights["Accent1"].obj.data.node_tree.nodes["Emission"].inputs[1].default_value = 1
    mol.render.lights["Accent1"].obj.data.energy = 1.25
    mol.render.lights["Accent2"].obj.data.angle = 0.0872665
    mol.render.lights["Accent2"].obj.data.node_tree.nodes["Emission"].inputs[1].default_value = 1
    mol.render.lights["Accent2"].obj.data.energy = 1.25

    # bpy.data.lights["batoms_light_Default"].node_tree.nodes["Emission"].inputs[1].default_value = 0.45
    # bpy.data.lights["batoms_light_Default"].angle
    #mol.render.lights["Default"].energy=10
    
    # Change view mode.
    if args.perspective == "perspective":
        mol.render.camera.type = "PERSP"
    
    else:
        mol.render.camera.type = "ORTHO"
    
    # Enable to add an outline.
    #bpy.context.scene.render.use_freestyle = True
    
    # We have plenty of memory to play with, use one tile.
    bpy.context.scene.cycles.tile_x = args.resolution
    bpy.context.scene.cycles.tile_y = args.resolution
    bpy.context.scene.cycles.tile_size = args.resolution

    # Performance options.
    bpy.context.scene.render.threads_mode = 'FIXED'
    bpy.context.scene.render.threads = args.cpus
    
    # Set our custom renderer so we can modify zoom etc.
    mol.render = Digichem_render()

    # We have two ways we can change which angle we render from.
    # 1) the viewport keyword arg (which places the camera in a certain location).
    # 2) rotate the molecule.
    #
    # We use option 2, because this gives us more control.

    # Work out how many angles we're rendering from.
    if args.multi == []:
        # Just one.
        targets = [[args.orientation[0], args.orientation[1], args.orientation[2], args.resolution, args.render_samples, args.output]]
    
    else:
        # More than one.
        targets = args.multi

    for x, y, z, resolution, samples, full_file_name in targets:
        # Add args.output and mini_file_name together (useful for --multi).
        orientation = (float(x), float(y), float(z))

        mol.render.resolution = [resolution, resolution]
        # Quality control, more = better and slower.
        bpy.context.scene.cycles.samples = int(samples)
        bpy.context.scene.eevee.taa_render_samples = int(samples)
        mol.obj.delta_rotation_euler = orientation
        
        if mol2 is not None:
            mol2.obj.delta_rotation_euler = orientation

        for arrow in arrows:
            arrow.delta_rotation_euler = orientation

        mol.get_image(viewport = [0,0,1], output = full_file_name, padding = args.padding)
    
    return 0
    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    try:
        sys.exit(main())
    
    except Exception as e:
        logging.error("Error", exc_info = True)
        sys.exit(1)
