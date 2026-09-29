# General imports.
from PIL import Image
from io import BytesIO

import rdkit.RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

# Digichem imports.
from digilib.image.base import Image_maker


class Skeletal_image_maker(Image_maker):
    """
    A class for rendering skeletal-style molecule structure images.
    """
    
    def __init__(self, output, atoms, *args, abs_resolution = None, rel_resolution = 100, numbering = "group", numbering_font_size = 0.7, explicit_h = False, **kwargs):
        """
        Constructor for Structure_image_maker objects.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg: png, jpeg).
        :param atoms: A list of atom of the molecule/system to render.
        :param resolution: The width and height of the rendered image.
        :param render_backend: The library to prefer to use for rendering, possible options are 'rdkit' or 'obabel'. If 'rdkit' is chosen but is not available, obabel will be used as a fallback. 
        :param numbering: Whether to show atom numberings, either False (off), "atomic" (for atom wise) or "group" (for group wise).
        :param explicit_h: Whether to show explicit H.
        """
        self.atoms = atoms
        self.abs_resolution = abs_resolution
        self.rel_resolution = rel_resolution
        self.numbering = numbering
        self.numbering_font_size = numbering_font_size
        self.explicit_h = explicit_h
        super().__init__(output, *args, **kwargs)
        
    @classmethod
    def from_options(self, output, *, atoms, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            atoms = atoms,
            abs_resolution = options['skeletal_image']['resolution']['absolute'],
            rel_resolution = options['skeletal_image']['resolution']['relative'],
            #render_backend = options['skeletal_image']['render_backend'],
            **kwargs
        )
        
    def make_files(self):
        """
        Make the image described by this object.
        """
        if self.abs_resolution:
            resolution = self.abs_resolution
        
        else:
            resolution = int(self.rel_resolution * self.atoms.X_length)
        
        molecule = self.atoms.to_rdkit_molecule()
        
        # Calculate atom groupings.
        atom_groups = self.atoms.groups
        
        for atom in molecule.GetAtoms():
            
            group = [group for group in atom_groups.values() if atom.GetIdx() +1 in [atom.index for atom in group.atoms]][0]
        
            if self.numbering == "both":
                # Add atom labelling.
                #atom.SetProp("molAtomMapNumber", str(atom.GetIdx()+1))
                atom.SetProp("atomNote",  "{} ({})".format(group.id[0], atom.GetIdx()+1))
            
            elif self.numbering == "atomic":
                atom.SetProp("atomNote",  "{}".format(atom.GetIdx()+1))
            
            elif self.numbering == "group":
                atom.SetProp("atomNote",  "{}".format(group.id[0]))
        
        # Remove C-H, if we've been asked to.
        if not self.explicit_h:
            edit_mol = rdkit.Chem.EditableMol(molecule)
            atoms = list(edit_mol.GetMol().GetAtoms())
            for atom_index, atom in enumerate(reversed(atoms)):
                # If this atom is a hydrogen, and it has a single bond to a carbon, delete it.
                if atom.GetSymbol() == "H":
                    bonds = list(atom.GetBonds())
                    if len(bonds) == 1 and bonds[0].GetOtherAtom(atom).GetSymbol() == "C":
                        # This is an implicit H.
                        edit_mol.RemoveAtom(atom.GetIdx())
                        
            molecule = edit_mol.GetMol()
        
        # Different libraries for generating 2D depictions.
        # Coordgen is supposedly superior.
        #rdkit.Chem.AllChem.Compute2DCoords(molecule)
        
        ps = rdkit.Chem.rdCoordGen.CoordGenParams()
        ps.minimizerPrecision = ps.sketcherBestPrecision
        rdkit.Chem.rdCoordGen.AddCoords(molecule, ps)
        rdkit.Chem.rdDepictor.NormalizeDepiction(molecule)
        
        # Then write the file, setting any options we need to.
        d = rdMolDraw2D.MolDraw2DCairo(resolution, resolution)
        d.drawOptions().annotationFontScale = self.numbering_font_size
        d.DrawMolecule(molecule)
        d.FinishDrawing()
        # To save directly, but we're going to open in PIL to crop.
        #d.WriteDrawingText(str(self.output))
        
        # Crop it to remove whitespace.
        with Image.open(BytesIO(d.GetDrawingText()), "r") as im:
            cropped_image = self.auto_crop_image(im)
        
        cropped_image.save(self.output)
