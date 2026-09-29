# General imports
import math
import periodictable
from itertools import zip_longest
import re

import digilib.log
from digilib.misc.base import dict_list_index
from digilib.exception.base import Result_unavailable_error, Digichem_exception
from digilib.result import Result_container
from digilib.result import Result_object
from digilib.result import Unmergeable_container_mixin
from digilib.file.prattle import Openprattle_converter

# Hidden import.


def get_chemical_group_mapping(rdkit_molecule):
    """
    Determine chemically equivalent atoms in this atom list.
    
    :return: A mapping between each group number and the atoms it contains.
    """
    from rdkit import Chem

    molecule = rdkit_molecule
    
    groupings = list(Chem.rdmolfiles.CanonicalRankAtoms(molecule, breakTies = False, includeChirality = False))
    
    atoms = list(molecule.GetAtoms())
    
    groups = {}
    for atom_index, atom in enumerate(atoms):
        # For most atoms, just use the assigned group number.
        group_num = groupings[atom_index]
        
        # If this atom is a hydrogen, and it has a single bond to a carbon, use that carbon's group number instead.
        if atom.GetSymbol() == "H":
            bonds = list(atom.GetBonds())
            if len(bonds) == 1 and bonds[0].GetOtherAtom(atom).GetSymbol() == "C":
                # This is an implicit H.
                other_atom = bonds[0].GetOtherAtom(atom)
                group_num = groupings[other_atom.GetIdx()]
        
        try:
            groups[(group_num, atom.GetSymbol())].append(atom_index +1)
        
        except KeyError:
            groups[(group_num, atom.GetSymbol())] = [atom_index +1]
            
    # We now need to re-assign group numbers to ensure they are:
    # 1) Sparse (no missing values, possible because of the hydrogens we removed).
    # 2) Consecutive in terms of the molecule's skeleton (so adjacent numbers are ideally placed
    # adjacent to each other).
    #
    # 2) is challenging to ensure, so for now we simply use the original atom position (because
    # molecules are normally written fairly consecutively).
    # TODO: Find a better way to ensure atom numbering is both deterministic and consecutive.
    
    order = []
    
    for atom_index, atom in enumerate(atoms):
        # Find the current group num of the group this atom is a part of.
        cur_group = dict_list_index(groups, atom_index +1)
        
        # If we've already processed this group, ignore.
        if cur_group[0] not in order:
            # Else, add it next.
            order.append(cur_group[0])
            
    new_groups = {(order.index(group_key[0]) +1, group_key[1]): group for group_key, group in groups.items()}
    
    return new_groups


class Nucleus():
    """
    A class that identifies a specific nucleus (having an element and isotope).
    """
    
    # A regex that can be used to identify a nucleus.
    regex = r'(\d+|\*)?([A-z][A-z]?)'
    
    def __init__(self, element):
        self.isotope = None
        
        if isinstance(element, type(self)):
            self.element = element.element
            self.isotope = element.isotope
        
        elif isinstance(element, list) or isinstance(element, tuple):
            # A list, first item is the element, second (if present) is the isotope).
            self.element = element[0]
            if len(element) > 1:
                self.isotope = element[1]
        
        elif isinstance(element, str):
            self.element, self.isotope = self.split_string(element)
        
        elif isinstance(element, int):
            self.element = periodictable.elements[element]
            
        else:
            raise TypeError("Unknown nucleus type '{}' {}".format(element, type(element)))
        
    @classmethod
    def split_string(self, string):
        """
        """
        match = re.search(self.regex, string)
        
        if match is None:
            raise ValueError("Unable to process nucleus string '{}'".format(string))
        
        match_groups = match.groups()
        isotope = int(match_groups[0]) if match_groups[0] != "*" else None
        
        if isinstance(match_groups[1], int):
            element = periodictable.elements[match_groups[1]]
            
        else:
            element = periodictable.elements.symbol(match_groups[1])
        
        return (element, isotope)
    
    def __str__(self):
        if self.isotope is not None:
            return "{}{}".format(self.isotope, self.element.symbol)
        
        else:
            return self.element.symbol
        

# Molecule_mixin is used for digilib.result.atom.Atom_list and digilib.input.digilib.Digichem_coords_ABC 
class Molecule_mixin():
    """
    Mixin for classes that represent molecules, compounds, or other molecular-like collections of atoms.
    
    Classes that inherit from this mixin must define 'element_dict' (as either an attribute or property).
    See the implementation in Atom_list for what this does.
    """
    
    # TODO: lots of opportunity to move stuff out of Atom_list and into this mixin.
    
    
    @property
    def formula(self):
        """
        Get a formula representation of this atom list.
        
        :return: The formula as a periodictable.formula object (which can be safely cast to string).
        """
        # A dictionary where each key is a type of atom (N, C, H etc) and the value is the number of that atom.
        atoms = self.element_dict
        # Build a string rep.
        form_string = ""
        for atom in atoms:
            form_string += "{}{}".format(atom, atoms[atom])

        # Get and return the formula object.
        return periodictable.formula(form_string)
        
    @property
    def formula_string(self):
        """
        Get a formula representation of this atom list as a string, including optional charge.
        """
        # Get the base formula
        #formula_string = str(self.formula)
        atoms = self.element_dict
        # Build a string rep.
        formula_string = ""
        for atom in atoms:
            formula_string += "{}{}".format(atom, atoms[atom])
        
        # Add charge, if we have one.
        if self.charge == 1:
            formula_string += " +"
        elif self.charge == -1:
            formula_string += " -"
        elif self.charge != 0:
            formula_string += " {}{}".format(abs(self.charge), "-" if self.charge < 0 else "+")
            
        return formula_string
        
    @property
    def molar_mass(self):
        """
        The molar mass of the molecule (takes into account different isotopes and relative isotope abundances, unlike the mass attribute).
        
        :return: The mass (in Daltons / gmol-1).
        """
        return self.formula.mass
    

class Atom_list(Result_container, Unmergeable_container_mixin, Molecule_mixin):
    """
    Class for representing a group of atoms.
    """
    
    # A warning issued when attempting to merge non-equivalent atom lists.
    MERGE_WARNING = "Attempting to merge lists of atoms that are not identical; non-equivalent atoms will be ignored"
    
    def __init__(self, *args, charge = None, assign_groups = True, **kwargs):
        """
        :param charge: The electronic charge of the system.
        :param assign_groups: If True (the default), assign groupings to all the atoms of this list. Call assign_groups() to re-calculate these groupings whenever atoms are added/removed.
        """
        super().__init__(*args, **kwargs)
        self.charge = charge if charge is not None else 0
        self._groups = {}
        
        if assign_groups:
            self.groups
        
    @property
    def mass(self):
        """
        The total mass of all the atoms in this set.
        
        :return: The mass (in Daltons).
        """
        try:
            return sum([atom.mass for atom in self])
        except TypeError:
            # Exact mass not available.
            raise Result_unavailable_error("Exact mass") from None
        
    @property
    def element_dict(self):
        """
        Get a dictionary where each key is one of the elements in the atom list (C, H, N etc) and the value is the number of that element that appears in the atom list.
        
        :return: The element dictionary.
        """
        atoms = {}
        for atom in self:
            # Try and increment the count of the atom.
            try:
                atoms[atom.element.symbol] += 1
            except KeyError:
                # Add the new atom.
                atoms[atom.element.symbol] = 1
        return atoms
    
    @property
    def groups(self):
        """
        Get a list of the groups of atoms in this list.
        
        Atom groups combine chemically equivalent atom positions into a single group.
        """
        if self._groups == {} and len(self) > 0:
            self.assign_groups()
            
        return self._groups
    
    @property
    def bond_matrix(self):
        """
        Get the bond distance matrix for this molecule, which indicates how many bonds separate each atom in the molecule.
        """
        from rdkit import Chem
        return Chem.rdmolops.GetDistanceMatrix(self.to_rdkit_molecule())
    
    def find(self, criteria = None, *, label = None, index = None):
        """
        Find an atom that matches a given criteria
        
        :raises ValueError: If the requested atom could not be found.
        :param criteria: Automatically determine which criteria to search by.
        :param label: The label (a string, such as C1), to find.
        :param index: The index (an int or string that looks like an int) of the atom
        :return: The requested atom.
        """
        if criteria is not None:
            if criteria.isdigit() or isinstance(criteria, int):
                index = int(criteria)
            else:
                label = criteria
        
        if index:
            # Fetch by index.
            try:
                return self[int(index)]
            
            except IndexError:
                raise Result_unavailable_error("Atom", "could not find atom with index {}".format(index)) from None
            
        else:
            # Search by label.
            try:
                return [atom for atom in self if atom.label == label][0]
            
            except IndexError:
                raise Result_unavailable_error("Atom", "could not find atom with label {}".format(label)) from None
    
    def assign_groups(self):
        """
        Assign groupings to all the atoms of this system.
        
        Atom groups combine chemically equivalent atom positions into a single group.
        """
        group_mappings = get_chemical_group_mapping(self.to_rdkit_molecule())            
        self._groups = {group_id: Atom_group(group_id[0], [self[atom_index -1] for atom_index in atom_indices]) for group_id, atom_indices in group_mappings.items()}
        
        for atom in self:
            group = [group for group in self.groups.values() if atom in group.atoms][0]
            
            atom.group = group
    
    @property
    def smiles(self):
        """
        Get this geometry in (canonical) SMILES format.
        """
        try:
            return self._smiles
        
        except AttributeError:
            pass
        
        # Cache miss, go do some work.

        from rdkit.Chem import MolToSmiles
        from rdkit.Chem.rdmolops import RemoveHs
        
        # TODO: Find some other way of generating SMILES.

        mol = self.to_rdkit_molecule()
        try:
            # TODO: rdkit is unreliable, this method can fail for lots of reasons...
            mol = RemoveHs(mol)
        except Exception:
            pass
            
        self._smiles = MolToSmiles(mol)
        return self._smiles

        # # TODO: Handle cases where obabel isn't available
        # conv = Openprattle_converter.get_cls("xyz")(input_file = self.to_xyz(), input_file_type = "xyz")
        # # Cache the result in case we need it again.
        # self._smiles = conv.convert("can").strip()
        # return self._smiles
        
    @property
    def X_length(self):
        return self.get_axis_length(0)
    
    @property
    def Y_length(self):
        return self.get_axis_length(1)
    
    @property
    def Z_length(self):
        return self.get_axis_length(2)
    
    def get_axis_length(self, axis):
        """
        Calculate the length of an axis, defined as the distance required in that axis to contain all the atoms of the set.
        
        :param axis: The axis to calculate for as an integer (0: X-axis, 1: Y-axis, 2: Z-axis).
        :return The length (in angstroms).
        """
        if not 0 <= axis <= 2:
            # Axis is invalid.
            raise ValueError("Axis '{}' is out of bounds. Possible  values are 0 (X), 1 (Y) or 2 (Z)")
        
        # First sort our list of atoms in terms of x, y or z coord.
        sorted_atoms = sorted(self, key = lambda atom: atom.coords[axis])
                
        # Now the axis length is simply the difference between the greatest and the smallest.
        try:
            return sorted_atoms[-1].coords[axis] - sorted_atoms[0].coords[axis]
        
        except IndexError:
            if len(sorted_atoms) == 0:
                # There are not atoms.
                return 0.0
            
            else:
                raise
    
    def get_linear_ratio(self):
        """
        Get the linear ratio of the molecule.
        
        The linear ratio is defined as 1 - (Y_length / X_length).
        
        :return: The ratio, from 0 (non-linear) to 1 (linear).
        """
        try:
            return 1- (self.Y_length / self.X_length)
        
        except (FloatingPointError, ZeroDivisionError):
            return 0
    
    def get_planar_ratio(self):
        """
        Get the planar ratio of the molecule.
        
        The planar ratio is defined as 1 - (Z_length / Y_length).
        
        :return: The ratio, from 0 (non-planar) to 1 (planar).
        """
        try:
            return 1- (self.Z_length / self.Y_length)
        
        except (FloatingPointError, ZeroDivisionError):
            return 0
    
    def get_X_axis_angle(self, start_coord, end_coord):
        """
        Get the angle between a line and the X axis.
        
        :param start_coord: A (X, Y, Z) tuple of coordinates of the start of the line.
        :param end_coord: A (X, Y, Z) tuple of coordinates of the end of the line.
        :return: The angle (in radians). 
        """
        return self.get_theta(math.sqrt( (end_coord[2] - start_coord[2])**2 + (end_coord[1] - start_coord[1])**2 ), end_coord[0] - start_coord[0])
    
    def get_XY_plane_angle(self, start_coord, end_coord):
        """
        Get the angle between a line and the XY plane.
        
        :param start_coord: A (X, Y, Z) tuple of coordinates of the start of the line.
        :param end_coord: A (X, Y, Z) tuple of coordinates of the end of the line.
        :return: The angle (in radians). 
        """
        # The 'secondary' axis is the opposite side of our triangle.
        secondary_axis = end_coord[2] - start_coord[2]
        # The 'primary' axis is the adjacent side of our triangle, which we can get with pythagoras.
        primary_axis = math.sqrt( (end_coord[0] - start_coord[0])**2 + (end_coord[1] - start_coord[1])**2 )
        return self.get_theta(secondary_axis, primary_axis)
    
    @classmethod
    def from_parser(self, parser, *args):
        """
        Get an Atom_list object from an output file parser.
        
        :param parser: An output file parser.
        :param charge: Charge of the system.
        :return: A list of TDM objects.
        """
        return self(Atom.list_from_parser(parser), *args, charge = parser.results.metadata.charge)
    
    @classmethod
    def from_dump(self, data, result_set, options, *args):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(Atom.list_from_dump(data['values'], result_set, options), *args,  charge = data['charge'])
    
    @classmethod
    def from_coords(self, coords):
        """
        Get an instance of this class from a Digichem input coordinates object.
        
        :param coords: Digichem input coords.
        """
        return self(Atom.list_from_coords(coords), charge = coords.charge)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        try:
            smiles = self.smiles
        
        except Exception:
            smiles = None
        
        dump_dict = {
            "formula": self.formula_string,
            "charge": self.charge,
            "smiles": smiles,
            "exact_mass": {
                "value": float(self.mass) if self.safe_get("mass") is not None else None,
                "units": "g mol^-1" 
            },
            "molar_mass": {
                "value": self.molar_mass,
                "units": "g mol^-1",
                },
            "num_atoms": len(self),
            "x-extension": {
                "value": float(self.X_length),
                "units": "Å"
            },
            "y-extension": {
                "value": float(self.Y_length),
                "units": "Å"
            },
            "z-extension": {
                "value": float(self.Z_length),
                "units": "Å"
            },
            "linearity_ratio": float(self.get_linear_ratio()),
            "planarity_ratio": float(self.get_planar_ratio()),
            "values": super()._dump_(digichem_options, all),
        }
        return dump_dict
    
    @classmethod
    def merge(self, *multiple_lists, charge):
        """
        Merge multiple lists of atoms into a single list.
         
        Note that it does not make logical sense to combine different list of atoms into one; hence the method only ensures that all given lists (which are not empty) are the same and then returns the first (non empty) given.
        If the atom lists are not equivalent, a warning will be issued.
        """
        return super().merge(*multiple_lists, charge = charge)
    
    def to_xyz(self):
        """
        Convert this list of atoms to xyz format.
        """
        # First, the number of atoms.
        xyz = "{}\n\n".format(len(self))
        
        # Then coordinates.
        # No effort is made here to truncate coordinates to a certain precision.
        for atom in self:
            xyz += "{}    {:f}    {:f}    {:f}\n".format(atom.element.symbol, float(atom.coords[0]), float(atom.coords[1]), float(atom.coords[2]))
        
        return xyz
    
    def to_mol(self):
        """
        Convert this list of atoms to mol format (useful for reading with rdkit).
        """
        return Openprattle_converter(input_file_buffer = self.to_xyz(), input_file_path = "internal atoms object", input_file_type = "xyz").convert("mol", charge = self.charge)
    
    def to_rdkit_molecule(self):
        """
        Convert this list of atoms to an rdkit molecule object.
        """
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        
        # RDKit has a lot of problems loading molecule information from other formats.
        # Loading from Mol specifies atom and bonding types, but rdkit is very fragile when it comes to
        # which bonds it considers allowed (tetravelent nitrogen is not allowed for example). This can
        # be 'fixed'/ignored by specifying sanitize = False, but then lots of other rdkit functions don't
        # work properly, and rdkit does not detect symmetry correctly.
        # Loading from SMILES is a bad idea because hydrogens are lost.
        # Loading from xyz is much more robust and maintains Hs, but no bond information is read
        # (because there isn't any). RDkit can try and determine bonding itself with DetermineBonds(),
        # but this function is also fragile.
        
        # NOTE: rdkit can fail silently and return None, best to check.
#        mol = Chem.MolFromMolBlock(self.to_mol(), removeHs = False, sanitize = False)
        
        
#         parse_settings = Chem.rdmolfiles.SmilesParserParams()
#         parse_settings.removeHs = False
#         mol = Chem.MolFromSmiles(self.smiles, parse_settings)
#         
        mol = Chem.MolFromXYZBlock(self.to_xyz())
        if mol is None:
            raise Exception("Failed to parse coordinates with rdkit")
        
        mol.UpdatePropertyCache()
        rdDetermineBonds.DetermineConnectivity(mol, charge = self.charge)
        try:
            rdDetermineBonds.DetermineBonds(mol, charge = self.charge)
        
        except Exception:
            #formula_string may also not be implemented...
            try:
                formula_string = self.formula_string

            except Exception:
                formula_string = None
            
            # This function is not implemented for some atoms (eg, Se).
            digilib.log.get_logger().warning(
                "Unable to determine bond ordering for '{}'; all bonds will be represented as single bonds only".format(formula_string)
                , exc_info = True
            )
        
        return mol


class Atom_ABC(Result_object):
    """
    ABC for atom-like result objects.
    """
    
    @property
    def label(self):
        return "{}{}".format(self.element, self.index)
    
    def __hash__(self):
        return hash(self.index)


class Atom_group(Atom_ABC):
    """
    A class that represents a group of atoms that are chemically equivalent.
    """
    
    def __init__(self, index, atoms):
        self.index = index
        self.atoms = atoms
        # Check all elements are the same.
        self.element
        
    @property
    def element(self):
        elements = list(set(atom.element for atom in self.atoms))
        
        if len(elements) > 1:
            raise Digichem_exception("Multiple element types found in atom group '{}'".format(elements))
        
        return elements[0]
    
    @property
    def id(self):
        return (self.index, self.element.symbol)
    
    @property
    def label(self):
        return "{}{}".format(self.element, self.index)
    
    def __str__(self):
        """
        Stringify this atom group.
        """
        return "G:{}".format(self.label)
        
    def __eq__(self, other):
        """
        Is this atom group equal to another?
        """
        return self.atoms == other.atoms
    
    def __hash__(self):
        return super().__hash__()


class Atom(Atom_ABC):
    """
    Class that represents an atom.
    """
    
    def __init__(self, index, atomic_number, coords, mass = None):
        """
        Construct an Atom class.
        
        :param index: The numerical index of this atom, starting from 1.
        :param atomic_number: The atomic/proton number of the atom. Conventional wisdom suggests this has to be an integer, but we make no check here.
        :param mass: The mass of the atom (in Daltons). This isn't always available for some reason (eg, Gaussian Freq calculations), but when it is it identifies the isotope of the given atom.
        :param coords: The coords of the atom, as an (x, y, z) tuple.        
        """
        # Just save each of our attributes.
        self.index = index
        self.mass = mass
        self.coords = coords
        
        # Get our element class.
        try:
            self.element = periodictable.elements[atomic_number]
        except KeyError:
            raise Digichem_exception("No element with atomic number {}".format(atomic_number)) from None
        
        # Atom groups are assigned by Atom_list objects.
        self.group = None
        
    def __str__(self):
        """
        Stringify this atom.
        """
        return "'{}' at '{}'".format(self.label, self.coords)
        
    def __eq__(self, other):
        """
        Is this atom equal to another?
        """
        # Atoms are considered equal if they are the same element in the same position.
        #return self.element == other.element and self.coords == other.coords
        return self.element == other.element and self.index == other.index
    
    def __hash__(self):
        return super().__hash__()
            
    def distance(self, foreign_atom):
        """
        Get the distance between this atom and another atom.
        
        :return: The distance. The units depend on the units of the atoms' coordinates. If the two atoms have coordinates of different units, then you will obviously get bizarre results.
        """
        return math.sqrt( (self.coords[0] - foreign_atom.coords[0])**2 + (self.coords[1] - foreign_atom.coords[1])**2 + (self.coords[2] - foreign_atom.coords[2])**2)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dict = {
            "index": self.index,
            "element": self.element.number,
            "label": self.label,
            "group": self.group.label if self.group else None,
            "coords": {
                "x": {
                    "value": float(self.coords[0]),
                    "units": "Å", 
                },
                "y": {
                    "value": float(self.coords[1]),
                    "units": "Å", 
                },
                "z": {
                    "value": float(self.coords[2]),
                    "units": "Å", 
                }
            },
            "mass": {
                "value": float(self.mass) if self.mass is not None else None,
                "units": "g mol^-1"
            }
        }
        return dump_dict
        
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return [self(atom_dict.get('index', index +1), atom_dict['element'], (float(atom_dict['coords']['x']['value']), float(atom_dict['coords']['y']['value']), float(atom_dict['coords']['z']['value'])), atom_dict['mass']['value']) for index, atom_dict in enumerate(data)]
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of Atom objects from an output file parser.
        
        :param parser: An output file parser.
        :result: A list of Atom objects. An empty list is returned if no atom data is available.
        """
        # First pack our data together to make is easier to loop through.
        try:
            atomnos = parser.data.atomnos
            # Atom coords contains a list for each iteration, we only want the last.
            atomcoords = parser.data.atomcoords[-1]
            atommasses = getattr(parser.data, 'atommasses', [])
            
            # Atommasses sometimes is longer than atomcoords or atomnos.
            # This might be a cclib bug, it is not clear.
            if len(atommasses) > len(atomnos):
                # Take the last values only.
                atommasses = atommasses[-len(atomnos):]
            
        except AttributeError:
            # No atom data available.
            return []
        
        if len(atomnos) != len(atomcoords):
            raise Digichem_exception("Cannot create a list of atoms from {} atomic numbers and {} atomic coordinates".format(len(atomnos), len(atomcoords)))
        
        # Zip.
        zip_data = zip_longest(atomnos, atomcoords, atommasses, fillvalue = None)
        
        # Loop through and rebuild our objects.
        return [self(index+1, atomic_number, (float(coords[0]), float(coords[1]), float(coords[2])), mass) for index, (atomic_number, coords, mass) in enumerate(zip_data)]
    
    @classmethod
    def list_from_coords(self, coords):
        """
        Get a list of Atom objects from a Digichem input coordinates object.
        
        :param coords: Digichem input coords.
        :result: A list of Atom objects. An empty list is returned if no atom data is available.
        """
        return [self(index+1, getattr(periodictable.elements, atom['atom']).number, (atom["x"], atom["y"], atom["z"])) for index, atom in enumerate(coords.atoms)]
