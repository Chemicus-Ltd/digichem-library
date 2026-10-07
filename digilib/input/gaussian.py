from digilib.exception import Digichem_exception


class Gaussian_input_parser():
    """
    Class for parsing Gaussian input files.
    """
    
    # Each 'section' in Gaussian is separated by a double newline; this seems to be the only consistent delimiter.
    SECTION_SEPARATOR = "\n\n"
    
    def __init__(self, file_str = None):
        """
        Constructor for Gaussian input files.
        
        :param file_str: String containing a Gaussian input file.
        """
        # Link 0 commands that appear at the top of the input file, we currently do no parsing on these (because we ignore them anyway)
        # Link 0 commands are stored as they appear (as a string, possibly containing newlines).
        self.link_0 = None
        
        # The route section, describes the calculation to be performed.
        ## Most of the route is ignored as we set it ourselves as part of the submission process, but some options (geom=connectivity, genECP etc) control the format of the rest of the file, so we do parse these.
        ## This is a dictionary of options; where 'key: value' translates to 'key=value' in the gaussian file. If value is None, then the translation is to 'key'.
        # The route section is currently not parsed (it is just a string).        
        self.route = None
        
        # The title of the input file.
        self.title = None
        
        # The multiplicity and charge of the molecule, appears at the top of the geometry section as 'charge, mult' or 'charge mult'.
        self.multiplicity = None
        self.charge = None
        
        # The geometry (atoms and charge etc) section as a string (almost certainly containing newlines).
        self.geometry = None
        
        # These are additional sections that can optionally appear in some input files (connectivity, basis set etc).
        self.additional_sections = []
        # If we've been given a file, load it now.
        if file_str is not None:
            self.load(file_str)
        
    @property
    def title(self):
        """
        Get the title section of this input file.
        
        None and empty string values are translated to a single whitespace character (because otherwise they will be interpreted wrong).
        """
        if self._title is not None and len(self._title.trim()) != 0:
            return self._title
        
        else:
            return "(title)"
    
    @title.setter
    def title(self, value):
        """
        Set the title section of this input file.
        """
        self._title = value
    
    def load(self, file_str):
        """
        Load a Gaussian input file.
        
        :param file_str: String containing a Gaussian input file.
        
        """
        # First, split on our delimeter.
        sections = file_str.split("\n\n")
        
        link_0_lines = []
        route_lines = []
        # Split the first section into link 0 and route (link 0 starts with %, the first line not to start without % is route).
        if len(sections) > 0:
            in_route = False
            for line in sections[0].split("\n"):
                # Check first char.
                if not in_route and line[:1] != "%":
                    in_route = True
                    
                # Add to one of our two lists.
                if not in_route:
                    link_0_lines.append(line)
                else:
                    route_lines.append(line)
        
        # Now set.
        self.link_0 = "\n".join(link_0_lines) if len(link_0_lines) > 0 else None
        self.route = "\n".join(route_lines) if len(route_lines) > 0 else None
        
        # Set title.
        self.title = sections[1] if len(sections) > 1 else None
        
        # Geometry (the first line contains charge and mult).
        try:
            geometry_section = sections[2].split("\n", 1)
        except IndexError:
            # No geometry available.
            raise Digichem_exception("Failed to read Gaussian geometry data from input file; is the file formatted correctly?") 
        
        # We'll first try to split on comma (,) for charge,mult.
        charge_mult = geometry_section[0].split(",", 1)
        # If we didn't get what we want, try on whitespace.
        if len(charge_mult) != 2:
            charge_mult = geometry_section[0].split(" ", 1)
        
        # Now try and set.
        try:
            self.charge = int(charge_mult[0])
            self.multiplicity = int(charge_mult[1])
        except Exception:
            raise Digichem_exception("Unable to determine charge and multiplicity from '{}'".format(charge_mult))
        
        # The rest of the geometry section contains atoms.
        if len(geometry_section) != 2 or geometry_section[1] == "":
            raise Digichem_exception("Gaussian input file does not appear to contain any atoms")
        
        self.geometry = geometry_section[1]
        
        # Run a sanity check on the geometry format.
        for geometry_line in self.geometry.split("\n"):
            num_columns = len(geometry_line.split())
            if num_columns > 4:
                # We have too many columns?
                raise Digichem_exception("Gaussian input file has too many columns ({}) in its geometry section".format(num_columns))
                
        
        # And anything else.
        self.additional_sections = sections[3:]
    
    def geometry_with_iso(self, isotopes = {}):
        """
        Get the geometry of this gaussian input file with specific isotopes.

        :param isotopes: Isotope information. Each key should be an atomic index, or an element symbol, and each value the corresponding isotope (as an integer).
        """
        lines = []
        for index, line in enumerate(self.geometry.split("\n")):
            split_line = line.split()
            atom = split_line[0]
            if "(" in atom and ")" in atom:
                atom = atom[atom.find("("):]

            iso = None
            if index in isotopes:
                iso = isotopes[index]
            
            elif atom in isotopes:
                iso = isotopes[atom]
            
            if iso is not None:
                atom = "{}(iso={})".format(atom, iso)
            

            lines.append("    ".join(
                [atom, *split_line[1:]]
            ))
        
        return "\n".join(lines)

        
    @property
    def xyz(self):
        """
        Get the geometry of this gaussian input file in XYZ format.
        """
        lines = []
        for line in self.geometry.split("\n"):
            split_line = line.split()
            atom = split_line[0]
            if "(" in atom and ")" in atom:
                atom = atom[atom.find("("):]
            lines.append("    ".join(
                [atom, *split_line[1:]]
            ))
        return "{}\n\n{}".format(len(self.geometry.split("\n")), "\n".join(lines))

