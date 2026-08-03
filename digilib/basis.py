import deepmerge

# Hidden imports.
# import basis_set_exchange.misc
# import basis_set_exchange.writers

class BSE_basis_set(dict):
    """
    A class for representing a (number of) basis sets fetched from the basis set exchange.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Constructor for BSE_basis_set objects.
        
        :param definition: A dictionary where each key is the name of a basis set and the value is the elements it applies to. The values 'all', '*', and None can be used to indicate the basis set applies to all atoms.
        """
        super().__init__( *args, **kwargs)
        
    @classmethod
    def is_filter_all(self, element_filter):
        """
        Determine whether a given element filter implicitly includes all elements.
        """
        return element_filter == "*" or element_filter == "all" or element_filter is None
    
    def __str__(self):
        return self.name
    
    @property
    def name(self):
        """
        Get a descriptive name of this basis set.
        """
        import basis_set_exchange.misc
        
        names = []
        
        for basis_set_name, basis_set_elements in self.items():
            name = basis_set_name
            if not self.is_filter_all(basis_set_elements):
                name += " []".format(basis_set_exchange.misc.compact_elements(basis_set_exchange.misc.expand_elements(basis_set_elements)))
                
        return ", ".join(names)
    
    def has_ECPs(self, elements_filter = None):
        """
        Determine whether this basis set would provide effective core potentials (ECPs) when applied to a list of elements.
        
        :param elements_filter: The elements to check against.
        """
        basis_sets = self.to_dict(elements_filter)
        
        return 'scalar_ecp' in basis_sets['function_types']
    
    def to_dict(self, elements_filter = None):
        """
        Convert the basis set information represented by this object to a basis set dict.
        
        :param elements_filter: A list/string of elements to filter by. Only elements given here will be printed in the final format. Each item can be an int or str representing a single element (1, '1', 'H' etc), or a range of elements ('1-5' etc).
        """
        import basis_set_exchange.misc
        
        if isinstance(elements_filter, int):
            elements_filter = str(elements_filter)

        if elements_filter is not None and not isinstance(elements_filter, str):
            elements_filter = ",".join((str(item) for item in elements_filter))

        if elements_filter is not None:
            elements_filter = basis_set_exchange.misc.expand_elements(elements_filter)
        
        basis_sets = {}
        
        # First, convert each definition to to a basis_set_exchange dict
        for basis_set_name, basis_set_elements in self.items():
            # Only elements that are actually present will have a basis set recorded.
            # The values of 'all', '*' and None are also accepted as meaning all atoms.
            if not self.is_filter_all(basis_set_elements):
                basis_set_elements = basis_set_exchange.misc.expand_elements(basis_set_elements)
                if elements_filter is not None:
                    elements = list(set(basis_set_elements).intersection(elements_filter))
                    
                else:
                    elements = basis_set_elements
            
            else:
                # We are using all elements from this basis set.
                elements = elements_filter
            
            if elements is None or len(elements) > 0:
                # Get the basis set.
                basis_set = basis_set_exchange.get_basis(basis_set_name, elements)
                
                # Merge it with our total basis set.
                # TODO: Check this is actually safe?
                deepmerge.always_merger.merge(basis_sets, basis_set)
            
        return basis_sets
    
    def to_format(self, fmt = None, elements_filter = None):
        """
        Convert the basis set information represented by this object to a basis set format (for example, one suitable for a CC program).
        
        :param fmt: The format to write to, see basis_set_exchange.get_formats().
        :param elements_filter: A list of elements to filter by. Only elements given here will be printed in the final format. Each item can be an int or str representing a single element (1, '1', 'H' etc), or a range of elements ('1-5' etc).
        :returns: The formatted basis set.
        """
        import basis_set_exchange.writers
        
        basis_sets = self.to_dict(elements_filter)
            
        # Output the merged basis set.
        if fmt is not None:
            return basis_set_exchange.writers.write_formatted_basis_str(basis_sets, fmt = fmt)
        
        else:
            return basis_sets
    