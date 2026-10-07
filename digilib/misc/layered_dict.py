from collections import OrderedDict

class Layered_dict(list):
    """
    An enhanced OrderedDict object that remembers layers/groups of positions.
    
    A layered dict is essentially a wrapper around a list of ordered dicts. They expose methods found both in list and in dict. Most methods are similar to those found in list (including setting an index), but iteration operates over each dict in turn.
    
    New dicts can either be updated to the last OrderedDict (the same as updating a normal OrderedDict object).
    Or appended to the end of the list, creating a new layer/group.
    
    Layered_dicts are useful for preserving order in cases where you might not know the dict keys in some instances of a layer. Layered_dicts can later be normalised; comparing the order of all keys from all layers (and intelligently skipping layers that contain None) to produce a normal OrderedDict.
    
    Layered_dicts can also be nested to support even more complex cases.
    
    Note that Layered_dict extends from list (and not from dict, at least not yet), be careful if type checking (isinstance(my_layered_dict, dict) == False ...).
    """
    
    def __init__(self, *iterables, **key_values):
        """
        Constructor for Layered_dict objects.
        
        Layered_dicts can be initialised from a single iterable in the same manor as normal dicts. Multiple iterables can also be given as multiple positional args; each iterable will be stored as its own layer. Alternatively, key/value pairs can be given directly as **keyword args (but note that this form is incompatible with the positional args form because the ordering of the layers would be ambiguous.
        
        If no args are given, an empty Layered_dict is created.
        
        :param *iterables: An optional number of iterables to add as new layers.
        :param **key_values: Alternatively, an optional number of keyword arguments which are taken to be key=value pairs to construct a single layer of dict.
        """
        # Init our listy self.
        super().__init__()
        
        if len(iterables) > 0 and len(key_values) > 0:
            raise ValueError("The positional and keyword argument forms are mutually exclusive")
        elif len(iterables) > 0:
            self.extend(iterables)
        elif len(key_values) > 0:
            self.append(key_values)
        else:
            # We'll always have at least one OrderedDict object.    
            self.append({})
    
    @property
    def active_dict(self):
        """
        Return the current active dict object.
        
        The active dict is merely the last dict appended, so use len() to get its index. This method exists only for convenience. 
        """
        return self[-1]
    
    def __setitem__(self, key, value):
        """
        Set the value of one of the layers in this Layered_dict.
        
        Layered_dicts can only contain OrderedDicts, Layered_dict and None as values. Any other values will first be used to construct a new OrderedDict object.
        """
        # First do some type checking.
        if not isinstance(value, (OrderedDict, Layered_dict, type(None))):
            value = OrderedDict(value)
        
        # Now insert.
        list.__setitem__(self, key, value)
        
        # If we inserted None, insert a new LayeredDict too.
        if value is None:
            self.append(OrderedDict)
        
    def __delitem__(self, *args, **kwargs):
        """
        Delete one of the layers of this Layered_dict.
        
        Note that it is impossible to delete all layers; a new empty OrderedDict layer will be created if no layers remain.
        """
        list.__delitem__(self, *args, **kwargs)
        if len(self) == 0:
            self[0] = OrderedDict()
            
    def get(self, layer_index, key, key_index = None):
        """
        Get a value from this Layered_dict.
        
        :raises IndexError: If the one of the given layer indices is not contained in this Layer_dict.
        :raises KeyError: If the given key is not found in the dict at the specified level.
        :param layer_index: A tuple of indices indicating the layer to get from. See items() for further explanation.
        :param key: The name of the key to get at the indicated layer.
        :param key_index: An optional constraint; if given, a KeyError exception is raised if the requested key is not at the the index given by key_index.
        """
        # If layer_index contains multiple parts, we can recurse.
        if len(layer_index) > 1:
            index = layer_index[0]
            if isinstance(self[index], type(self)):
                return self[index].get(layer_index[1:], key, key_index = key_index)
            else:
                # We'be been given a multi-depth layer_index, but our object at that position is not another Layered_dict.
                raise IndexError("Bad layer_index, object ({}) at index {} is not of type {}".format(type(self[index]).__name__, index, type(self).__name__))
            
        # We only have one index to worry about.
        index = layer_index[0]
        
        # Get the key from the dict.
        # But first check to make sure our object at index is not another Layered_dict.
        if isinstance(self[index], type(self)):
            # It is a Layered_dict, raise IndexError.
            raise IndexError("Bad layer_index, object ({}) at index {} is of type {}".format(type(self[index]).__name__, index, type(self).__name__))
        
        # Also check to make sure the dict isn't None.
        elif self[index] is None:
            raise KeyError("Cannot retrieve key '{}' at layer_index '{}'; dict is None".format(key, index))
        
        # Now check the given key_index is the same (if a key_index was given).
        if key_index is not None and list(self[index].keys()).index(key) != key_index:
            raise KeyError("Key/index mismatch; the key '{}' is not at requested index '{}'".format(key, key_index))
    
        return self[index][key]
        
        
    def items(self):
        """
        Generator to support iteration over this Layered_dict.
        
        Each iteration over a Layered_dict returns a 3-membered tuple of the form (layer_index, key_name, value) where:
            layer_index is  a tuple of indices of the current layer being iterated through. For example, the second item of the first layer would be returned as (1,). The third item of the second layer would be (1,2). Each item in layer_index is an index to a nested layer, so the length of layer_index is variable and depends on the degree of nesting. A Layered_dict that contains no other Layered_dict objects will have len(layer_index) == 1 for all items.
            key_name is the name of the key. Note that unlike real dicts, the same key can appear multiple times in a Layered_dict.
            value is the value of the current key in the current dict being iterated through.
            
        Although layered dicts can contain None as placeholder values, iteration will skip these automatically.
        """
        for layer_index, dict_layer in enumerate(self):
            # Some layers can be None place holders; skip these.
            if dict_layer is None:
                continue
            
            # Other layers can be nested Layered_dict objects.
            elif isinstance(dict_layer, Layered_dict):
                for sub_path, key, value in dict_layer.items():
                    # Add our child's index hierarchy to ours.
                    index_path = [layer_index]
                    index_path.extend(sub_path)
                    
                    yield (tuple(index_path), key, value)
            
            # The remaining layers are normal dicts.
            else:
                for key in dict_layer:
                    yield ((layer_index,), key, dict_layer[key])
            
        #while True:
        #    # Nothing more to iterate over.
        #    raise StopIteration()
        
    def update(self, *iterables, **key_values):
        """
        Update the current active dict with new values.
        
        The order of the keys is preserved as for normal OrderedDicts (the first insertion order is preserved). 
        
        :param *iterables: An optional number of iterables to update with. Each iterable will be traversed in the order given and used to update the same dictionary.
        :param **key_values: Alternatively, an optional number of keyword arguments which are used to update the active dict directly.
        :return: This Layered_dict object, for convenience.
        """
        if len(iterables) > 0 and len(key_values) > 0:
            raise ValueError("The positional and keyword argument forms are mutually exclusive")
        elif len(iterables) > 0:
            # Positional arg form.
            # Add each iterable in turn.
            for iterable in iterables:
                if isinstance(iterable, type(self)):
                    # The iterable is another Layered_dict, update with flatten().
                    self.active_dict.update(iterable.flatten())
                else:
                    self.active_dict.update(iterable)
        else:
            # Keyword arg form, update from given keywords.
            self.active_dict.update(**key_values)
            
        return self
        
    def flatten(self):
        """
        Flatten this Layered_dict, returning a single OrderedDict object.
        
        The order of keys in the returned dict will honour the ordering of this Layered_dict object.
        If the same key is present more than once, the first instance determines the order the key appears in, but the last determines its value.
        
        Note that unlike iterating directly over this Layered_dict, the same key cannot appear more than once in the returned flat dict.
        """
        # The implementation here is simple, just update a new dict from each dict we contain.
        return OrderedDict({key:value for index, key, value in self.items()})
    
    def normalise(self, mapping):
        """
        Normalise this Layered_dict according to a mapping so that is aligns with a number of other Layered_dict.
        
        A mapping describes the order in which key/values should be aligned. Although the ordering of keys is defined for each individual Layered_dict, this ordering can be ambiguous if two Layered_dicts contain different keys. Similarly, an individual Layered_dict cannot properly align None placeholder keys because it cannot know how many keys are expected. A mapping resolves this ambiguity.
        
        A mapping is a list of tuples of the form (layer_index, key_name), where:
            layer_index has the same meaning as in items().
            #key_index is the index of the key in dict to get. It is required to differentiate between duplicate keys with different orderings.
            key_name is the name of the key to get (same as in items()).
        
        :param mapping: A list of tuples of the form (layer_index, key_name).
        :return: A list of values as dictated by mapping. Each index of the returned list is controlled by the corresponding index of mapping (and so the two lists will be of the same length). If the index in mapping is None or specifies a key not contained in this dict, None is inserted into the list.
        """
        normal_list = []
        
        for mapped_index in mapping:
            # If the map is None, so is the value we insert.
            mapped_value = None
            
            # Try and get a value from ourself.
            if mapped_index is not None:
                try:
                    mapped_value = self.get(mapped_index[0], mapped_index[1])
                except (IndexError, KeyError):
                    mapped_value = None
            
            # Add to the list.
            normal_list.append(mapped_value)
                
        # And return.
        return normal_list
        
    @classmethod
    def common_mapping(self, layered_dicts):
        """
        Generate a mapping from a list of Layered_dict objects, suitable for passing to normalise().
        
        A mapping is a list of tuples of the form (layer_index, key_name), with the same meaning as in normalise(). It is used to align the keys of multiple dicts.
        
        The mapping generated will contain an entry for each unique (layer_index, key_name) in all of the given layered_dicts. The order of the mapping is well defined, respecting the order of layer_index in the first instance, then the order of keys as found in the respective OrderedDict.
        
        The order is poorly defined, however, for non-identical keys at the same layer_index (because the ordering is ambiguous).
        
        :param layered_dicts: An iterable of Layered_dict objects to generate a mapping from.
        :return: A list describing the mapping.
        """
        # To get our mapping, we'll add non duplicate flattened forms of each layered_dict to our list.
        mapping = []
        for layered_dict in layered_dicts:
            # We can discard the value part as we're only interested in the order for now.
            for index, key, value in layered_dict.items():
                # Add this index_key pair to mapping if not already present.
                # In addition to the normal layer_index, we'll also include the key_index (to support multiple identical dicts with different key ordering).
                if (index, key) not in mapping:
                    mapping.append((index, key))
        
        
        # Finally, we partially re-sort in case any layer indices got added in the wrong order.
        # We won't re-order key names though, so these will still respect their relative orders from their OrderedDicts.
        # If two layer index pairs have OrderedDicts with different keys, the order in which these two OrderedDicts will be combined is undefined (but the keys within them will remember their relative order).
        mapping.sort(key = lambda index_key: index_key[0])
        
        # Done.
        return mapping
    
    @classmethod
    def tabulate(self, layered_dicts):
        """
        Tabulate a list of layered_dicts.
        
        tabulate() is a convenience function, first calling common_mapping() on the list of layered_dicts followed by normalise().
        
        :param layered_dicts: A list of layered_dict objects to tabulate. The order of rows will be the same as the order of this list. Alternatively, a number of OrderedDicts can be included in the list, which will automatically be converted to Layered_dicts.
        :return: A tuple of (headers, table), where headers is a list of table headers (derived from the names of the keys in layered_dicts) and table is the tabulated data (each row will be of the same length, possible padded with None).
        """
        # First convert any OrderedDicts
        layered_dicts = [self(current_dict) if isinstance(current_dict, OrderedDict) else current_dict for current_dict in layered_dicts]
        
        # First, get our mapping.
        mapping = self.common_mapping(layered_dicts)
        
        # Now get our list.
        table = [layered_dict.normalise(mapping) for layered_dict in layered_dicts]
        
        # Return both.
        return ([header[1] for header in mapping], table)
        
        
        
        
        
        
        
        