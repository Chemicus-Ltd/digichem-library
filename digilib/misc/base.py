from itertools import chain, combinations
import math

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def regular_range(median, number, spacing):
    """
    Generate a range of numbers with a given median and spacing.
    """
    if number % 2 == 0:
        # Even number of peaks
        peaks = []
        for magnitude in range(int(number /2)):
            peaks.append(median + (spacing /2 + magnitude * spacing))
            peaks.append(median - (spacing /2 + magnitude * spacing))
        
    else:
        # Odd number of peaks.
        peaks = [median]
        
        for magnitude in range(int((number -1) /2)):
            magnitude = magnitude+1
            peaks.append(median + magnitude * spacing)
            peaks.append(median - magnitude * spacing)
        
    return sorted(peaks)

def round_sig(number, sig_figs):
    """Round a number to a given number of significant figures."""
    try:
        return round(number, sig_figs-int(math.floor(math.log10(abs(number))))-1)

    except ValueError:
        if number == 0:
            return 0
        else:
            raise
    

def dict_list_index(dictionary, item):
    "Find the key in a dictionary which contains the list which contains an item."
    for dict_key, dict_value in dictionary.items():
        if item in dict_value:
            return dict_key
    
def dict_get(dict_obj, *fields):
    """
    Get a value from a nested dict of arbitrary depth.
    
    :param dict_obj: The nested dict.
    :param fileds: Names of nested keys.
    """
    if len(fields) == 1:
        return dict_obj[fields[0]]
    else:    
        return dict_get(dict_obj[fields[0]], *fields[1:])
    
def transpose(nested_list, dimensions):
    """
    Transpose a nested list, returning separate lists of the nested items.
    
    eg, transpose([(1, "a"), (2, "b"), (3, "c")])
    -> [[1, 2, 3], ["a", "b", "c"] 
    """
    if len(nested_list) == 0:
        return [list() for i in range(dimensions)]
    
    else:
        return list(map(list, zip(*nested_list)))
    