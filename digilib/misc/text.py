from itertools import groupby

def andjoin(listlike):
    """
    Joins a list of string-like variables with commas, except for the last item which is joined with 'and'.
    
    For example:
    >>> andjoin(["item1"])
    "item1"
    >>> andjoin(["item1", "item2"])
    "item1 and item2"
    >>> andjoin(["item1", "item2", "item3"])
    "item1, item2 and item3"
    """
    # Now combine.
    joined = ", ".join([str(item) for item in listlike[:-1]])
    
    if len(listlike) > 1:
        joined += " and " + str(listlike[-1])
    elif len(listlike) == 1:
        joined = listlike[0]
    
    return joined

def text_integer(number):
    """
    Convert an integer to a text form suitable for inline inclusion in the body of some text.
    
    If the integer is 0 <= x <= 9, the name of the number ("zero", "one", "two" etc) is returned, otherwise the value ("10", "11", "12" etc) is returned.
    """
    if number == 0:
        return "zero"
    elif number == 1:
        return "one"
    elif number == 2:
        return "two"
    elif number == 3:
        return "three"
    elif number == 4:
        return "four"
    elif number == 5:
        return "five"
    elif number == 6:
        return "six"
    elif number == 7:
        return "seven"
    elif number == 8:
        return "eight"
    elif number == 9:
        return "nine"
    else:
        return str(number)
    
def text_float(number, decimals = 2):
    """
    Convert a float to a text form suitable for inline inclusion in the body of some text.
    
    If the rounded number would appear to be 0.0, then the next smallest number that can be displayed with the given number of decimal points is returned instead, along with a 'less-than sing <'.
    For example:
    >>> text_float(0.567)
    "0.57"
    >>> text_float(0.001)
    "< 0.01"
    >>> text_float(0.007)
    "0.01"
    >>> text_float(0.0)
    "0.00"
    """
    if number == 0 or round(number, decimals) != 0:
        return "{:.{}f}".format(number, decimals)
    
    else:
        # The number is too small to diplay with this sig fig,
        # show that the number is not actually zero but is close.
        return "< {}".format(1 / 10 ** decimals)

    
def isare(listlike):
    if len(listlike) == 1:
        return "is"
    else:
        return "are"
    
def werewas(listlike):
    if len(listlike) == 1:
        return "was"
    else:
        return "were"
    
def ordinal_suffix(number):
    """
    Retuns the suffix of the ordinal (st, nd, rd or th) of a number.
    """
    last_digit = abs(int(number)) % 10
    last_two = abs(int(number)) % 100
    if last_two == 11 or last_two == 12 or last_two == 13:
        return "th"
    elif last_digit == 1:
        return "st"
    elif last_digit == 2:
        return "nd"
    elif last_digit == 3:
        return "rd"
    else:
        return "th"
    

def text_join(index, total, second_join = None):
    """
    """
    if index == total -1:
        # The last item, full stop.
        return "."
    elif second_join is not None and index == 0:
        return ", {}".format(second_join)
    elif index == total -2:
        # Second to last item, and.
        return " and"
    else:
        # Other item, comma.
        return ","


def listjoin(numbers):
    """
    
    """
    sorted_numbers = sorted(numbers)
    joined = []
    for keys, group in groupby(enumerate(sorted_numbers), lambda i: i[0]-i[1]):
        group = [number for index, number in group]        
        if len(group) == 1:
            joined.append("{}".format(group[0]))
        else:
            joined.append("{}-{}".format(group[0], group[-1]))
        
    return andjoin(joined)
    
    
