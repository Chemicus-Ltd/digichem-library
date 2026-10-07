# Miscellaneous functions for handling date and time.

# General imports.
import math

def latest_datetime(*datetimes):
    """
    Return the latest of an iterable of datetimes.
    
    :param datetimes:  Datetimes, of which the latest will be returned.
    :returns: The latest of the given datetimes, or None if no datetimes were given.
    """
    try:
        datetimes = list(datetimes)
        datetimes.sort()
        return datetimes[-1]
    except IndexError:
        # No datetimes.
        return None

def total_timedelta(*timedeltas):
    """
    Return the sum of multiple timedeltas.
    
    :param timedeltas: Timedelta objects to be summed.
    :returns: A new timedelta object that is the sum of the objects given, or None if no timedeltas were given.
    """
    try:
        return sum(timedeltas[1:], timedeltas[0])
    except IndexError:
        return None


def date_to_string(datetime_object):
    """
    Convert a datetime object to a standard string representation.
    
    :param datetime_object: The date and time to convert
    :return: A string representation of datetime_object.
    """
    return datetime_object.strftime("%d/%m/%Y %H:%M:%S")

def timedelta_to_string(timedelta_object):
    """
    Convert a timedelta object to a standard string representation.
    
    :param timedelta_object: The time difference to convert
    :return: A string representation of timedelta_object.
    """
    hours = math.floor(timedelta_object.seconds / 3600)
    minutes = math.floor((timedelta_object.seconds - hours * 3600) / 60)
    # Seconds lacks microsecond accuracy, but we can add it on.
    seconds = round((timedelta_object.seconds + timedelta_object.microseconds /1000000) - (hours * 3600 + minutes * 60))
    
    date_str = ""
    add = False
    
    # Only add each section if not zero.
    if timedelta_object.days != 0:
        date_str += "{} d, ".format(timedelta_object.days)
        add = True
        
    if hours != 0 or add:
        date_str += "{} h, ".format(hours)
        add = True
        
    if minutes != 0 or add:
        date_str += "{} m, ".format(minutes)
        add = True
        
    date_str += "{} s".format(seconds)
    
    return date_str