import argparse


class Extend_action(argparse.Action):
    """
    An argparse action class to replace the "extend" action for versions of python (3.6 etc) that do not support it natively.
    """

    def __call__(self, parser, namespace, values, option_string = None):
        """
        """
        try:
            attr = getattr(namespace, self.dest)
            
        except AttributeError:
            attr = []
            setattr(namespace, self.dest, attr)
            
        attr.extend(values)

# TODO: Unused.
class List_grouper(argparse.Action):
    """
    Custom action class that groups lists together so we know in what order they were specified.
    """    
    
    def __init__(self, option_strings, *args, **kwargs):
        """
        """
        argparse.Action.__init__(self, option_strings, *args, **kwargs)
        # We'll give ourself a name so we also group no matter which of our option_strings is used.
        self.name = [name for name in option_strings if name[:2] == "--"]
    
    def __call__(self, parser, namespace, values, option_string=None):
        """
        """
        grouped_list = {
            'group': self,
            'values': values
        }
        try:
            getattr(namespace, self.dest).append(grouped_list)
        except AttributeError:
            setattr(namespace, self.dest, [grouped_list])