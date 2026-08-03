from deepmerge import Merger

from digilib.config.base import Digilib_options
from digilib.config.parse import Config_file_parser, Config_parser
from digilib.config.locations import master_config_path, system_config_location, user_config_location

from digilib.log import get_logger


# The main digichem options object.
# When running as a program, this will be merged with run-time options.
_options = None

def get_config(
        extra_config_files = None,
        extra_config_strings = None,
        cls = Digilib_options,
        clear_cache = False,
        sources = (
            master_config_path,
            system_config_location,
            user_config_location
        )
    ):
        """
        Get a Digichem options object.
        
        The returned object will take options from three sources:
        1) Any config files found in the default locations (see digilib.config.locations).
        2) Any additional config files specified as an argument.
        3) Any additional config strings specified as an argument.
        
        IMPORTANT: The returned options object will be merged with any additional options without prior copying.
        Hence the object returned by this function will be the same as the options attribute of this module.
        
        :param extra_config_files: An iterable of additional file paths to read from.
        :param extra_config_strings: An iterable of additional config options to parse.
        :param cls: The type of object to return.
        :param clear_cache: If True, and previously cached options will be discarded.
        :param sources: An iterable of file locations to read from.
        :return: A Digilib_options object (a fancy dict).
        """
        global _options
        if clear_cache:
            _options = None

        if _options is not None and extra_config_files is None and extra_config_strings is None:
            # Config has already been loaded (and we have nothing new to add).
            # Return the config object.
            return _options

        # Either this is the first time we've loaded the config (cache miss)
        # or we've been given extra options to add in.
        log_level = get_logger().level
        get_logger().setLevel("DEBUG")
        
        # First, load options if not already done so.
        if _options is None:
            # Load config options from given sources.
            # These objects are simple dicts.
            config = Config_file_parser(sources[0]).load(True)
            for source in sources[1:]:
                merge_dict(Config_file_parser(source).load(True), config)
                #  config.merge(Config_file_parser(source).load(True))
            
            # No need to validate here, we're going to do it later anyway.
            _options = cls(validate_now = False, **config)
        
        if extra_config_files is None:
            extra_config_files = []
            
        if extra_config_strings is None:
            extra_config_strings = []
        
        # Load any additional config files.
        for extra_config_file in extra_config_files:
            _options.deep_merge(Config_file_parser(extra_config_file).load())
            
        # Then load any additional config strings.
        for extra_config_string in extra_config_strings:
            _options.deep_merge(Config_parser(extra_config_string).load())
        
        # Check everything is valid.
        _options.validate()
        
        # Only restore log level if there was no problem.
        get_logger().setLevel(log_level)

        # And return.
        return _options

def merge_dict(new, old):
        """
        Recursively merge two dictionaries
         
        Any keys specified in new will overwrite those specified in old.
         
        :param new: A new dictionary to merge into the old.
        :param old: An old dictionary to be overwritten by new.
        """
        # Taken from deepmerge docs: https://deepmerge.readthedocs.io/en/latest/
        merger = Merger(
            # pass in a list of tuple, with the
            # strategies you are looking to apply
            # to each type.
            [
                (list, ["override"]),
                (dict, ["merge"]),
                (set, ["union"])
            ],
            # next, choose the fallback strategies,
            # applied to all other types:
            ["override"],
            # finally, choose the strategies in
            # the case where the types conflict:
            ["override"]
        )
        return merger.merge(old, new)