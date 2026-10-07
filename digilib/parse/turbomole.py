# General imports.
import re
from datetime import timedelta, datetime
import glob, pathlib
import warnings

# Digichem imports.
from digilib.parse.cclib import Cclib_parser

# Hidden imports.
#from cclib.io.ccio import sort_turbomole_outputs


class Turbomole_parser(Cclib_parser):
    """
    Top level class for parsing output from Turbomole files.
    """
    DAYS_REGEX = re.compile(r"([0-9.]*) days")
    HOURS_REGEX = re.compile(r"([0-9.]*) hours")
    MINUTES_REGEX = re.compile(r"([0-9.]*) minutes")
    SECONDS_REGEX = re.compile(r"([0-9.]*) seconds")
    
    @classmethod
    def sort_log_files(self, log_files):
        """
        Sort a list of log files into a particular order, if required for this parser.
        """
        from cclib.parser.turbomoleparser import Turbomole
        
        return Turbomole.sort_input(log_files)
    
    def duration_to_timedelta(self, duration_str):
        """
        Convert a Turbomole duration string into an equivalent timedelta object.
        """
        time_parts = {'days': 0, 'hours': 0, 'minutes': 0, 'seconds': 0}
        
        for time_part in time_parts:
            # Use regex to look for each part in the string.
            match = getattr(self, time_part.upper() + '_REGEX').search(duration_str)
            if match:
                time_parts[time_part] = float(match.group(1))
                
        # Build a timedelta from our parts.
        duration = timedelta(days = time_parts['days'], hours = time_parts['hours'], minutes = time_parts['minutes'], milliseconds = time_parts['seconds'] * 1000)
        
        # All done.
        return duration    
    
    def pre_parse(self):
        """
        Perform any setup before line-by-line parsing.
        """
        super().pre_parse()
        # Look for duration information.
        # Only bother doing this if we don't have timings from cclib.
        if 'wall_time' not in self.data.metadata:
            self.data.metadata['wall_time'] = []
        
        if 'cpu_time' not in self.data.metadata:
            self.data.metadata['cpu_time'] = []
    
    def parse_output_line(self, log_file, line):
        """
        Perform custom line-by-line parsing of an output file.
        """
        # Only bother doing this if we don't have timings from cclib.
        if 'wall_time' not in self.data.metadata or 'cpu_time' not in self.data.metadata:
            # Look for duration.
            if "total  cpu-time :" in line:
                self.data.metadata['cpu_time'].append(self.duration_to_timedelta(line))
                
            elif "total wall-time :" in line:
                self.data.metadata['wall_time'].append(self.duration_to_timedelta(line))
            
        # And also end date.
        if ": all done  ****" in line:
            # Skip 2 lines.
            line = next(log_file)
            line = next(log_file)
            line = next(log_file)
            try:
                self.data.metadata['date'] = datetime.strptime(line.strip(), "%Y-%m-%d %H:%M:%S.%f").timestamp()
            except ValueError:
                # We couldn't parse.
                pass
    
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """
        super().post_parse()
            
        # Delete our durations if they are zero.
        if len(self.data.metadata['wall_time']) == 0:
            del(self.data.metadata['wall_time'])
            
        if len(self.data.metadata['cpu_time']) == 0:
            del(self.data.metadata['cpu_time'])

    def process(self):
        """
        Get a Result set object from this parser.
        
        :return: The populated result set.
        """
        super().process()
        
        # After processing is complete, have a look for excited state density files.
        # These have the general file name:
        # adcp2-xsdn-1a-001-total.cao
        #  ^    ^     ^  ^------------ Excited state number (1, 2, 3 etc).
        #  |    |     ---------------- Irrep (multiplicity and symmetry)
        #  |    ---------------------- Excited state
        #  --------------------------- Method (MP2, ADC(2), CC2). 
        #
        # Each found density file will be stored under auxiliary_files['excited_state_cao_files'][state_symbol] where state_symbol is the corresponding state, eg S(1).
        # If we have no excited states we can skip this altogether.
        if len(self.results.excited_states) != 0:
            for log_file in self.log_file_paths:
                excited_densities = {}
                state_num = 1
                # Look for each numbered excited state until we run out of files.
                while True:
                    found_densities = glob.glob(str(pathlib.Path(log_file.parent, "*xsdn*{:03}*total*.cao".format(state_num))))
                    
                    # Sort results in case there is more than one result the same behaviour is encountered between multiple runs.
                    found_densities.sort()
                    
                    if len(found_densities) > 0:
                        try:
                            # Get the state that corresponds to this file.
                            excited_state = self.results.excited_states[state_num -1]
                            excited_densities[excited_state.state_symbol] = pathlib.Path(found_densities[0])
                            
                            # Print a warning if there's more than one (because we don't know how to handle that scenario).
                            if len(found_densities) > 1:
                                warnings.warn("Found multiple excited state density files for state '{}' in Turbomole calculation directory '{}'; using file '{}' and ignoring '{}'".format(excited_state.state_symbol, log_file.parent, pathlib.Path(excited_densities[excited_state.state_symbol]).name, ", ".join([pathlib.Path(density).name for density in found_densities[1:]])))
                        
                        except IndexError:
                            warnings.warn("Could not find excited state that corresponds to excited state density file '{}' with index {}; this file will be ignored".format(found_densities[0], state_num -1))
                
                    if len(found_densities) == 0:
                        # All done.
                        break
                    
                    state_num += 1
                    
                # Update auxiliary_files with new files.
                if len(excited_densities) > 0:
                    try:
                        self.results.metadata.auxiliary_files['excited_state_cao_files'].update(excited_densities)
                        
                    except KeyError:
                        if 'excited_state_cao_files' not in self.results.metadata.auxiliary_files:
                            self.results.metadata.auxiliary_files['excited_state_cao_files'] = excited_densities
                        else:
                            raise
        
        return self.results
            
    @classmethod
    def find_auxiliary_files(self, hint, basename):
        """
        Find auxiliary files from a given hint.
        
        :param hint: A path to a file to use as a hint to find additional files.
        :returns: A dictionary of found aux files.
        """
        auxiliary_files = super().find_auxiliary_files(hint, basename)
        
        parent = pathlib.Path(hint).parent if not hint.is_dir() else hint
            
        # Find .cao density files.
        # First look for ground state density files.
        # These have the general file name:
        # mp2-gsdn-1a-000-total.cao
        #  ^    ^   ^----------------- Irrep (multiplicity and symmetry)
        #  |    ---------------------- Ground state
        #  --------------------------- Method (MP2, ADC(2), CC2). 
        #
        ground_densities = glob.glob(str(pathlib.Path(parent, "*gsdn*total*.cao")))
        
        # Sort results in case there is more than one result the same behaviour is encountered between multiple runs.
        ground_densities.sort()
        
        if len(ground_densities) > 0:
            auxiliary_files['ground_state_cao_file'] = pathlib.Path(ground_densities[0])
        
        # Print a warning if there's more than one (because we don't know how to handle that scenario).
        if len(ground_densities) > 1:
            warnings.warn("Found multiple ground state density files in Turbomole calculation directory '{}'; using file '{}' and ignoring '{}'".format(parent, pathlib.Path(auxiliary_files['ground_state_cao_file']).name, ", ".join([pathlib.Path(density).name for density in ground_densities[1:]])))
        
        #################################################################################################
        # Note that excited state densities are also located, but this is done in the process() method. #
        #################################################################################################
        
        return auxiliary_files

    