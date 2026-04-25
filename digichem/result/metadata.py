"""Code for metadata regarding finished calculation results."""

# General imports.
from datetime import timedelta, datetime
import math
import itertools
from deepmerge import conservative_merger
from pathlib import Path
import copy
import warnings
from scipy import integrate

from digichem.misc.time import latest_datetime, total_timedelta, date_to_string,\
    timedelta_to_string
from digichem.misc.text import andjoin
from digichem.exception import Result_unavailable_error
from digichem.result import Result_object
import digichem
from digichem import translate
from digichem.memory import Memory


class Solvent(Result_object):
    """
    Class for storing solvent metadata.
    """
    
    def __init__(self, model = None, name = None, params = None):
        self.model = model
        self.name = name
        # A raw dictionary of solvent related parameters.
        self.params = params if params is not None else {}
        
        # If we are missing only one of name and params['epsilon'], look up the missing value.
        if self.name is None and "epsilon" in self.params:
            try:
                self.name = translate.Solvent.epsilon_to_name(self.params['epsilon'])
            
            except ValueError:
                # Nothing close.
                pass
        
        if self.name is not None:    
            self.name = self.name.capitalize()
        
    @classmethod
    def from_parser(self, parser):
        """
        Construct a Solvent object from an output file parser.
        
        :param parser: Output data parser.
        :return: A populated Solvent object.
        """
        return self(
            name = parser.data.metadata.get('solvent_name', None),
            model = parser.data.metadata.get('solvent_model', None),
            params = parser.data.metadata.get('solvent_params', {}),
        )
        
    @property
    def description(self):
        """
        A common-name description of this solvent.
        """
        # If we have no solvent, say so.
        if self.model is None:
            return "Gas-phase"
        
        # If we have a name, just use that.
        elif self.name is not None:
            return self.name
        
        # If we only have epsilon, use that.
        elif "epsilon" in self.params:
            return "ε = {}".format(self.params['epsilon'])
        
        # Give up.
        else:
            return "Unknown"
        
    def _dump_(self, digichem_options, all):
        return {
            "model": self.model,
            "name": self.name,
            "params": self.params
        }
        
    @classmethod
    def merge(self, *multiple_objects):
        """
        Merge multiple solvent implementations together.
        """
        # Check all items are the same.
        if not all(obj == multiple_objects[0] for obj in multiple_objects if obj is not None):
            warnings.warn("Refusing to merge different solvent methods")
            return self()
            
        return multiple_objects[0]
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(
            model = data.get('model', None),
            name = data.get('name', None),
            params = data.get('params', {})
        )
        
    def __eq__(self, other):
        """Is this solvent implementation the same as another one?"""
        return (
            # Might want to do this case-insensitive.
            self.model == other.model and
            abs(self.params.get("epsilon", math.inf) - self.params.get("epsilon", math.inf)) < 0.001
        )


class Metadata(Result_object):
    """
    Class for storing calculation metadata.
    """
    
    # A dictionary of known computational orders, where the key is name of the method (upper case) and the value is an integer describing the ordered 'level of theory'.
    METHODS = {
        "HF":       1,
        "DFT":      2,
        "LMP2":     3,
        "DF-MP2":   4,
        "MP2":      5,
        "MP3":      6,
        "MP4":      7,
        "MP5":      8,
        "CCSD":     9,
        "CCSD(T)":  10,
        "CCSD-T":   11
    }
    
    def __init__(
            self,
            jobId = None,
            name = None,
            user = None,
            userID = None,
            log_files = None,
            auxiliary_files = None,
            history = None,
            date = None,
            insert_date = None,
            duration = None,
            package = None,
            package_version = None,
            calculations = None,
            success = None,
            methods = None,
            functional = None,
            basis_set = None,
            charge = None,
            multiplicity = None,
            optimisation_converged = None,
            temperature = None,
            pressure = None,
            orbital_spin_type = None,
            digichem_version = None,
            silico_version = None,
            solvent = None,
            num_cpu = None,
            memory_available = None,
            memory_used = None,
            performance = None,
            
            # Deprecated.
            solvent_model = None,
            solvent_name = None,):
        """
        Constructor for result Metadata objects.
        
        :param jobId: If this result was generated from a digichem calculation, the relevant jobID.
        :param name: Optional name of this calculation result.
        :param user: The username of the user who parsed this result.
        :param userID: The unique ID of the user who parsed this result.
        :param log_files: An optional list of text-based calculation log files from which this result was parsed.
        :param auxiliary_files: An optional dict of auxiliary files associated with this calculation result.
        :param history: Optional SHA of the calculation from which the coordinates of this calculation were generated.
        :param num_calculations: Optional number of individual calculations this metadata represents.
        :param date: Optional date (datetime object) of this calculation result.
        :param insert_date: Optional date (datetime object) of when this calculation result was stored (normally in a DB).
        :param duration: Optional duration (timedelta object) of this calculation.
        :param package: Optional string identifying the computational chem program that performed the calculation.
        :param package_version: Optional string identifying the version of the computational chem program that performed the calculation.
        :param calculations: A list of strings of the different calculations carried out (Opt, Freq, TD, TDA, SP etc).
        :param success: Whether the calculation completed successfully or not.
        :param methods: List of methods (DFT, HF, MPn etc) used in the calculation.
        :param functional: Functional used in the calculation.
        :param basis_set: Basis set used in the calculation.
        :param charge: Charge (positive or negative integer) of the system studied.
        :param multiplicity: The multiplicity of the system.
        :param optimisation_converged: Whether the optimisation converged or not.
        :param temperature: The temperature used in the calculation (not always relevant).
        :param pressure: The pressure used in the calculation (not always relevant).
        :param orbital_spin_type: The types of orbitals that have been calculated, either 'restricted' or 'unrestricted'.
        :param num_cpu: The number of CPUs used to perform this calculation.
        :param memory_available: The maximum amount of memory available to this calculation (the amount requested by the user).
        :param memory_used: The maximum amount of memory used by the calculation (the amount requested by the user).
        """
        self.jobId = jobId
        self.num_calculations = 1
        self.name = name
        self.user = user
        self.userID = userID
        self.log_files = log_files if log_files is not None else []
        self.auxiliary_files = auxiliary_files if auxiliary_files is not None and len(auxiliary_files) != 0 else {}
        self.history = history
        self.date = date
        self.insert_date = insert_date
        self.duration = duration
        self.package = package
        self.package_version = package_version
        self.calculations = calculations if calculations is not None else []
        self.success = success
        self.methods = methods if methods is not None else []
        self.functional = functional
        self.basis_set = basis_set
        # TODO: charge and mult should be deprecated here, they are available in ground_state.
        self.charge = charge
        self.multiplicity = multiplicity
        self.optimisation_converged = optimisation_converged
        self.temperature = temperature
        self.pressure = pressure
        self.orbital_spin_type = orbital_spin_type
        # TOOD: Ideally this would be parsed from the calculation output somehow, but this is fine for now.
        if silico_version:
            self.digichem_version = silico_version
        self.digichem_version = digichem.__version__ if digichem_version is None else digichem_version
        self.solvent = solvent
        self.num_cpu = num_cpu
        self.memory_available = memory_available
        self.memory_used = memory_used
        self.performance = performance
        
        # Deprecated solvent system.
        if solvent_model is not None:
            self.solvent.model = solvent_model
        
        if solvent_name is not None:
            self.solvent.name = solvent_name
    
    # TODO: This is more than a bit clumsy and in general the handling of names should be improved.
    @property
    def molecule_name(self):
        return Path(self.name).name if self.name is not None else None
    
    @property
    def level_of_theory(self):
        """
        A short-hand summary of the methods and basis sets used.
        """
        theories = []
        if len(self.converted_methods) > 0:
            #theories.extend(self.converted_methods)
            theories.append(self.converted_methods[-1])
            
        if self.basis_set is not None:
            theories.append(self.basis_set)
            
        return("/".join(theories))
        
    @classmethod
    def merge(self, *multiple_metadatas):
        """
        Merge multiple metadata objects into a single metadata.
        """
        return Merged_metadata.from_metadatas(*multiple_metadatas)
    
    @property
    def converted_methods(self):
        """
        Similar to the methods attribute but where DFT is replaced with the actual functional used.
        """
        return [self.functional if method == "DFT" and self.functional is not None else method for method in self.methods if method is not None]
    
    @property
    def description(self):
        desc = []
        if self.name is not None:
            desc.append(self.name)
        
        desc.append(self.identity_string)
        
        return ", ".join(desc)
        
    @property
    def identity(self):
        """
        A dictionary of critical attributes that identify a calculation.
        """
        # Get our list of methods, replacing 'DFT' with the functional used if available.
        return {
            "package": self.package,
            "calculations": self.calculations_string,
            "methods": ", ".join(self.converted_methods) if len(self.converted_methods) > 0 else None,
            "basis": self.basis_set,
        }
        
    @property
    def identity_string(self):
        """
        A string that identifies this calculation.
        """
        return " ".join([value for value in self.identity.values() if value is not None])
        
    
    @classmethod
    def sorted_methods(self, methods):
        """
        Order a list of methods (HF, DFT, MP2 etc) in terms of 'level of theory', with the lowest level (HF or DFT probably) first and highest (MP4, CCSD etc) last.
        """
        return sorted((method.upper() for method in methods), key = lambda method: self.METHODS[method] if method in self.METHODS else math.inf)
    
    @property
    def package_string(self):
        """
        The comp chem package and version combined into one string.
        """
        package_string = self.package if self.package is not None else ""
        
        # Add version string if we have it.
        package_string += " " if package_string != "" else ""
        package_string += "({})".format(self.package_version)
        
        # Done.
        return package_string 
        
    @property
    def calculations_string(self):
        """
        Get the list of calculation types as a single string, or None if there are no calculations set.
        """
        return ", ".join(self.calculations) if len(self.calculations) != 0 else None
    
    @property
    def methods_string(self):
        """
        Get the list of calculation methods as a single string, or None if there are no methods set.
        """
        return ", ".join(self.methods) if len(self.methods) != 0 else None
    
    @property
    def human_calculations_string(self):
        """
        A version of calculations_string that is more pleasant to read.
        """
        calculations = []
        if "Single Point" in self.calculations:
            calculations.append("single point energy")

        if "Screening" in self.calculations:
            calculations.append("conformer screening")
        
        if "Optimisation" in self.calculations:
            calculations.append("optimised structure")
            
        if "Frequencies" in self.calculations:
            calculations.append("vibrational frequencies")
            
        if "Excited States" in self.calculations:
            calculations.append("excited states")
            
        if "NMR" in self.calculations:
            calculations.append("NMR properties")
            
        if len(calculations) == 0:
            # Use a generic term.
            calculations.append("properties")
        
        return andjoin(calculations)
            
    @classmethod
    def get_calculation_types_from_cclib(self, ccdata):
        """
        Determine what types of calculation we did from the data provided by cclib.
        
        :param ccdata: The data from cclib.
        :return: A list of strings.
        """
        calcs = []
        # Need to think of a better way to determine what is an optimisation, as frequency calcs contain an optdone for some reason.
        #if hasattr(ccdata, 'optdone'):
        if len(getattr(ccdata, 'scfenergies', [])) > 1:
            calcs.append('Optimisation')
        if hasattr(ccdata, 'vibfreqs'):
            calcs.append('Frequencies')
        if hasattr(ccdata, 'etenergies'):
            calcs.append('Excited States')
        if hasattr(ccdata, 'nmrtensors'):
            calcs.append('NMR')
        # If our list is empty, assume we did an SP.
        if len(calcs) == 0 and ( hasattr(ccdata, 'scfenergies') or hasattr(ccdata, 'mpenergies') or hasattr(ccdata, 'ccenergies') ):
                calcs = ['Single Point']
        # Return our list.
        return calcs
    
    @classmethod
    def get_methods_from_cclib(self, ccdata):
        """
        Get a list of method types (DFT, MP etc.) from the data provided by cclib.
        
        :param ccdata: The data from cclib.
        :return: A list of strings of the methods used (each method appears once only; the order has no special significance). The list may be empty.
        """
        return self.sorted_methods((set(ccdata.metadata.get('methods', []))))
        #return list(dict.fromkeys(ccdata.metadata.get('methods', [])))
    
    @classmethod
    def get_orbital_spin_type_from_cclib(self, ccdata):
        """
        Determine the orbital type (restricted, unrestricted) from the data provided by cclib.
        
        :return ccdata: The data from cclib.
        :return: A string describing the orbital type ('restricted', 'unrestricted' or 'none').
        """
        #mo_len = len(ccdata.get('moenergies', []))
        mo_len = len(getattr(ccdata, 'moenergies', []))
        if mo_len == 1:
            return "restricted"
        elif mo_len == 2:
            return "unrestricted"
        else:
            return "none"
    
    @classmethod
    def from_parser(self, parser):
        """
        Construct a Metadata object from an output file parser.
        
        :param parser: Output data parser.
        :return: A populated Metadata object.
        """
        try:
            # Do some processing on time objects.
            if parser.data.metadata.get('wall_time', None) is not None:
                duration = sum(parser.data.metadata['wall_time'], timedelta())
                
            elif parser.data.metadata.get('cpu_time', None) is not None and parser.data.metadata.get('num_cpus', None) is not None:
                # We can estimate the wall time from the CPU time and number of cpus.
                duration = sum(parser.data.metadata['cpu_time'], timedelta()) / parser.data.metadata['num_cpus']
                
            else:
                duration = None
                
            if parser.data.metadata.get('date', None) is not None:
                date = datetime.fromtimestamp(parser.data.metadata['date'])
            else:
                date = None
                
            memory_used = Memory(parser.data.metadata['memory_used']) if "memory_used" in parser.data.metadata else None
            memory_available = Memory(parser.data.metadata['memory_available']) if "memory_available" in parser.data.metadata else None
            
            # TODO: This doesn't seem to make sense; the parser already contains a metadata object...
            return self(
                jobId = parser.data.metadata.get('jobId', None),
                name = parser.data.metadata.get('name', None),
                user = parser.data.metadata.get('user', None),
                userID = parser.data.metadata.get('userID', None),
                log_files = parser.data.metadata.get('log_files', None),
                auxiliary_files = parser.data.metadata.get('auxiliary_files', None),
                date = date,
                duration = duration,
                package = parser.data.metadata.get('package', None),
                package_version = parser.data.metadata.get('package_version', None),
                calculations = self.get_calculation_types_from_cclib(parser.data),
                success = parser.data.metadata.get('success', None),
                methods = self.get_methods_from_cclib(parser.data),
                functional = parser.data.metadata.get('functional', None),
                basis_set = parser.data.metadata.get('basis_set', None),
                
                charge = getattr(parser.data, 'charge', None),
                multiplicity = getattr(parser.data, 'mult', None),
                optimisation_converged = getattr(parser.data, 'optdone', None),
                temperature = getattr(parser.data, 'temperature', None),
                pressure = getattr(parser.data, 'pressure', None),
                orbital_spin_type = self.get_orbital_spin_type_from_cclib(parser.data),
                
                solvent = Solvent.from_parser(parser),
                
                num_cpu = parser.data.metadata.get('num_cpu', None),
                memory_available = memory_available,
                memory_used = memory_used,
                performance = Performance.from_parser(parser) if "performance" in parser.data.metadata else None
            )
        except AttributeError:
            # There is no metadata available, give up.
            raise Result_unavailable_error("Metadata", "no metadata is available")
        
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        # Most attributes we can just dump as is.
        attr_dict = {
            "name": self.name,
            "log_files": [str(log_file) for log_file in self.log_files],
            "auxiliary_files": {aux_file_name: str(aux_file) for aux_file_name, aux_file in self.auxiliary_files.items()}
        }
        
        attrs = [
            "jobId",
            "history",
            "charge",
            "multiplicity",
            "user",
            "userID",
            "package",
            "package_version",
            "digichem_version",
            "calculations",
            "methods",
            "functional",
            "basis_set",
            "orbital_spin_type",
            "success",
            "optimisation_converged",
        ]
        attr_dict.update({attr: getattr(self, attr) for attr in attrs})
        
        # Add some more complex stuff.
        attr_dict['date'] = {
            "value": self.date.timestamp() if self.date is not None else None,
            "units": "s",
            "string": date_to_string(self.date) if self.date is not None else None
        }
        attr_dict['insert_date'] = {
            "value": self.insert_date.timestamp() if self.insert_date is not None else None,
            "units": "s",
            "string": date_to_string(self.insert_date) if self.insert_date is not None else None
        }
        attr_dict['duration'] = {
            "value": self.duration.total_seconds() if self.duration is not None else None,
            "units": "s",
            "string": timedelta_to_string(self.duration) if self.duration is not None else None
        }
        attr_dict["temperature"] = {
            "value": self.temperature,
            "units": "K"
        }
        attr_dict["pressure"] = {
            "value": self.pressure,
            "units": "atm"
        }
        attr_dict["solvent"] = self.solvent.dump(digichem_options, all)
        
        attr_dict['num_cpu'] = self.num_cpu
        for attr_name in ("memory_used", "memory_available"):
            if getattr(self, attr_name) is not None:
                value, unit = getattr(self, attr_name).auto_units
                attr_dict[attr_name] = {
                    "value": value,
                    "units": unit
                }
            else:
                attr_dict[attr_name] = {
                    "value": None,
                    "units": None
                }
        
        attr_dict['performance'] = self.performance.dump(digichem_options, all) if self.performance else None

        return attr_dict
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        # Assemble our args to pass to the constructor.
        # Most of these can be used from the metadata dict as is.
        kwargs = copy.deepcopy(data)
        
        # For more complex fields, use the data item.
        for attr in ['insert_date', 'date', 'duration', 'temperature', "pressure"]:
            if attr in data:
                kwargs[attr] = data[attr]['value']
            
            else:
                kwargs[attr] = None
        
        kwargs['insert_date'] = datetime.fromtimestamp(kwargs['insert_date']) if kwargs['insert_date'] is not None else None
        kwargs['date'] = datetime.fromtimestamp(kwargs['date']) if kwargs['date'] is not None else None
        kwargs['duration'] = timedelta(seconds = kwargs['duration'])  if kwargs['duration'] is not None else None
        
        kwargs['solvent'] = Solvent.from_dump(data.get('solvent', {}), result_set, options)
        
        for attr_name in ("memory_used", "memory_available"):
            if attr_name in data and data[attr_name]['value'] is not None:
                kwargs[attr_name] = Memory.from_units(data[attr_name]["value"], data[attr_name]["units"])
            
            else:
                kwargs[attr_name] = None

        kwargs['performance'] = Performance.from_dump(data['performance'], result_set, options) if "performance" in data and data['performance'] is not None else None
        
        return self(**kwargs)


class Merged_metadata(Metadata):
    """
    A modified metadata class for merged calculation results.
    """
    
    
    def __init__(self, num_calculations, *args, **kwargs):
        """
        :param num_calculations: The number of merged calculations this metadata represents.
        """
        super().__init__(*args, log_files = None, auxiliary_files = None, **kwargs)
        self.num_calculations = num_calculations
        
    @classmethod
    def from_metadatas(self, *multiple_metadatas):
        """
        Create a merged metadata object from multiple metadata objects.
        
        :param multiple_metadatas: A list of metadata objects to merge.
        :returns: A new Merged_metadata object.
        """
        # Our merged metadata.
        merged_metadata = self(num_calculations = len(multiple_metadatas))
        for attr in ("name", "user", "userID", "package", "package_version", "functional", "basis_set"):
            setattr(merged_metadata, attr, self.merged_attr(attr, multiple_metadatas))
            
        # We take the latest of the two dates.
        merged_metadata.date = latest_datetime(*[metadata.date for metadata in multiple_metadatas if metadata.date is not None])
        # And the total duration.
        merged_metadata.duration = total_timedelta(*[metadata.duration for metadata in multiple_metadatas if metadata.duration is not None])
        
        # Merge methods and calculations (but keep unique only).
        merged_metadata.calculations = list(dict.fromkeys(itertools.chain(*(metadata.calculations for metadata in multiple_metadatas))))
        merged_metadata.methods = self.sorted_methods(set(itertools.chain(*(metadata.methods for metadata in multiple_metadatas))))
        
        # Keep the solvent if it's the same for all, otherwise discard.
        merged_metadata.solvent = multiple_metadatas[0].solvent.merge(*[other.solvent for other in multiple_metadatas[1:]])
        
        # Combine performance data.
        merged_metadata.num_cpu = sum([meta.num_cpu for meta in multiple_metadatas if meta.num_cpu is not None])
        merged_metadata.memory_available = Memory(sum([int(meta.memory_available) for meta in multiple_metadatas if meta.memory_available is not None]))
        merged_metadata.memory_used = Memory(sum([int(meta.memory_used) for meta in multiple_metadatas if meta.memory_used is not None]))
        
        # We are only successful if all calcs are successful.
        merged_metadata.success = all((metadata.success for metadata in multiple_metadatas))
        converged = [metadata.optimisation_converged for metadata in multiple_metadatas]
        merged_metadata.optimisation_converged = False if False in converged else True if True in converged else None
        
        # Only keep these if all the same.
        for attr in ("temperature", "pressure"):
            # Get a unique list of the attrs.
            attr_set = set(getattr(metadata, attr) for metadata in multiple_metadatas)
            
            # Only keep if we have exactly one entry.
            if len(attr_set) == 1:
                setattr(merged_metadata, attr, tuple(attr_set)[0])
                
        # We take the first orbital spin type charge and mult, as this is what our orbital list will actually be.
        merged_metadata.orbital_spin_type = multiple_metadatas[0].orbital_spin_type
        merged_metadata.charge = multiple_metadatas[0].charge
        merged_metadata.multiplicity = multiple_metadatas[0].multiplicity
        
        # CAREFULLY merge aux files, so that later files do not overwrite earlier ones.
        # This is useful behaviour because it matches how other results are merged, so certain aux files
        # (turbomole density files) will still match their respective results.
        for metadata in multiple_metadatas:
            merged_metadata.auxiliary_files = conservative_merger.merge(merged_metadata.auxiliary_files, metadata.auxiliary_files)
        
        return merged_metadata

class Performance(Result_object):
    """
    Performance metrics and profiling data for the calculation.
    """

    def __init__(
        self,
        duration = [],
        memory_used = [],
        memory_used_percent = [],
        memory_available = [],
        memory_available_percent = [],
        cpu_used = [],
        output_available = [],
        scratch_used = [],
        scratch_available = [],
        memory_allocated = None,
        cpu_allocated = None,
    ):
        self.duration = duration
        self.memory_used = memory_used
        self.memory_used_percent = memory_used_percent
        self.memory_available = memory_available
        self.memory_available_percent = memory_available_percent
        self.cpu_used = cpu_used
        self.output_available = output_available
        self.scratch_used = scratch_used
        self.scratch_available = scratch_available

        self.memory_allocated = memory_allocated if memory_allocated is not None else self.max_mem
        self.cpu_allocated = cpu_allocated if cpu_allocated is not None else math.ceil(max(cpu_used) / 100)

    @property
    def output_space(self):
        warnings.warn("output_space is deprecated, use output_available instead", DeprecationWarning)
        return self.output_available

    @property
    def scratch_space(self):
        warnings.warn("scratch_space is deprecated, use scratch_available instead", DeprecationWarning)
        return self.scratch_available
    

    @classmethod
    def from_parser(self, parser):
        """
        Construct a Performance object from an output file parser.
        
        :param parser: Output data parser.
        :return: A populated Performance object.
        """
        return self(
            duration = parser.data.metadata['performance']['duration'].tolist(),
            memory_used = parser.data.metadata['performance']['memory_used'].tolist(),
            memory_allocated = Memory(parser.data.metadata['memory_available']) if "memory_available" in parser.data.metadata else None,
            memory_used_percent = parser.data.metadata['performance']['memory_used_percent'].tolist(),
            memory_available = parser.data.metadata['performance']['memory_available'].tolist(),
            memory_available_percent = parser.data.metadata['performance']['memory_available_percent'].tolist(),
            cpu_used = parser.data.metadata['performance']['cpu_used'].tolist(),
            cpu_allocated = parser.data.metadata.get('num_cpu', None),
            output_available = parser.data.metadata['performance']['output_available'].tolist(),
            scratch_used = parser.data.metadata['performance']['scratch_used'].tolist() if 'scratch_used' in parser.data.metadata['performance'] else [0] * len(parser.data.metadata['performance']['duration']),
            scratch_available = parser.data.metadata['performance']['scratch_available'].tolist()
        )
    
    @property
    def max_mem(self):
        """
        The maximum amount of memory used in the calculation (in bytes)
        """
        return max(self.memory_used)
    
    @property
    def memory_margin(self):
        max_memory = float(self.memory_allocated) if self.memory_allocated is not None else self.max_mem
        
        return max_memory - self.max_mem
    
    @property
    def memory_efficiency(self):
        """
        Calculate the memory efficiency of this calculation.

        :param max_memory: The amount of allocated memory (in bytes), this will be guestimated automatically if not available.
        """
        # Integrate to find the number of byte seconds used.
        area = integrate.trapezoid(self.memory_used, self.duration)

        # How much we should/could have used.
        total_area = (self.duration[-1] - self.duration[0]) * float(self.memory_allocated)

        # Return as %.
        try:
            return area / total_area * 100
        
        except Exception:
            return 0
    
    @property
    def cpu_efficiency(self):
        """
        Calculate the CPU efficiency of this calculation.

        :param max_cpu: The number of allocated CPUs, this will be guestimated automatically if not available.
        """
        # Integrate to find the number of CPU seconds used.
        area = integrate.trapezoid(self.cpu_used, self.duration)

        # How much we should/could have used.
        total_area = (self.duration[-1] - self.duration[0]) * self.cpu_allocated * 100

        # Return as %.
        try:
            return area / total_area * 100

        except Exception:
            # Div zero
            return 0
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        duration = [0.0] * len(data['values'])
        memory_used = [0.0] * len(data['values'])
        memory_allocated = Memory(data['memory_allocated']['value'])
        memory_used_percent = [0.0] * len(data['values'])
        memory_available = [0.0] * len(data['values'])
        memory_available_percent = [0.0] * len(data['values'])
        cpu_used = [0.0] * len(data['values'])
        cpu_allocated = data['cpu_allocated']
        output_available = [0.0] * len(data['values'])
        scratch_used = [0.0] * len(data['values'])
        scratch_available = [0.0] * len(data['values'])

        for i, value in enumerate(data['values']):
            duration[i] = value['duration']['value']
            memory_used[i] = value['memory_used']['value']
            memory_used_percent[i] = value['memory_used_percent']['value']
            memory_available[i] = value['memory_available']['value']
            memory_available_percent[i] = value['memory_available_percent']['value']
            cpu_used[i] = value['cpu_used']['value']
            output_available[i] = value['output_space']['value']
            if 'scratch_used' in value:
                scratch_used[i] = value['scratch_used']['value']
            scratch_available[i] = value['scratch_space']['value']
        
        return self(
            duration = duration,
            memory_used = memory_used,
            memory_allocated = memory_allocated,
            memory_used_percent = memory_used_percent,
            memory_available = memory_available,
            memory_available_percent = memory_available_percent,
            cpu_used = cpu_used,
            cpu_allocated = cpu_allocated,
            output_available = output_available,
            scratch_used = scratch_used,
            scratch_available = scratch_available
        )

    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "cpu_allocated": self.cpu_allocated,
            "cpu_efficiency": {
                "units": "%",
                "value": float(self.cpu_efficiency),
            },
            "memory_allocated": {
                "units": "bytes",
                "value": float(self.memory_allocated)
            },
            "maximum_memory": {
                "units": "bytes",
                "value": self.max_mem,
            },
            "memory_margin": {
                "units": "bytes",
                "value": self.memory_margin
            },
            "memory_efficiency": {
                "units": "%",
                "value": float(self.memory_efficiency)
            },
            "values":[
                {
                    'duration': {
                        "units": "s",
                        "value": self.duration[i]
                    },
                    'memory_used': {
                        "units": "bytes",
                        "value": self.memory_used[i]
                    },
                    'memory_used_percent': {
                        "units": "%",
                        "value": self.memory_used_percent[i]
                    },
                    'memory_available': {
                        "units": "bytes",
                        "value": self.memory_available[i]
                    },
                    'memory_available_percent': {
                        "units": "bytes",
                        "value": self.memory_available_percent[i]
                    },
                    'cpu_used': {
                        "units": "%",
                        "value": self.cpu_used[i]
                    },
                    'output_space': {
                        "units": "bytes",
                        "value": self.output_space[i]
                    },
                    'scratch_used': {
                        "units": "bytes",
                        "value": self.scratch_used[i]
                    },
                    'scratch_space': {
                        "units": "bytes",
                        "value": self.scratch_space[i]
                    }
                } for i in range(len(self.duration))
            ]
        }

        return {
            'duration': {
                "units": "s",
                "values": self.duration.tolist()
            },
            'memory_used': {
                "units": "bytes",
                "values": self.memory_used.tolist()
            },
            'memory_used_percent': {
                "units": "%",
                "values": self.memory_used_percent.tolist()
            },
            'memory_available': {
                "units": "bytes",
                "values": self.memory_available.tolist()
            },
            'memory_available_percent': {
                "units": "bytes",
                "values": self.memory_available_percent.tolist()
            },
            'cpu_used': {
                "units": "%",
                "values": self.cpu_used.tolist()
            },
            'output_space': {
                "units": "bytes",
                "values": self.output_space.tolist()
            },
            'scratch_space': {
                "units": "bytes",
                "values": self.scratch_space.tolist()
            }
        }

