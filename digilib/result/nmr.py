import numpy
from itertools import filterfalse
import periodictable
from fractions import Fraction
import re
import statistics
import math
from configurables.misc import is_int

from digilib.misc.base import regular_range, powerset
from digilib.exception.base import Result_unavailable_error
from digilib.result.base import Result_object, Result_container, Floatable_mixin
import digilib.log

# Hidden.
#from digilib.result.spectroscopy import Combined_graph


# TODO: NMR tensors are currently not re-orientated according to the alignment method used.
# This needs implementing.
class NMR_spectrometer(Result_object):
    """
    A class for generating NMR spectra on-demand.
    """
    
    def __init__(self, nmr_results, frequency = 300, fwhm = 0.001, resolution = 0.001, cutoff = 0.01, y_filter = 1e-6, coupling_filter = 0.1, pre_merge = 0.01, post_merge = None, isotope_options = None):
        """
        Constructor for NMR_spectrometer.
        
        :param nmr_results: A list of NMR_group result objects.
        :param frequency: The frequency of the simulated spectrometer.
        :param fwhm: The full-width at half-maximum of the simulated peaks (in ppm).
        :param resolution: The resolution/step-size of the plotted peaks (in ppm). Decreasing this value may increase computational time.
        :param reference: An optional reference isotropic value to correct this shielding by.
        """
        self.nmr_results = nmr_results
        self.frequency = frequency
        self.pre_merge = pre_merge
        self.post_merge = post_merge
        self.fwhm = fwhm
        self.gaussian_resolution = resolution
        self.gaussian_cutoff = cutoff
        self.y_filter = y_filter
        self.coupling_filter = coupling_filter
        self._isotope_options = isotope_options if isotope_options is not None else {}
        
    def isotope_options(self, element, isotope):
        options = {
            "frequency": self.frequency,
            "pre_merge": self.pre_merge,
            "post_merge": self.post_merge,
            "fwhm": self.fwhm,
            "y_filter": self.y_filter,
            "gaussian_resolution": self.gaussian_resolution,
            "coupling_filter": self.coupling_filter,
            "gaussian_cutoff": self.gaussian_cutoff,
        }
        options.update(self._isotope_options.get("{}{}".format(isotope, periodictable.elements[element].symbol), {}))
        return options
        
    @classmethod
    def from_options(self, nmr_results, *, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            nmr_results,
            frequency = options['nmr']['frequency'],
            fwhm = options['nmr']['fwhm'],
            resolution = options['nmr']['gaussian_resolution'],
            cutoff = options['nmr']['gaussian_cutoff'],
            y_filter = options['nmr']['y_filter'],
            coupling_filter = options['nmr']['coupling_filter'],
            pre_merge = options['nmr']['pre_merge'],
            post_merge = options['nmr']['post_merge'],
            isotope_options = options['nmr']['isotopes'],
            **kwargs
        )
    
    @property
    def available(self):
        """
        Determine which nuclei NMR can be simulated for.
        
        :return: A set of the available nuclei, each as a tuple of (proton_number, isotope). Note that the special value of 0 is used to indicate isotope-independent data is available. Eg, (6, 0) would indicate isotope-independent (no coupling) carbon NMR. 
        """
        nuclei = set()
        
        for nmr_result in self.nmr_results.values():
            # Add the generic (non-isotope specific) element to our list of supported types.
            nuclei.add((nmr_result.group.element.number, 0, ()))
            
            couplings = {}
            
            # Also consider specific isotope coupling.
            for coupling in nmr_result.couplings:
                main_index = coupling.groups.index(nmr_result.group)
                second_index = abs(1 - main_index)
                
                if coupling.isotopes[main_index] not in couplings:
                    couplings[coupling.isotopes[main_index]] = set()
                
                couplings[coupling.isotopes[main_index]].add( (coupling.groups[second_index].element.number, coupling.isotopes[second_index]) )
            
            for isotope, isotope_couplings in couplings.items():
                for combination in powerset(sorted(isotope_couplings)):
                    nuclei.add((nmr_result.group.element.number, isotope, tuple(combination)))
                
        return {
            self.make_shortcode(element, isotope, decoupling):
            {
                "element": element,
                "isotope": isotope,
                #"decoupling": [list(decouple) for decouple in decoupling],
                "decoupling": list(decoupling)
            } for element, isotope, decoupling in sorted(nuclei, key = lambda nmr: (nmr[0], nmr[1]))
        }
    
    @classmethod
    def parse_shortcode(self, code):
        """
        Parse an NMR shortcode.
        
        Each shortcode specifies an element and one of its isotopes to investigate, optionally followed by a list of isotopes to decouple (not show coupling for).
        For example:
            13C{1H}    - Carbon-13 spectrum with protons decoupled.
            1H{1H,13C} - Hydrogen-1 spectrum with protons and 13C decoupled.
            *C or 0C   - Carbon NMR spectrum with no coupling considered.
        """
        match = re.search(r'(\d+|\*)([A-z][A-z]?){?((\d+[A-z][A-z]?,?)*)}?', code)
        
        if match is None:
            raise ValueError("Unable to process NMR shortcode '{}'".format(code))
        
        match_groups = match.groups()
        isotope = int(match_groups[0]) if match_groups[0] != "*" else 0
        element = getattr(periodictable, match_groups[1]).number
        
        decoupling = [tuple(re.search(r"(\d+)([A-z][A-z]?)", decouple).groups()) for decouple in match_groups[2].split(",") if re.search(r"(\d+)([A-z][A-z]?)", decouple) is not None]
        decoupling = sorted([(getattr(periodictable, ele).number, int(iso)) for iso, ele in decoupling])
        
        return element, isotope, decoupling
    
    @classmethod
    def make_shortcode(self, element, isotope, decoupling):
        """
        Create an NMR shortcode from a given NMR experiment.
        """
        return "{}{}".format(isotope, periodictable.elements[element].symbol) + (
                "{{{}}}".format(",".join(["{}{}".format(decouple[1], periodictable.elements[decouple[0]]) for decouple in decoupling]))
                if len(decoupling) > 0 else "")
        
    
    def __contains__(self, item):
        """
        Can this spectrometer generate a given spectrum.
        """
        element, isotope, decoupling = self.parse_shortcode(item)
        
        #for available_element, available_isotope, available_decoupling in self.available.values():
        for available in self.available.values():
            if element == available['element'] and isotope == available['isotope'] and sorted(decoupling) == sorted(available['decoupling']):
                return True
        
        return False
    
    def __iter__(self):
        return iter(self.available.keys())
    
    def __getitem__(self, item):
        """
        Generate a spectrum for a given shortcode.
        """
        # We return a lambda function here because of the slightly strange dumping mechanism.
        # This object is setup as a generate_for_dump() target. generate_for_dump() returns a
        # dict of callables. We don't actually require this behaviour, but we still need to 
        # match the interface.
        return lambda digichem_options: self.spectrum(*self.parse_shortcode(item))
    
    def _get_dump_(self):
        """
        Method used to get a dictionary used to generate on-demand values for dumping.
        
        This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
        
        Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
        """
        return self
    
    def _dump_(self, digichem_options, all):
        # Return a list of possible spectra we can generate.
        return {
            "codes": list(self.available.keys()),
            #"available": list(self.available.values()),
        }
        
    def get_graph(self, element, isotope, decoupling = None):
        """
        Return a result.spectroscopy.Combined_graph object that can be used to generate a spectrum (with Gaussian broadened peaks) for a given experiment.
        
        If you are only interested in the spectrum itself (not the machinery that generates it), try spectrum().
        """
        from digilib.result.spectroscopy import Combined_graph
        
        
        decoupling = decoupling if decoupling is not None else []
         
        # First, simulate vertical peaks.
        # Each peaks is grouped by the atom-group that causes it.
        # This allows us to plot merged peaks separately, for example if two atoms result in overlapping chemical shift.
        grouped_peaks = self.simulate(element, isotope, decoupling)
         
        # Get options specific to the isotope we're looking at.
        isotope_options = self.isotope_options(element, isotope)
        
        # The total spectrum takes all simulated peaks.
        # These are grouped by atom_group, flatten this list before passing to spectroscopy.
        graph = Combined_graph.from_nmr(grouped_peaks, isotope_options['fwhm'], isotope_options['gaussian_resolution'], isotope_options['gaussian_cutoff'], filter = isotope_options['y_filter'])
        
        return graph
    
    def spectrum(self, element, isotope, decoupling = None):
        """
        Simulate a full NMR spectrum.
        
        :param element: The element to simulate (as a proton number).
        :param isotope: The isotope of the element to simulate.
        :param decoupling: A list of elements to 'decouple' (not consider couplings to). Each element should be specified as a tuple of (proton_number, isotope).
        """
        # Get options specific to the isotope we're looking at.
        isotope_options = self.isotope_options(element, isotope)
        
        graph = self.get_graph(element, isotope, decoupling)
         
        digilib.log.get_logger().info("Simulating broadened NMR peaks for {}{} at {} MHz with {} decoupling".format(
            isotope,
            periodictable.elements[element].symbol,
            isotope_options['frequency'],
            ",".join(["{}{}".format(decouple_iso, periodictable.elements[decouple_ele].symbol) for decouple_ele, decouple_iso in decoupling] if len(decoupling) != 0 else ["no"])
        ))
         
        values = {
            "values": [
                {"x":{"value": float(x), "units": "ppm"}, "y": {"value": float(y), "units": "arb"}}
                for x,y in graph.plot_cumulative_gaussian()
            ],
            "groups": {
                atom_group.label: [
                    {"x":{"value": float(x), "units": "ppm"}, "y": {"value": float(y), "units": "arb"}}
                    for x,y in group_spectrum.plot_cumulative_gaussian()
                ]
                for atom_group, group_spectrum in graph.graphs.items()
            }
        }
        digilib.log.get_logger().info("Done simulating NMR spectrum")
        return values
    
    def coupling(self, element, isotope, decoupling = None):
        """
        Return couplings for a given experiment.
        
        Be aware that small couplings may be filtered out (not returned) depending on the value of coupling_filter.
        
        :returns: The relevant couplings, as a dictionary of dictionary of dicts. The keys of the two outer dict correspond to the atom groups this coupling is between. The inner dict is a coupling dict, containing 'total', 'isotopes' and 'groups'.
        """
        isotope_options = self.isotope_options(element, isotope)
        outer_coupling = {}
        for nmr_result in self.nmr_results.values():
            if nmr_result.group.element.number != element:
                # Wrong element, skip.
                continue
                
            # Matches our element.
            inner_coupling = {}
            for coupling in nmr_result.couplings:
                # Check the main isotope matches.
                main_index = coupling.groups.index(nmr_result.group)
                second_index = 1 if main_index == 0 else 0
                
                # Only include this coupling if it involves our isotope,
                # and if it's value is above our filter.
                if coupling.isotopes[main_index] == isotope and abs(coupling.total) > isotope_options['coupling_filter']:
                    no_couple = False
                    # Check we haven't been asked to de-couple this coupling.
                    for decouple_element, decouple_isotope in decoupling:
                        if coupling.groups[second_index].element.number == decouple_element \
                        and coupling.isotopes[second_index] == decouple_isotope:
                            # This coupling is good.
                            no_couple = True
                            break
                    
                    if not no_couple:
                        #inner_coupling[(coupling.groups[second_index], coupling.isotopes[second_index])] = coupling
                        
                        # Create a dictionary for this foreign atom.
                        if coupling.groups[second_index] not in inner_coupling:
                            inner_coupling[coupling.groups[second_index]] = {}
                        
                        # Add another for the isotope of this foreign atom.
                        inner_coupling[coupling.groups[second_index]][coupling.isotopes[second_index]] = coupling
            
            outer_coupling[nmr_result.group] = inner_coupling
        
        return outer_coupling
    
    def split_peaks(self, nmr_result, coupling, peaks, isotope_options):
        """
        Split a list of peaks according to a given coupling.
        
        :param nmr_result: The NMR result which gives rise to the peak(s) being split.
        :param coupling: The coupling that is being used to cause the splitting.
        :param peaks: An existing list of peaks to split.
        :param isotope_options: Isotope options for this atom.
        """
        main_index = coupling.groups.index(nmr_result.group)
        second_index = 1 if main_index == 0 else 0
        
        # Each atom in the group we are coupling to.
        for atom in range(coupling.num_coupled_atoms(nmr_result.group)):
            new_peaks = {}
            for old_peak in peaks.values():
                # Calculate the shift of the new peaks (in ppm).
                coupling_constant = coupling.total / isotope_options['frequency']
                
                # If we're merging peaks, round the value appropriately.
                if isotope_options['pre_merge']:
                    coupling_constant = round(coupling_constant / isotope_options['pre_merge']) * isotope_options['pre_merge']
                
                # Bafflingly, calling 'neutron' here is necessary to make nuclear_spin available.
                ele = getattr(periodictable, coupling.groups[second_index].element.symbol)
                iso = ele[coupling.isotopes[second_index]]
                iso.neutron
                spin = float(Fraction(iso.nuclear_spin))
                #abundance = iso.abundance /100
                num_peaks = 2 * spin +1
                
                # Add the new peaks to the shifts originating from coupling between these two groups.
                for new_peak in ([sub_peak, old_peak[1] / num_peaks] for sub_peak in regular_range(old_peak[0], num_peaks, coupling_constant)):
                    if new_peak[0] in new_peaks:
                        new_peaks[new_peak[0]][1] += new_peak[1]
                    
                    else:
                        new_peaks[new_peak[0]] = new_peak
            
            peaks = new_peaks
            
        return peaks

    def simulate(self, element, isotope, decoupling = None):
        """
        Simulate vertical NMR peaks.
        
        :param element: The element to simulate (as a proton number).
        :param isotope: The isotope of the element to simulate.
        :param decoupling: A list of elements to 'decouple' (not consider couplings to). Each element should be specified as a tuple of (proton_number, isotope).
        """
        element = int(element)
        isotope = int(isotope)
        decoupling = [(int(decouple_ele), int(decouple_iso)) for decouple_ele, decouple_iso in decoupling] if decoupling is not None else None
        isotope_options = self.isotope_options(element, isotope)
        
        group_peaks = {}
        
        digilib.log.get_logger().info("Simulating NMR peaks for {}{} at {} MHz with {} decoupling".format(
            isotope,
            periodictable.elements[element].symbol,
            isotope_options['frequency'],
            ",".join(["{}{}".format(decouple_iso, periodictable.elements[decouple_ele].symbol) for decouple_ele, decouple_iso in decoupling] if len(decoupling) != 0 else ["no"])
        ))
        
        # Get relevant couplings.
        all_coupling = self.coupling(element, isotope, decoupling)
        
        for nmr_result in self.nmr_results.values():
            if nmr_result.group.element.number != element:
                # Wrong element, skip.
                continue
                
            # Matches our element.
            
            group_coupling = all_coupling.get(nmr_result.group, {})
            
            # Make some peaks.
            # Start with a single shift peak.
            # 0: chemical shift, 1: intensity
            #
            # If we are merging to a nearest point, do so now (to preserve symmetry).
            if isotope_options['pre_merge']:
                # FYI: Round breaks ties by rounding to the nearest even number (banker's rounding), not by rounding upwards.
                initial_shielding = round(nmr_result.shielding / isotope_options['pre_merge']) * isotope_options['pre_merge']
                
            else:
                initial_shielding = nmr_result.shielding
            
            peaks = {initial_shielding: [initial_shielding, len(nmr_result.group.atoms)]}
            
            # Now split it by each coupling.
            for foreign_atom_group, isotopes in group_coupling.items():
                # Each atom-group that we're going to couple to may have couplings for a number of isotopes.
                # Additionally, any percentage abundance not taken up by the isotopes for which coupling is available
                # will result in an unsplit peak.
                total_abundance = sum([foreign_atom_group.element[isotope].abundance for isotope in isotopes.keys()])
                
                new_peaks = {}
                
                for isotope, coupling in isotopes.items():
                    # Decrease the apparent peak intensity by the natural abundance.
                    
                    abundance = foreign_atom_group.element[isotope].abundance / 100
                    isotope_peaks = {peak[0]: (peak[0], peak[1] * abundance) for peak in peaks.values()}
                    
                    for new_peak in self.split_peaks(nmr_result, coupling, isotope_peaks, isotope_options).values():
                        if new_peak[0] in new_peaks:
                            new_peaks[new_peak[0]] = (new_peak[0], new_peaks[new_peak[0]][1] + new_peak[1])
                        else:
                            new_peaks[new_peak[0]] = new_peak
                            
                # After splitting by all isotopes, if we have any abundance not accounted for,
                # we'll add back the original unsplit peaks (with their intensity decreased by the remaining
                # abundance).
                if 100 - total_abundance > 0.0:
                    # The starting height of any remaining unsplit peaks.
                    abundance = (100 - total_abundance) / 100
                    
                    for new_peak in ((peak[0], peak[1] * abundance) for peak in peaks.values()):
                        if new_peak[0] in new_peaks:
                            new_peaks[new_peak[0]] = (new_peak[0], new_peaks[new_peak[0]][1] + new_peak[1])
                        
                        else:
                            new_peaks[new_peak[0]] = new_peak
                            
                peaks = new_peaks
            
            peaks = list(peaks.values())
            peaks.sort(key = lambda peak: peak[0])
            
            # If we've been asked to, merge similar peaks.
            # We do this last so as to not carry forward rounding and averaging errors.
            if isotope_options['post_merge']:
                # Each generated peak will be aligned to a 'grid', the spacing of which is
                # given by post_merge.
                # Start by generating our grid.
                # We want to make sure there is a grid point exactly on our mid-point, this should preserve symmetry.
                median = statistics.median((peak[0] for peak in peaks))
                start = median - math.ceil((median - peaks[0][0]) / isotope_options['post_merge']) * isotope_options['post_merge']
                stop =  median + math.ceil((peaks[-1][0] - median) / isotope_options['post_merge']) * isotope_options['post_merge']
                steps = round(((stop - start) / isotope_options['post_merge'])) +1
                new_shifts = numpy.linspace(start, stop, steps)
                
                new_peaks = {}
                new_shift_index = 0
                
                for peak in peaks:
                    # Find where this peak best aligns to our grid.
                    # Because we are going through in order, we don't need to start from the beginning of our new peaks.
                    for index, new_shift in enumerate(new_shifts[new_shift_index:]):
                        if peak[0] < (new_shift + isotope_options['post_merge'] /2):
                            if new_shift not in new_peaks:
                                # New peaks at this shift.
                                new_peaks[new_shift] = [new_shift,  peak[1]]
                            
                            else:
                                # Existing peak here, add to intensity.
                                new_peaks[new_shift][1] += peak[1]
                            
                            # Update our index.
                            new_shift_index = index
                            
                            break
                        
                peaks = list(new_peaks.values())
                
            # Done for this group of atoms.
            group_peaks[nmr_result.group] = peaks
        
        digilib.log.get_logger().info("Done simulating NMR peaks")
        
        return group_peaks
    

class NMR_list(Result_container):
    """
    A container to hold a list of NMR results.
    
    For practical applications, see the output of the group() method or the spectrometer attribute.
    """
    
    def __init__(self, *args, atoms, options, **kwargs):
        """
        :param *args: A list of NMR objects.
        :param atoms: An Atom_list object.
        :param options: A digichem options object.
        """
        super().__init__(*args, **kwargs)
        self.atoms = atoms
        self.options = options
        self.groups = self.group()
        self.spectrometer = NMR_spectrometer.from_options(self.groups, options = options)
        # This is used as part of the dump mechanic.
        self.spectrum = self.spectrometer
    
    @classmethod
    def from_parser(self, parser):
        return self(NMR.list_from_parser(parser), atoms = parser.results.atoms, options = parser.options)
    
    def find(self, criteria = None, *, label = None, index = None):
        return self.search(criteria = criteria, label = label, index = index, allow_empty = False)[0]
    
    def search(self, criteria = None, *, label = None, index = None, allow_empty = True):
        """
        """
        if label is None and index is None and criteria is None:
            raise ValueError("One of 'criteria', 'label' or 'index' must be given")
        
        if criteria is not None:
            if str(criteria).isdigit():
                index = int(criteria)
            
            else:
                label = criteria
        
        # Now get our filter func.
        if label is not None:
            filter_func = lambda nmr: nmr.atom.label != label
        
        elif index is not None:
            filter_func = lambda nmr: nmr.atom.index != index
        
        # Now search.
        found = type(self)(filterfalse(filter_func, self), atoms = self.atoms, options = self.options)
        
        if not allow_empty and len(found) == 0:
            if label is not None:
                criteria_string = "label = '{}'".format(label)
            
            elif index is not None:
                criteria_string = "index = '{}'".format(index)
            
            raise Result_unavailable_error("NMR", "could not find NMR data for atom '{}'".format(criteria_string))
        
        return found
        
    def group(self, no_self_coupling = True):
        """
        """
        # First, decide which atoms are actually equivalent.
        # We can do this by comparing canonical SMILES groupings.
        atom_groups = self.atoms.groups

        # Also get the bond distance matrix for the molecule, so we can work out which couplings are actually equivalent.
        bond_matrix = self.atoms.bond_matrix
        
        nmr_groups = {}
        # Next, assemble group objects.
        for group_id, atom_group in atom_groups.items():
            nmr_results = [nmr_result for nmr_result in self if nmr_result.atom in atom_group.atoms]
            
            if len(nmr_results) == 0:
                continue
            
            # Shieldings.
            shieldings = [nmr_result.shielding for nmr_result in nmr_results]
            # Only keep couplings in which at least one of the two atoms is not in this group (discard self coupling)
            couplings = [coupling for nmr_result in nmr_results for coupling in nmr_result.couplings if not no_self_coupling or len(set(atom_group.atoms).intersection(coupling.atoms)) != 2]
            
            nmr_groups[group_id] = {"group": atom_group, "shieldings": shieldings, "couplings": couplings}
        
        # Now everything is assembled into groups, re-calculate couplings based on groups only.
        # We need to do this after initial group assembly in order to discard self coupling.
        # Get unique couplings (so we don't consider any twice).
        group_couplings = {}
        unique_couplings = list({(coupling.atoms, coupling.isotopes): coupling for group in nmr_groups.values() for coupling in group['couplings']}.values())
        for coupling in unique_couplings:
            # Find the group numbers that correspond to the two atoms in the coupling.
            coupling_groups = tuple(
                [atom_group.id for atom_group in atom_groups.values() if atom in atom_group.atoms][0] for atom in coupling.atoms
            )

            # We need to ensure that coupling_groups is a unique representation of the coupling,
            # the ordering should be fixed.
            indices = [
                int(coupling_groups[0] >= coupling_groups[1]),
                int(coupling_groups[0] < coupling_groups[1])
            ]
            coupling_groups = (
                coupling_groups[indices[0]],
                coupling_groups[indices[1]]
            )
            isotopes = (
                coupling.isotopes[indices[0]],
                coupling.isotopes[indices[1]]
            )

            # The group key contains the two atom groups, and the distance between them.
            coupling_groups = (
                coupling_groups[0],
                coupling_groups[1],
                float(bond_matrix[coupling.atoms[0].index -1][coupling.atoms[1].index -1])
            )
            
            # Append the isotropic coupling constant to the group.
            if coupling_groups not in group_couplings:
                group_couplings[coupling_groups] = {}
                
            if isotopes not in group_couplings[coupling_groups]:
                group_couplings[coupling_groups][isotopes] = []
                
            group_couplings[coupling_groups][isotopes].append(coupling)
            
        # Average each 'equivalent' coupling.
        average_couplings = {
            group_key: {
                isotope_key: NMR_group_spin_coupling(
                    # TODO: Add in bond distance.
                    groups = [atom_groups[group_sub_key] for group_sub_key in group_key[:2]],
                    isotopes = isotope_key,
                    couplings = isotope_couplings,
                    distance = group_key[2]
                ) for isotope_key, isotope_couplings in  isotopes.items()}
            for group_key, isotopes in group_couplings.items()
        }
        
        # Assemble the final group objects.
        nmr_object_groups = {}
        for group_id, raw_group in nmr_groups.items():
            # Get appropriate couplings.
            
            coupling = [
                isotope_coupling
                for group_key, group_coupling in average_couplings.items()
                    for isotope_coupling in group_coupling.values()
                        if group_id in group_key[:2]
            ]
            nmr_object_groups[raw_group['group']] = NMR_group(raw_group['group'], raw_group['shieldings'], coupling)
        
        return nmr_object_groups
    
    def _dump_(self, digichem_options, all):
        """
        Dump this list of NMR results to a list of primitive types.
        """
        grouping = self.groups
        dump_dict = {
            "values": super()._dump_(digichem_options, all),
            "groups": {group_id.label: group.dump(digichem_options, all) for group_id, group in grouping.items()},
            "spectrum": self.spectrometer.dump(digichem_options, all)
        }
        return dump_dict
    
    @classmethod
    def from_dump(self, data, result_set, options):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return self(NMR.list_from_dump(data['values'], result_set, options), atoms = result_set.atoms, options = options)
    
    # def _get_dump_(self):
    #     """
    #     Method used to get a dictionary used to generate on-demand values for dumping.
         
    #     This functionality is useful for hiding expense properties from the normal dump process, while still exposing them when specifically requested.
         
    #     Each key in the returned dict is the name of a dumpable item, each value is a function to call with digichem_options as its only param.
    #     """
    #     return {
    #         "spectrum": lambda digichem_options: self.spectrometer
    #     }


class NMR_group(Result_object, Floatable_mixin):
    """
    A result object containing all the NMR related data for a group of equivalent nuclei.
    """
    
    def __init__(self, group, shieldings, couplings):
        self.group = group
        self.shieldings = shieldings
        self.couplings = sorted(
            couplings,
            key = lambda coupling: abs(coupling.total)
        )
        
        # Calculate average shieldings and couplings.
        self.shielding = float(sum([shielding.isotropic("total") for shielding in shieldings]) / len(shieldings))
        
    def __float__(self):
        return float(self.shielding)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "group": self.group.label,
            "atoms": [atom.label for atom in self.group.atoms],
            "shielding": {
                "units": "ppm",
                "value": self.shielding
            },
            #"couplings": [{"groups": [group.label for group in coupling['groups']], "isotopes": list(coupling["isotopes"]), "total": coupling["total"]} for coupling in self.couplings],
            "couplings": [coupling.dump(digichem_options, all) for coupling in self.couplings]
        }
        
class NMR_group_spin_coupling(Result_object):
    """
    A result object containing the average coupling between two different groups of nuclei.
    """
    
    def __init__(self, groups, isotopes, couplings, distance = None):
        """
        :param groups: The two atom groups that this coupling is between.
        :param isotopes: The isotopes of the two groups (the order should match that of groups).
        :param couplings: A list of individual coupling constants between the atoms of these two groups.
        :param distance: The bond distance between the two atoms.
        """
        self.groups = groups
        self.isotopes = isotopes
        self.couplings = couplings
        self.distance = int(distance) if distance is not None and is_int(distance) else distance
        
    @property
    def total(self):
        """
        The total (average coupling) between the atoms of the groups contributing to this coupling.
        """
        return sum([coupling.isotropic('total') for coupling in self.couplings]) / len(self.couplings)
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "groups": [group.label for group in self.groups],
            "isotopes": list(self.isotopes),
            "total": {
                "units": "Hz",
                "value": float(self.total),
            },
            "distance": {
                "units": "bonds",
                "value": self.distance
            }
            #"couplings": [coupling.dump(digichem_options, all) for coupling in self.couplings]
        }
        
    def other(self, atom_group):
        """
        Get the index of the other atom_group involved in this coupling.
        
        :param atom_group: One of the two atom groups.
        """
        return abs(1 - self.groups.index(atom_group))
    
    def num_coupled_atoms(self, atom_group):
        """
        Calculate the number of atoms one of the atom groups is coupled to.
        """
        # second_index = self.other(atom_group)

        # if [group.label for group in self.groups] == ["H9", "H10"] or [group.label for group in self.groups] == ["H10", "H9"]:
        #     print()

        # atoms = []
        # for coupling in self.couplings:
        #     for atom in coupling.atoms:
        #         if atom not in atom_group.atoms:
        #             atoms.append(atom)
        
        # return len(list(set(atoms)))

        return int(len(self.couplings) / len(atom_group.atoms))
        
    def multiplicity(self, atom_group):
        """
        Calculate the multiplicity (number of peaks generated) by this coupling.
        
        :param atom_group: The atom_group who's corresponding peak is to be split. This should be one of the two groups in self.groups.
        """
        second_index = self.other(atom_group)
        # Calculate how many peaks are going to be generated.
        # This is the number of equivalent nuclei * (2 * spin) + 1
        
        # First, determine the number of atoms beings coupled to.
        # Note that this is not simply the number of atoms in the other atom group, because
        # not all atoms of the main group are necessarily coupled to all atoms of the foreign group
        num_coupled_atoms = self.num_coupled_atoms(atom_group)
        
        ele = getattr(periodictable, self.groups[second_index].element.symbol)
        iso = ele[self.isotopes[second_index]]
        iso.neutron
        spin = float(Fraction(iso.nuclear_spin))
        number = num_coupled_atoms * 2 * spin + 1
        
        # Multiplicity label
        if number == 1:
            multiplicity = "singlet"
            symbol = "s"
        elif number == 2:
            multiplicity = "doublet"
            symbol = "d"
        elif number == 3:
            multiplicity = "triplet"
            symbol = "t"
        elif number == 4:
            multiplicity = "quartet"
            symbol = "q"
        elif number == 5:
            multiplicity = "pentet"
            symbol = "p"
        elif number == 6:
            multiplicity = "sextet"
            symbol = "sext"
        elif number == 7:
            multiplicity = "septet"
            symbol = "sept"
        elif number == 8:
            multiplicity = "octet"
            symbol = "oct"
        elif number == 9:
            multiplicity = "nonet"
            symbol = "non"
        else:
            multiplicity = "10"
            symbol = "10"
        
        return {"symbol": symbol, "number": number, "multiplicity": multiplicity}


class NMR(Result_object, Floatable_mixin):
    """
    A result object containing all the NMR related data for a single atom.
    
    For a given atom, this class will contain:
        - The chemical shielding of this atom (including a breakdown by tensors, if available).
        - The spin-spin coupling constants between this atom and all other atoms for which couplings were calculated (also with a breakdown by tensor, if available).
    """
    
    def __init__(self, atom, shielding, couplings):
        """
        :param atom: The atom these NMR parameters relate to.
        :param shielding: The chemical shielding object for this atom.
        :param couplings: A dictionary of all the coupling interactions calculated for this atom.
        """
        self.atom = atom
        self.shielding = shielding
        self.couplings = couplings
        
    def __float__(self):
        return float(self.shielding.isotropic())
    
#     def __eq__(self, other):
#         return self.atom == other.atom and self.shielding == other.shielding and self.couplings == other.couplings
    
    @classmethod
    def list_from_parser(self, parser):
        """
        """
        return [self(atom, parser.results.nmr_shieldings[atom], parser.results.nmr_couplings.between(atom)) for atom in parser.results.atoms if atom in parser.results.nmr_shieldings]
    
    @classmethod
    def list_from_dump(self, data, result_set, options):
        """
        Get a list of instances of this class from its dumped representation.
          
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        return [
            self(
                result_set.atoms.find(nmr_dict['atom']),
                NMR_shielding.from_dump(nmr_dict['shielding'], result_set, options),
                NMR_spin_couplings_list.from_dump(nmr_dict['couplings'], result_set, options)
            ) for nmr_dict in data
        ]
        
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "atom": self.atom.label,
            "shielding": self.shielding.dump(digichem_options, all),
            "couplings": self.couplings.dump(digichem_options, all)
        }
        

class NMR_tensor_ABC(Result_object):
    """ABC for classes that contain dicts of NMR tensors."""
    
    tensor_names = ()
    units = ""
    
    def __init__(self, tensors, total_isotropic = None):
        self.tensors = tensors
        
        # This is unused.
        self.total_isotropic = total_isotropic
    
    def eigenvalues(self, tensor = "total", real_only = True):
        """
        Calculate the eigenvalues for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        """
        try:
            return numpy.array([val.real for val in numpy.linalg.eigvals(self.tensors[tensor])])
        
        except KeyError:
            if tensor not in self.tensor_names:
                raise ValueError("The tensor '{}' is not recognised".format(tensor)) from None
            
            elif tensor not in self.tensors:
                raise ValueError("The tensor '{}' is not available".format(tensor)) from None
        
    def isotropic(self, tensor = "total"):
        """
        Calculate the isotropic value for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        """
        try:
            eigenvalues = self.eigenvalues(tensor)
            return sum(eigenvalues) / len(eigenvalues)

        except ValueError:
            if tensor == "total" and self.total_isotropic is not None:
                # Use the fallback.
                return self.total_isotropic

            else:
                raise
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        return {
            "tensors": {t_type: {"value": list(list(map(float, dim)) for dim in tensor), "units": self.units} for t_type, tensor in self.tensors.items()},
            "eigenvalues": {t_type: {"value": list(list(map(float, self.eigenvalues(t_type)))), "units": self.units} for t_type in self.tensors},
            "isotropic": {t_type: {"value": float(self.isotropic(t_type)), "units": self.units} for t_type in self.tensors}
        }


class NMR_shielding(NMR_tensor_ABC):
    """
    A result object to represent the chemical shielding of an atom.
    """
    
    tensor_names = ("paramagnetic", "diamagnetic", "total")
    units = "ppm"
    
    def __init__(self, tensors, reference = None, total_isotropic = None):
        """
        :param tensors: A dictionary of tensors.
        :param reference: An optional reference isotropic value to correct this shielding by.
        """
        super().__init__(tensors, total_isotropic)
        self.reference = reference
        
    def isotropic(self, tensor = "total", correct = True):
        """
        Calculate the isotropic value for a given tensor.
        
        :param tensor: The name of a tensor to calculate for (see tensor_names). Use 'total' for the total tensor.
        :param correct: Whether to correct this shielding value by the reference.
        """
        try:
            eigenvalues = self.eigenvalues(tensor)
            absolute = sum(eigenvalues) / len(eigenvalues)
        
        except ValueError:
            if tensor == "total" and self.total_isotropic is not None:
                # Use the fallback.
                absolute = self.total_isotropic
            
            else:
                raise None

        if correct and self.reference is not None:
            return self.reference - absolute
        else:
            return absolute
        
    @classmethod
    def dict_from_parser(self, parser):
        """
        Create a dict of NMR shielding objects from an output file parser.
        Each key is the atom being shielded.
        
        :param parser: An output file parser.
        :return: A list of NMR_shielding objects. The list will be empty if no NMR data is available.
        """
        shieldings = {}
        try:
            for atom_index, tensors in parser.data.nmrtensors.items():
                total_isotropic = tensors.pop("isotropic")
                shieldings[parser.results.atoms[atom_index]] = self(
                    tensors,
                    reference = parser.options['nmr']['standards'].get(parser.results.atoms[atom_index].element.symbol, None),
                    total_isotropic = total_isotropic
                )
        
        except AttributeError:
            return {}
        
        return shieldings
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dic = {
            "reference": self.reference,
        }
        
        dump_dic.update(super()._dump_(digichem_options, all))
        return dump_dic
    
    @classmethod
    def from_dump(self, data, result_set, options):
        tensors = {
            t_type: tensor_dict['value']
            for t_type, tensor_dict
            in data['tensors'].items()
        }
        return self(
            tensors,
            reference = data['reference']
        )


# We could look at some more advanced type of container for couplings.
# A dictionary might make sense, but it's difficult to choose how exactly the keys should be arranged and ordered.
class NMR_spin_couplings_list(Result_container):
    """A collection of NMR spin-spin couplings."""
    
    @classmethod
    def from_parser(self, parser):
        """
        Get an NMR_spin_couplings_list object from an output file parser.
        
        :param parser: An output file parser.
        :return: A NMR_spin_couplings_list object. The list will be empty if no NMR data is available.
        """
        return self(NMR_spin_coupling.list_from_parser(parser))
    
    @classmethod
    def from_dump(self, data, result_set, options):
        return self(NMR_spin_coupling.list_from_dump(data, result_set, options))
    
    def find(self, criteria = None, *, label = None, index = None):
        """
        """
        if label is None and index is None and criteria is None:
            raise ValueError("One of 'criteria', 'label' or 'index' must be given")
        
        if criteria is not None:
            if str(criteria).isdigit():
                index = int(criteria)
            
            else:
                label = criteria
        
        # Now get our filter func.
        if label is not None:
            filter_func = lambda coupling: label not in [atom.label for atom in coupling.atoms]
        
        elif index is not None:
            filter_func = lambda coupling: index not in [atom.index for atom in coupling.atoms]
        
        # Now search.
        found = type(self)(filterfalse(filter_func, self))
        
        if len(found) == 0:
            if label is not None:
                criteria_string = "label = '{}'".format(label)
            
            elif index is not None:
                criteria_string = "index = '{}'".format(index)
            
            raise Result_unavailable_error("NMR", "could not NMR data for atom '{}'".format(criteria_string))
        
        return found
    
    def between(self, atom1, atom1_isotopes = None, atom2 = None, atom2_isotopes = None):
        """
        Return a list containing all couplings involving either one or two atoms.
        """
        return type(self)([coupling for coupling in self
            if atom1 in coupling.atoms and
            (atom2 is None or atom2 in coupling.atoms) and
            (atom1_isotopes is None or coupling.isotopes[coupling.atoms.index(atom1)] in atom1_isotopes) and
            (atom2_isotopes is None or atom2 is None or coupling.isotopes[coupling.atoms.index(atom2)] in atom2_isotopes)
        ])

class NMR_spin_coupling(NMR_tensor_ABC):
    """
    A result object to represent spin-spin NMR couplings.
    """
    
    tensor_names = ("paramagnetic", "diamagnetic", "fermi", "spin-dipolar", "spin-dipolar-fermi", "total")
    units = "Hz"
    
    def __init__(self, atoms, isotopes, tensors, total_isotropic = None):
        """
        :param atoms: Tuple of atoms that this coupling is between.
        :param isotopes: Tuple of the specific isotopes of atoms.
        :param tensors: A dictionary of tensors.
        """
        super().__init__(tensors, total_isotropic = total_isotropic)
        self.atoms = atoms
        self.isotopes = isotopes
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a list of NMR coupling objects from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of NMR_spin_coupling objects. The list will be empty if no NMR data is available.
        """
        couplings = []
        try:
            for atom_tuple, isotopes in parser.data.nmrcouplingtensors.items():
                for isotope_tuple, tensors in isotopes.items():
                    total_isotropic = tensors.pop("isotropic")
                    couplings.append(
                        self(
                            (parser.results.atoms[atom_tuple[0]], parser.results.atoms[atom_tuple[1]]),
                            isotope_tuple,
                            tensors,
                            total_isotropic = total_isotropic
                        )
                    )
        
        except AttributeError:
            return []
        
        return couplings
    
    @classmethod
    def list_from_dump(self, data, result_set, options):
        return [
            self(
                tuple(result_set.atoms.find(atom_label) for atom_label in dump_dict['atoms']),
                tuple(dump_dict['isotopes']),
                {
                    t_type: tensor_dict['value']
                    for t_type, tensor_dict
                    in dump_dict['tensors'].items()
                }
            )
            for dump_dict
            in data
        ]
    
    def _dump_(self, digichem_options, all):
        """
        Get a representation of this result object in primitive format.
        """
        dump_dic = {
            "atoms": (self.atoms[0].label, self.atoms[1].label),
            "isotopes": (self.isotopes[0], self.isotopes[1]),
        }
        
        dump_dic.update(super()._dump_(digichem_options, all))
        return dump_dic