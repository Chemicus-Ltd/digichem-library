"""Tests for checking the values of parsed results."""

import pytest
from pathlib import Path
import scipy.constants

from digichem.parse import parse_calculation
from digichem.parse.util import parse_and_merge_calculations
from digichem.test.util import data_directory, digichem_options
from digichem.test.util import check_float_list, check_dipole, check_orbitals


@pytest.fixture(scope="module")
def gaussian_SP_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Single Point (Singlet) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_opt_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_radical_anion_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Benzene Anion/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Gas Phase 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_ES_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_emission_result(digichem_options):
    return parse_and_merge_calculations(
        Path(data_directory(), "Naphthalene/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)"),
        Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
        options = digichem_options
    )

@pytest.fixture(scope="module")
def gaussian_opt_ES_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_PDM_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_TDM_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def gaussian_PDM_ES_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options)

#############
# Turbomole #
#############

@pytest.fixture(scope="module")
def turbomole_sp_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Turbomole Single Point DFT PBE0 (GD3BJ) Toluene Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_opt_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ES_singlets_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ES_triplets_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ES_result(digichem_options):
    return parse_and_merge_calculations(
        Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**.tar.gz"),
        Path(data_directory(), "Naphthalene/Turbomole Excited States TDA 10 Triplets PBE0 (GD3BJ) 6-31G**.tar.gz"),
        options = digichem_options
    )

@pytest.fixture(scope="module")
def turbomole_PDM_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_TDM_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Turbomole Excited States TDA 10 Singlets PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_radical_anion_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Benzene Anion/Turbomole Optimisation Frequency PBE0 (GD3BJ) 6-31G**.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ADC2_opt_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Optimisation ADC(2) cc-pVDZ.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ADC2_singlets_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States ADC(2) S(1) and S(2) cc-pVDZ.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def turbomole_ADC2_triplets_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Naphthalene/Turbomole Excited States ADC(2) T(1) and T(2) cc-pVDZ.tar.gz"), options = digichem_options)

########
# ORCA #
########

@pytest.fixture(scope="module")
def orca_SP_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Single Point PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_solvent_SP_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Single Point PBE0 (GD3BJ) Toluene Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_grad_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Gradient PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_opt_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Optimisation PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_opt_freq_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Optimisation Frequencies PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_freq_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Frequencies PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_nmr_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca NMR PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)

@pytest.fixture(scope="module")
def orca_ES_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Pyridine/Orca Excited States TDA 10 Singlets 10 Triplets PBE0 (GD3BJ) Gas Phase Pople Basis Sets STO-3G.tar.gz"), options = digichem_options)


@pytest.mark.parametrize("result_set, num, final", [
        (pytest.lazy_fixture("gaussian_SP_result"), 1, -10488.990333747),
        (pytest.lazy_fixture("gaussian_opt_result"), 5, -10488.990333747),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, -10488.990333747),
        (pytest.lazy_fixture("gaussian_opt_ES_result"), 4, -10488.883505795338),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 7, -10432.311381999756),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 1, -10432.311386936175),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 1, -10432.311386936175)
    ])
def test_SCF_energy(result_set, num, final):
    """Test the parsed energy is correct"""
    assert result_set.energies.scf.final == pytest.approx(final, abs=1e-5)
    
    # Check length, which will be 1 for SP, and >1 for the opts.
    assert len(result_set.energies.scf) == num


@pytest.mark.parametrize("result_set, num, final", [
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 7, -10467.158186939303),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 1, -10467.158181211305),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 1, -10467.158181211305)
    ])
def test_MP_energy(result_set, num, final):
    """Test the parsed energy is correct"""
    assert result_set.energies.mp.final == pytest.approx(final, abs=1e-5)
    
    # Check length, which will be 1 for SP, and >1 for the opts.
    assert len(result_set.energies.scf) == num
    assert len(result_set.energies.mp) == num


@pytest.mark.parametrize("result_set, charge, mult, energy", [
        (pytest.lazy_fixture("gaussian_opt_result"), 0, 1, -10488.990333747),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 0, 1, -10467.158186939303),
        # Radical anions.
        (pytest.lazy_fixture("gaussian_radical_anion_result"), -1, 2, -6310.5380531853125),
        (pytest.lazy_fixture("turbomole_radical_anion_result"), -1, 2, -6310.4572047072043),
    ])
def test_ground_state(result_set, charge, mult, energy):
    """Test the ground state properties"""
    
    assert result_set.ground_state.charge == charge
    assert result_set.ground_state.multiplicity == mult
    # Ground states have zero excited state energy by definition.
    assert result_set.ground_state.energy == 0.0
    # This is the total energy.
    assert result_set.ground_state.absolute_energy == pytest.approx(energy, abs=1e-5)


@pytest.mark.parametrize("result_set", [
        pytest.lazy_fixture("gaussian_SP_result"),
        pytest.lazy_fixture("gaussian_opt_result"),
        pytest.lazy_fixture("gaussian_ES_result"),
        pytest.lazy_fixture("gaussian_opt_ES_result"),
        pytest.lazy_fixture("turbomole_ADC2_opt_result"),
        pytest.lazy_fixture("turbomole_ADC2_singlets_result"),
        pytest.lazy_fixture("turbomole_ADC2_triplets_result"),
    ])
def test_atoms(result_set):
    """Test the parsed (unaligned) atoms are correct."""
    # The molecule is Naphthalene; C10H8
    # Check num atoms.
    assert len(result_set.raw_atoms) == 18
    
    # Check atom types.
    assert result_set.raw_atoms.element_dict['C'] == 10
    assert result_set.raw_atoms.element_dict['H'] == 8
    
    # Check mass.
    assert result_set.raw_atoms.molar_mass == pytest.approx(128.1706)
    
    # We don't check positions here because these can vary from calc to calc due to reorientation...


@pytest.mark.parametrize("result_set, length, width, height", [
        (pytest.lazy_fixture("gaussian_SP_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("gaussian_opt_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("gaussian_ES_result"), 6.739394, 4.972394, 0.000004),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 6.8000171599, 5.0161488999, 0.0000475622),
        (pytest.lazy_fixture("turbomole_ADC2_singlets_result"), 6.8000171599, 5.0161488999, 0.0000475622),
        (pytest.lazy_fixture("turbomole_ADC2_triplets_result"), 6.8000171599, 5.0161488999, 0.0000475622),
    ])
def test_alignment(result_set, length, width, height):
    """Test the aligned atom positions are correct."""
    # Rounding errors plague parsing and conversion of atomic coordinates.
    # Much of this error comes from the accuracy printed by the CC programs themselves,
    # so there is not much we can do about it.
    tolerance = 1e-5
    
    # Dimensions in angstrom.
    try:
        linear_ratio = 1-(width/length)
    except (FloatingPointError, ZeroDivisionError):
        linear_ratio = 1.0
        
    try:
        planar_ratio = 1-(height/width)
    except (FloatingPointError, ZeroDivisionError):
        planar_ratio = 1.0
    
    # Check dimensions.
    assert result_set.atoms.X_length == pytest.approx(length, abs=tolerance)
    assert result_set.atoms.Y_length == pytest.approx(width, abs=tolerance)
    assert result_set.atoms.Z_length == pytest.approx(height, abs=tolerance)
    
    # Check params.
    assert result_set.atoms.get_linear_ratio() == pytest.approx(linear_ratio, abs=tolerance)
    assert result_set.atoms.get_planar_ratio() == pytest.approx(planar_ratio, abs=tolerance)


@pytest.mark.parametrize("result_set, num_occ, num_unocc, homo, lumo", [
        (pytest.lazy_fixture("gaussian_opt_result"), 34, 156, -6.13072481, -0.92437071),
        (pytest.lazy_fixture("turbomole_ADC2_opt_result"), 34, 146, -7.78345591, 2.37048311)
    ])
def test_orbitals(result_set, num_occ, num_unocc, homo, lumo):
    """Test the energies of the parsed MOs"""
    check_orbitals(result_set.orbitals, num_occ, num_unocc, homo, lumo)
    
    # Check spin labels
    assert result_set.orbitals.spin_type == "none"


@pytest.mark.parametrize("result_set, num_occ, num_unocc, homo, lumo", [
        (pytest.lazy_fixture("gaussian_radical_anion_result"), (22, 21), (98, 99), (3.832451319, -0.0157826027), (6.268142298, 6.34270149)),
        (pytest.lazy_fixture("turbomole_radical_anion_result"), (22, 21), (92, 93), (3.825336056180260069, 0.0034865153882969327284), (6.275672672175910627, 6.387468734183322283))
    ])
def test_unrestricted_orbitals(result_set, num_occ, num_unocc, homo, lumo):
    """Test the energies of unrestricted orbitals"""
    # Check overall length.
    assert len(result_set.orbitals) == len(result_set.beta_orbitals)
    
    # Check each set of orbitals.
    check_orbitals(result_set.orbitals, num_occ[0], num_unocc[0], homo[0], lumo[0])
    check_orbitals(result_set.beta_orbitals, num_occ[1], num_unocc[1], homo[1], lumo[1])
    
    # Check spin labels
    assert result_set.orbitals.spin_type == "alpha"
    assert result_set.beta_orbitals.spin_type == "beta"


@pytest.mark.parametrize("result_set, num, index, frequency, intensity", [
        (pytest.lazy_fixture("gaussian_opt_result"), 48, 13, 806.2853, 121.0015),
        (pytest.lazy_fixture("gaussian_SP_result"), 0, None, None, None),
        (pytest.lazy_fixture("turbomole_opt_result"), 48, 13, 805.28, 103.77),
    ])
def test_frequencies(result_set, num, index, frequency, intensity):
    """Test the vibrational frequencies are correct"""
    
    # Check lengths.
    assert len(result_set.vibrations) == num
    assert len(result_set.vibrations.negative_frequencies) == 0

    if index != None:
        # Check one of the energies.
        assert result_set.vibrations[index].frequency == pytest.approx(frequency)
        assert result_set.vibrations[index].intensity == pytest.approx(intensity)
        
        # Check the frequencies are ordered.
        indices = [freq.level for freq in result_set.vibrations]
        assert indices == list(range(1, num+1))
        
    else:
        with pytest.raises(IndexError):
            result_set.vibrations[-1].frequency


@pytest.mark.parametrize("result_set, coords, axis_angle, plane_angle", [
        (pytest.lazy_fixture("gaussian_PDM_result"), (0.0, 2.5103, 0.0), 90.0, 0.0),
        (pytest.lazy_fixture("gaussian_PDM_ES_result"), (0.0001, -0.6147, 0.0001), 90.0, 0.0),
        (pytest.lazy_fixture("turbomole_PDM_result"), (-0.000743401, -2.20405438, -0.000494493), 90.0, 0.0),
    ])
def test_pdm(result_set, coords, axis_angle, plane_angle):
    """Test the permanent dipole moment"""
     
    check_dipole(result_set.pdm, coords)
    
    # Check angles
    assert float(result_set.pdm.X_axis_angle) == pytest.approx(axis_angle, abs=1e-1)
    assert float(result_set.pdm.XY_plane_angle) == pytest.approx(plane_angle, abs=1e-1)


@pytest.mark.parametrize("result_set, number, TDM", [
        (pytest.lazy_fixture("gaussian_TDM_result"), 20, (0.0, -0.0003, 0.5879)),
        (pytest.lazy_fixture("turbomole_TDM_result"), 10, (-0.000628, -0.000132, 0.558096))
    ])
def test_tedm(result_set, number, TDM):
    """Test transition dipole moments"""
    # Check number of dipoles.
    assert len([excited_state.transition_dipole_moment for excited_state in result_set.excited_states if excited_state.transition_dipole_moment.electric is not None]) == number
    
    # Check the S1 moments are correct.
    S1_TDM = result_set.excited_states.get_state("S(1)").transition_dipole_moment
    
    # Check coords and magnitude.
    check_dipole(S1_TDM.electric, TDM)


@pytest.mark.parametrize("result_set, number, TDM", [
        (pytest.lazy_fixture("gaussian_TDM_result"), 20, (0.6909, 0.0, 0.0)),
        (pytest.lazy_fixture("turbomole_TDM_result"), 10, (0.673218, -0.000274, -0.000274))
    ])
def test_tmdm(result_set, number, TDM):
    """Test transition dipole moments"""
    
    # Check number of dipoles.
    assert len([excited_state.transition_dipole_moment for excited_state in result_set.excited_states if excited_state.transition_dipole_moment.magnetic is not None]) == number
    
    # Check the S1 moments are correct.
    S1_TDM = result_set.excited_states.get_state("S(1)").transition_dipole_moment
    
    # Check coords and magnitude.
    check_dipole(S1_TDM.magnetic, TDM, abs = 1e-3)


@pytest.mark.parametrize("result_set, state, g_lum", [
        (pytest.lazy_fixture("gaussian_TDM_result"), "S(5)", -0.995)
    ])
def test_g_lum(result_set, state, g_lum):
    """Test calculation of G(lum) values for circularly polarized light."""
    
    TDM = result_set.excited_states.get_state(state).transition_dipole_moment
    
    TEDM = TDM.electric.gaussian_cgs_vector
    TMDM = TDM.magnetic.gaussian_cgs_vector
    
    # Test values in Gaussian-CGS units.))
    check_float_list(TEDM, tuple(254.35e-20 * i for i in TEDM))
    check_float_list(TMDM, tuple(-0.97401e-20 * i for i in TMDM))
    
    # Check angle.
    angle = (TEDM[0] * TMDM[0] + TEDM[1] * TMDM[1] + TEDM[2] * TMDM[2]) / ( ((TEDM[0]**2 + TEDM[1] **2 + TEDM[2] **2)**0.5) * ((TMDM[0]**2 + TMDM[1] **2 + TMDM[2] **2)**0.5) )
    assert TDM.cos_angle() == pytest.approx(angle)
    
    # Check g_lum
    assert TDM.g_value == pytest.approx(g_lum, abs=1e-3)


@pytest.mark.parametrize("result_set, num_singlets, num_triplets, dest, state_label, state_index, state_symmetry, state_energy, state_f", [
        (pytest.lazy_fixture("gaussian_ES_result"), 10, 10, 1.6231, "S(2)", 7, "Singlet-B1u", 4.7387, 0.1168),
        (pytest.lazy_fixture("turbomole_ES_result"), 10, 10, 1.6298, "S(2)", 7, "Singlet-A", 4.790868037954573, 0.0805)
    ])
def test_excited_states(result_set, num_singlets, num_triplets, dest, state_label, state_index, state_symmetry, state_energy, state_f):
    """Test excited states"""
    
    # Check number of states.
    assert result_set.excited_states.num_singlets == num_singlets
    assert result_set.excited_states.num_triplets == num_triplets
    
    # Check singlet-triplet splitting energy (dE(ST)).
    assert result_set.excited_states.singlet_triplet_energy == pytest.approx(dest, abs=1e-4)
    
    # Now check a specific state.
    state = result_set.excited_states.get_state(state_label)
    assert state.level == state_index
    assert state.symmetry == state_symmetry
    assert state.energy == pytest.approx(state_energy, abs=1e-4)
    assert state.oscillator_strength == pytest.approx(state_f, abs=1e-4)
    
    # Check wavelength.
    # e in J
    joule_energy = state_energy * scipy.constants.physical_constants['electron volt'][0]
    # wavelength in m
    meter_wavelength = (scipy.constants.Planck * scipy.constants.c) / joule_energy
    # wavelength in nm
    wavelength = meter_wavelength *1000000000
    
    assert state.wavelength == pytest.approx(wavelength, abs=1e-4)


@pytest.mark.parametrize("result_set, state_index, transition_index, start_orbital, end_orbital, coefficient", [
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 1, 34, 35, 0.957238734),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 2, 32, 37, 0.186520627),
        (pytest.lazy_fixture("gaussian_ES_result"), 1, 3, 33, 36, 0.175928167),
        (pytest.lazy_fixture("turbomole_ES_result"), 2, 1, 33, 35, 0.772657751),
        (pytest.lazy_fixture("turbomole_ES_result"), 2, 2, 34, 36, 0.618869938),
    ])
def test_excited_state_transitions(result_set, state_index, transition_index, start_orbital, end_orbital, coefficient):
    """Test excited state transitions."""
    
    state = result_set.excited_states[state_index -1]
    transition = state.transitions[transition_index -1]
    
    assert transition.level == transition_index
    assert transition.starting_mo.level == start_orbital
    assert transition.ending_mo.level == end_orbital
    assert transition.coefficient == pytest.approx(coefficient)
    assert transition.probability == pytest.approx(coefficient **2)


@pytest.mark.parametrize("result_set, transition_type, multiplicity, emission_type, energy, excited_energy, oscillator_strength, rate", [
        (pytest.lazy_fixture("gaussian_emission_result"), "adiabatic", 1, "fluorescence", 4.542828007, -10484.452882954, 0.0, 312658.74889),
        (pytest.lazy_fixture("gaussian_emission_result"), "vertical", 1, "fluorescence", 4.4360, -10484.452882954, 0.0, 291116.171714),
    ])
def test_emission(result_set, transition_type, multiplicity, emission_type, energy, excited_energy, oscillator_strength, rate):
    # TODO: Testing of fluorescence rate should be moved to another module (where other calculated properties can be tested).
    if transition_type == "adiabatic":
        emission = result_set.emission.adiabatic[multiplicity]
        
    else:
        emission = result_set.emission.vertical[multiplicity]
        
    # Test some general properties.
    assert emission.transition_type == transition_type
    assert emission.emission_type == emission_type
    
    # Test numerical properties.
    assert emission.energy == pytest.approx(energy)
    assert emission.excited_energy == pytest.approx(excited_energy)
    assert emission.oscillator_strength == pytest.approx(oscillator_strength)
    assert emission.emission_rate == pytest.approx(rate)


def test_h_nmr_standard(orca_nmr_result):
    reference = 31.68766666666667
    # Check the reference standard is available.
    assert all([nmr_result.shielding.reference == pytest.approx(reference, abs = 1e-4) for nmr_result in orca_nmr_result.nmr if nmr_result.atom.element.symbol == "H"])

def test_c_nmr_standard(orca_nmr_result):
    reference = 197.90316666666664
    # Check the reference standard is available.
    assert all([nmr_result.shielding.reference == pytest.approx(reference, abs = 1e-4) for nmr_result in orca_nmr_result.nmr if nmr_result.atom.element.symbol == "C"])

def test_nmr_h_isotope_options(orca_nmr_result):
    # Check we have the right options set for H.
    options = orca_nmr_result.nmr.spectrometer.isotope_options(1, 1)
    assert options['frequency'] == 400
    assert options['fwhm'] == 0.0015
    assert options['gaussian_resolution'] == 0.0001
    assert options['coupling_filter'] == 0.001
    assert options['pre_merge'] == 0.0005

def test_nmr_c_isotope_options(orca_nmr_result):
    # Check we have the right options set for H.
    options = orca_nmr_result.nmr.spectrometer.isotope_options(6, 13)
    assert options['frequency'] == 100.6

@pytest.mark.parametrize("result_set", [
        pytest.lazy_fixture("gaussian_SP_result"),
        pytest.lazy_fixture("gaussian_emission_result"),
        pytest.lazy_fixture("turbomole_sp_result"),
        pytest.lazy_fixture("orca_solvent_SP_result"),
    ])
def test_metadata_solvent(result_set):
    """Can we parse solution metadata correctly?"""
    assert result_set.metadata.solvent.name == "Toluene"
    # A bit of a rudimentary test, but there's already more advanced testing in cclib
    assert result_set.metadata.solvent.model is not None
    assert result_set.metadata.solvent.params['epsilon'] == pytest.approx(2.374, abs = 1e-3)
