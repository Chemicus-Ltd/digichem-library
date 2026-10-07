import pytest
from pathlib import Path

from digilib.image.excited_states import Excited_states_diagram_maker
from digilib.image.graph import Convergence_graph_maker
from digilib.image.orbitals import Orbital_diagram_maker
from digilib.image.spectroscopy import Absorption_graph_maker, Emission_graph_maker, \
    Frequency_graph_maker, NMR_graph_maker, NMR_graph_zoom_maker
from digilib.image.structure import Skeletal_image_maker
from digilib.image.vmd import Structure_image_maker, Orbital_image_maker, Density_image_maker, \
    Combined_orbital_image_maker, Permanent_dipole_image_maker, Transition_dipole_image_maker
from digilib.file.cube import Fchk_to_cube, Fchk_to_density_cube, Fchk_to_nto_cube, Fchk_to_spin_cube
from digilib.parse.util import open_for_parsing, parse_calculation

from digilib.test.util import digichem_options, data_directory
from digilib.test.test_result import gaussian_ES_result, turbomole_ES_result, orca_ES_result, \
    gaussian_opt_result, turbomole_opt_result, orca_opt_result, orca_opt_freq_result, \
    orca_nmr_result

pytest.skip(allow_module_level=True)

@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_es_diagram(result_set, tmp_path, digichem_options):
    """Can we make an excited states diagram?"""
    maker = Excited_states_diagram_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        ground_state = result_set.ground_state,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_result")
])
def test_convergence_diagram(result_set, tmp_path, digichem_options):
    """Can we make an optimisation graph?"""
    maker = Convergence_graph_maker.from_options(
        tmp_path / "tmp.png",
        energies = result_set.energies.scf,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_result")
])
def test_orbital_diagram(result_set, tmp_path, digichem_options):
    """Can we make an orbital diagram?"""
    maker = Orbital_diagram_maker.from_options(
        tmp_path / "tmp.png",
        orbitals = result_set.orbitals,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_abs_diagram(result_set, tmp_path, digichem_options):
    """Can we make an absorption spectrum?"""
    maker = Absorption_graph_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_ES_result"),
    pytest.lazy_fixture("turbomole_ES_result"),
    pytest.lazy_fixture("orca_ES_result")
])
def test_pl_diagram(result_set, tmp_path, digichem_options):
    """Can we make an emission spectrum?"""
    maker = Emission_graph_maker.from_options(
        tmp_path / "tmp.png",
        excited_states = result_set.excited_states,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_freq_result")
])
def test_freq_diagram(result_set, tmp_path, digichem_options):
    """Can we make an IR spectrum?"""
    maker = Frequency_graph_maker.from_options(
        tmp_path / "tmp.png",
        vibrations = result_set.vibrations,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("orca_nmr_result"),
])
def test_nmr_diagram(result_set, tmp_path, digichem_options):
    """Can we make an NMR spectrum?"""
    # Get a suitable spectrometer object.
    spectrometer = result_set.nmr.spectrometer

    maker = NMR_graph_maker.from_options(
        tmp_path / "tmp.png",
        graph = spectrometer.get_graph(1, 1),
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("orca_nmr_result"),
])
def test_focus_nmr_diagram(result_set, tmp_path, digichem_options):
    """Can we make an NMR spectrum?"""
    # Get a suitable spectrometer object.
    spectrometer = result_set.nmr.spectrometer

    maker = NMR_graph_zoom_maker.from_options(
        tmp_path / "tmp.png",
        graph = spectrometer.get_graph(1, 1),
        focus = result_set.atoms.groups[(2, "H")],
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_set", [
    pytest.lazy_fixture("gaussian_opt_result"),
    pytest.lazy_fixture("turbomole_opt_result"),
    pytest.lazy_fixture("orca_opt_freq_result")
])
def test_2d_diagram(result_set, tmp_path, digichem_options):
    """Can we make an IR spectrum?"""
    maker = Skeletal_image_maker.from_options(
        tmp_path / "tmp.png",
        atoms = result_set.atoms,
        options = digichem_options
    )

    maker.get_image()
    assert Path(tmp_path, "tmp.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    Structure_image_maker,
    Orbital_image_maker,
])
def test_3d_image(result_path, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = maker_cls.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


# @pytest.mark.parametrize("result_path", [
#     Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
# ], ids = ["Gaussian"])
# @pytest.mark.parametrize("density", [
#     "SCF"
# ])

@pytest.mark.parametrize("result_path, density", [
    [Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), "SCF"]
], ids = ["Gaussian"])
def test_density_image(result_path, density, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Density_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_density_cube.from_options(
                tmp_path / "tmp.SCF.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                density_type = density,
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path, density", [
    [Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), "SCF"]
], ids = ["Gaussian"])
def test_spin_density_image(result_path, density, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Density_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_spin_cube.from_options(
                tmp_path / "tmp.Spin.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                spin_density = density,
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
], ids = ["Gaussian"])
def test_combined_orbital_image(result_path, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Combined_orbital_image_maker.from_options(
            tmp_path / "tmp.png",
            HOMO_cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.HOMO.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                orbital = "LUMO",
                options = digichem_options
            ),
            LUMO_cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.LUMO.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                orbital = "LUMO",
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Benzene Anion/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Gas Phase 6-31G(d,p).tar.gz")
], ids = ["Gaussian"])
def test_unrestricted_orbital_image(result_path, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Orbital_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                cubegen_type = "AMO",
                orbital = "HOMO",
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz")
])
def test_pdm_image(result_path, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Permanent_dipole_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                options = digichem_options
            ),
            dipole_moment = result_set.dipole_moment,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
])
def test_tdm_image(result_path, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Transition_dipole_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_cube.from_options(
                tmp_path / "tmp.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                options = digichem_options
            ),
            dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.electric,
            magnetic_dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.magnetic,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
])
def test_nto_image(result_path, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Orbital_image_maker.from_options(
            tmp_path / "tmp.png",
            cube_file = Fchk_to_nto_cube.from_options(
                tmp_path / "tmp.cube",
                fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
                options = digichem_options
            ),
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()