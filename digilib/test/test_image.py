import pytest
from pathlib import Path
import shutil

from digilib.image.excited_states import Excited_states_diagram_maker
from digilib.image.graph import Convergence_graph_maker
from digilib.image.orbitals import Orbital_diagram_maker
from digilib.image.spectroscopy import Absorption_graph_maker, Emission_graph_maker, \
    Frequency_graph_maker, NMR_graph_maker, NMR_graph_zoom_maker
from digilib.image.structure import Skeletal_image_maker
import digilib.image.vmd
import digilib.image.render


from digilib.parse.util import open_for_parsing, parse_calculation

from digilib.test.util import digichem_options, data_directory
from digilib.test.test_result import gaussian_ES_result, turbomole_ES_result, orca_ES_result, \
    gaussian_opt_result, turbomole_opt_result, orca_opt_result, orca_opt_freq_result, \
    orca_nmr_result

HAS_VMD = shutil.which(digilib.config.get_config()['render']['vmd']['executable'])
HAS_BLENDER = shutil.which(digilib.config.get_config()['render']['batoms']['blender'])


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


@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.HOMO.cube")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Structure_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.vmd.Orbital_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Structure_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
    pytest.param(digilib.image.render.Orbital_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
])
def test_3d_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.SCF.cube")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Density_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Density_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_density_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""

    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Benzene.anion.spin.cube")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Spin_density_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Spin_density_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_spin_density_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("homo_cube, lumo_cube", [
    [Path(data_directory(), "Cubes/Pyridine.HOMO.cube"), Path(data_directory(), "Cubes/Pyridine.LUMO.cube")]
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Combined_orbital_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Combined_orbital_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_combined_orbital_image(homo_cube, maker_cls, lumo_cube, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        HOMO_cube_file = homo_cube,
        LUMO_cube_file = lumo_cube,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Benzene.anion.AHOMO.cube")
], ids = ["Gaussian"])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Alpha_orbital_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Alpha_orbital_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_unrestricted_orbital_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path, cube_file", [
    [  
        Path(data_directory(), "Pyridine/Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
        Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
    ]
])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Permanent_dipole_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Permanent_dipole_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_pdm_image(result_path, cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = maker_cls.from_options(
            tmp_path / "tmp.png",
            cube_file = cube_file,
            dipole_moment = result_set.dipole_moment,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("result_path, cube_file", [
    [
        Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
        Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
    ]
])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Transition_dipole_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Transition_dipole_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_tdm_image(result_path, cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = maker_cls.from_options(
            tmp_path / "tmp.png",
            cube_file = cube_file,
            dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.electric,
            magnetic_dipole_moment = result_set.excited_states.find("S(1)").transition_dipole_moment.magnetic,
            options = digichem_options
        )

        maker.get_file('x0y0z0')
        assert Path(tmp_path, "tmp.x0y0z0.png").exists()


@pytest.mark.parametrize("cube_file", [
    Path(data_directory(), "Cubes/Pyridine.HOMO.cube"),
])
@pytest.mark.parametrize("maker_cls", [
    pytest.param(digilib.image.vmd.Orbital_image_maker, marks = pytest.mark.skipif(not HAS_VMD, reason = "VMD is not available")),
    pytest.param(digilib.image.render.Orbital_image_maker, marks = pytest.mark.skipif(not HAS_BLENDER, reason = "Blender is not available")),
], ids = ["vmd", "batoms"])
def test_nto_image(cube_file, maker_cls, tmp_path, digichem_options):
    """Can we make a 3D structure image?"""
    maker = maker_cls.from_options(
        tmp_path / "tmp.png",
        cube_file = cube_file,
        options = digichem_options
    )

    maker.get_file('x0y0z0')
    assert Path(tmp_path, "tmp.x0y0z0.png").exists()