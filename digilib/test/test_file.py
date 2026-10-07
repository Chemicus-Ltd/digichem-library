import pytest
from pathlib import Path
import shutil

from digilib.file.fchk import Chk_to_fchk
from digilib.file.cube import Fchk_to_cube
from digilib.parse.util import open_for_parsing, parse_calculation
from digilib.exception import File_maker_exception

from digilib.test.util import digichem_options, data_directory
from digilib.test.test_result import gaussian_ES_result, turbomole_ES_result, orca_ES_result, \
    gaussian_opt_result, turbomole_opt_result, orca_opt_result, orca_opt_freq_result, \
    orca_nmr_result

@pytest.mark.skipif(not shutil.which("cubegen"),
                    reason="No cubegen available")
@pytest.mark.parametrize("result_path", [
    Path(data_directory(), "Pyridine/Gaussian 16 Excited States TDA Optimised S(1) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"),
], ids=["Pyridine"])
@pytest.mark.parametrize("sanity", [False, True], ids=["Normal", "Sanitize"])
def test_cube_from_fchk(result_path, sanity, tmp_path, digichem_options):
    """Can we make a cube file from an fchk file?"""
    with open_for_parsing(result_path) as open_files:
        result_set = parse_calculation(*open_files, options = digichem_options)
        maker = Fchk_to_cube.from_options(
            tmp_path / "tmp.cube",
            fchk_file = result_set.metadata.auxiliary_files['fchk_file'],
            options = digichem_options,
        )
        maker.sanitize = sanity

        assert maker.get_file().exists()

        # Check the cube is formatted correctly.
        with open(maker.get_file()) as cube:
            line = next(cube)
            line = next(cube)
            line = next(cube)
            split_line = line.split()
            if sanity:
                assert len(split_line) == 4
            
            else:
                assert len(split_line) == 5


def test_no_cubegen_error(digichem_options, tmp_path):
    """Do we throw the right exception if cubegen doesn't exist?"""
    maker = Fchk_to_cube.from_options(
        tmp_path / "tmp.cube",
        fchk_file = "foobar",
        options = digichem_options,
    )
    maker.cubegen_executable = "NOTREAL"

    with pytest.raises(File_maker_exception):
        maker.get_file()

def test_no_fchk_error(digichem_options, tmp_path):
    """Do we throw the right exception if we don't give an fchk file?"""
    maker = Fchk_to_cube.from_options(
        tmp_path / "tmp.cube",
        fchk_file = "foobar",
        options = digichem_options,
    )

    with pytest.raises(File_maker_exception):
        maker.get_file()

@pytest.mark.skipif(not shutil.which("formchk"),
                    reason="No formchk available")
def test_chk_to_fchk(digichem_options, tmp_path):
    """Can we convert a chk to an fchk file?"""
    maker = Chk_to_fchk.from_options(
        tmp_path / "tmp.fchk",
        chk_file = Path(data_directory(), "Chk/Pyridine.opt.chk"),
        options = digichem_options
    )
    assert maker.get_file().exists()


def test_no_formchk_error(digichem_options, tmp_path):
    """Do we throw the right exception if formchk doesn't exist?"""
    maker = Chk_to_fchk.from_options(
        tmp_path / "tmp.fchk",
        chk_file = Path(data_directory(), "Input/Benzene.chk"),
        options = digichem_options
    )
    maker.formchk_executable = "NOTREAL"

    with pytest.raises(File_maker_exception):
        maker.get_file()


def test_no_chk_error(digichem_options, tmp_path):
    """Do we throw the right exception if we don't give a chk file?"""
    maker = Chk_to_fchk.from_options(
        tmp_path / "tmp.fchk",
        chk_file = "foobar",
        options = digichem_options
    )

    with pytest.raises(File_maker_exception):
        maker.get_file()