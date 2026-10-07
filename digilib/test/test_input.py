"""Tests for parsing input files"""

import pytest
from pathlib import Path

from digilib.test.util import pyridine_si_v2, pyridine_si_v1, pyridine_cml,\
    result_files, ethane_xyz, benzene_cdx, cyclopentane_com, cyclopentane_iso_com, pyridinium_si
from digilib.input.digichem_input import si_from_file

@pytest.mark.parametrize("file_path", [
        pyridine_si_v2,
        pyridine_si_v1
     ])
def test_si_reading(file_path):
    """
    Test whether we can correctly read from an existing si file.
    """
    si_file = si_from_file(file_path)
    
    # Check charge and mult.
    assert si_file.charge == 0
    assert si_file.multiplicity == 1
    
    # Check geometry.
    # All files should contain data for pyridine.
    assert len(si_file.atoms) == 11
    assert str(si_file.formula) == "NC5H5"


def test_si_writing(tmp_path):
    """
    Test whether we can correctly write a new .si file.
    """
    # Read in an input file.
    si_file = si_from_file(pyridine_cml)
    
    # Write to our temp dir.
    new_si_path = Path(tmp_path, "Pyrdine.si")
    with open(new_si_path, "wt") as out_file:
        si_file.to_file(out_file)
        
    # Parse the newly written file.
    new_si_file = si_from_file(new_si_path)
    
    # Check they're the same.
    assert new_si_file == si_file


@pytest.mark.parametrize("file_path, sha", [
        [result_files['gaussian'][0], "ebd4f4e9f81cdb57ec2d2f2e1fba9cef0698f4c0"],
        [result_files['turbomole'][0], "60a8ebd9e701b849cfccd9cbb41684519a7fdf0b"],
        [result_files['orca'][0], "e48e7f653f4e67c1bd4c5c4bb76405fad2d441d0"],
     ])
def test_si_history(file_path, sha):
    """
    Test whether the history attribute is set properly.
    """
    si_file = si_from_file(file_path)
    
    assert si_file.history == sha
    
    dump = si_file.dict
    
    assert dump['history'] == sha


@pytest.mark.parametrize("file_path, format", [
    (cyclopentane_com, "mol"),
    (cyclopentane_iso_com, "mol"),
    (ethane_xyz, "cml"),
    (benzene_cdx, "xyz"),
    (result_files['gaussian'][0], "com")
])
def test_si_file_convert(file_path, format, tmp_path):
    """
    Test whether we can convert from arbitrary formats to another format.
    """
    si_file = si_from_file(file_path)
    si_file.to_format(format, (tmp_path / "file").with_suffix("." + format))


def test_si_charge_preserved():
    """
    Can we preserve charge and multiplicity when converting?
    """
    si_file = si_from_file(pyridinium_si)

    # Check charge and mult.
    assert si_file.charge == 1
    assert si_file.multiplicity == 1

    # Convert, using the in-built charge.
    com = si_file.to_format('com')
    assert com.splitlines()[5] == "1  1"


def test_si_charge_override():
    """
    Can we preserve charge and multiplicity when converting?
    """
    si_file = si_from_file(pyridinium_si, charge = 0, multiplicity = 2)

    # Check charge and mult.
    assert si_file.charge == 0
    assert si_file.multiplicity == 2

    # Convert, using the in-built charge.
    com = si_file.to_format('com')
    assert com.splitlines()[5] == "0  2"