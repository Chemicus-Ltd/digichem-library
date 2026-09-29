"""Tests for functioning of oprattle"""
from pathlib import Path

import pytest

from digichem.file.prattle import Oprattle_formats, Openprattle_converter
from digichem.test.util import ethane_xyz, benzene_cdx, c6h6_cml, benzene_implicit_cml, benzene_explicit_cml, cyclohexane_implicit_cml, cyclohexane_explicit_cml, acetone_BH3_implicit_cml, acetone_BH3_explicit_cml, benzene_explicit_cdx, benzene_implicit_cdx

@pytest.mark.parametrize("readwrite", ["read", "write"])
def test_formats(readwrite):
    if readwrite == "read":
        assert len(Oprattle_formats().read()) > 0
    
    else:
        assert len(Oprattle_formats().write()) > 0


@pytest.mark.parametrize("file_path, size", [
    (ethane_xyz, 8),
    (benzene_cdx, 12)
])
def test_from_path(file_path, size):
    """
    Can we convert a file found on the filesystem?
    """
    output_data = Openprattle_converter(input_file_path = file_path).convert("com")
    assert len(output_data.splitlines()) == size +7


@pytest.mark.parametrize("file_path, mode, size", [
    (ethane_xyz, "r", 8),
    (benzene_cdx, "rb", 12)
])
def test_from_file(file_path, mode, size):
    """
    Can we convert from an open file?
    """
    with open(file_path, mode) as input_file:
        output_data = Openprattle_converter(input_file = input_file, input_file_type = Path(file_path).suffix[1:]).convert("com")
    
    assert len(output_data.splitlines()) == size +7


@pytest.mark.parametrize("file_path, mode, size", [
    (ethane_xyz, "r", 8),
    (benzene_cdx, "rb", 12)
])
def test_from_buffer(file_path, mode, size):
    """
    Can we convert from a string
    """
    with open(file_path, mode) as input_file:
        data = input_file.read()
    
    output_data = Openprattle_converter(input_file_buffer = data, input_file_type = Path(file_path).suffix[1:]).convert("com")
    assert len(output_data.splitlines()) == size +7


@pytest.mark.parametrize("file_path, mode, size", [
    (ethane_xyz, "r", 8),
    (benzene_cdx, "rb", 12)
])
def test_file_to_buffer(file_path, mode, size):
    """
    Can we convert to a string
    """
    with open(file_path, mode) as input_file:
        output_data = Openprattle_converter(input_file = input_file, input_file_type = Path(file_path).suffix[1:]).convert("xyz")
    assert isinstance(output_data, str)
    assert len(output_data.splitlines()) == size +2


@pytest.mark.parametrize("file_path, mode, size", [
    (ethane_xyz, "r", 8),
    (benzene_cdx, "rb", 12)
])
def test_buffer_to_buffer(file_path, mode, size):
    """
    Can we convert to a string
    """
    with open(file_path, mode) as input_file:
        input_data = input_file.read()
    
    output_data = Openprattle_converter(input_file_buffer = input_data, input_file_type = Path(file_path).suffix[1:]).convert("xyz")
    assert isinstance(output_data, str)
    assert len(output_data.splitlines()) == size +2


@pytest.mark.parametrize("file_path, size", [
    (ethane_xyz, 8),
    (benzene_cdx, 12)
], ids = ["ethane", "benzene"])
def test_gen3d(file_path, size):
    first = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = False)

    same = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = False)
    different = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = True)

    assert len(first.splitlines()) == size +2
    assert len(same.splitlines()) == size +2
    assert len(different.splitlines()) == size +2
    assert first == same
    assert first != different


@pytest.mark.parametrize("file_path, size", [
    (c6h6_cml, 12),
    (benzene_explicit_cml, 12),
    (benzene_implicit_cml, 12),
    (cyclohexane_explicit_cml, 18),
    (cyclohexane_implicit_cml, 18),
    (acetone_BH3_explicit_cml, 14),
    # Obabel can't handle this one sadly
    #(acetone_BH3_implicit_cml, 14),
    (benzene_explicit_cdx, 12),
    (benzene_implicit_cdx, 12),
], ids = ["c6h6", "benzene_explicit_cml", "benzene_implicit_cml", "cyclohexane_explicit", "cyclohexane_implicit", "acetone_BH3_explicit", "benzene_explicit_cdx", "benzene_implicit_cdx"])
def test_hidden_h(file_path, size):
    """Can we correctly determine the number of implicit Hs in a structure"""
    output_data = Openprattle_converter(input_file_path = file_path).convert("xyz", gen3D = False)

    # Check we haven't suddenly grown more Hs.
    assert len(output_data.splitlines()) == size +2