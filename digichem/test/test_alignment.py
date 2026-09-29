"""Tests for checking the molecular alignment procedures"""
from pathlib import Path
import pytest

from digichem.parse import parse_calculation
from digichem.test.util import data_directory, result_files, digichem_options

@pytest.mark.parametrize(
    "method, x, y, z",
    [
        ("MIN",  6.7394, 4.9724, 5.650116197822169e-16),
        ("FAP",  7.183965325890709, 5.521277116003083, 4.126692319086836e-16),
        ("AA",   6.7394, 4.9724, 1.6950348593466506e-15),
        ("AAA",  7.069039617795138, 5.32635947371929, 3.7109659835279796e-16),
        ("GRID", 6.7394, 4.9724, 5.650116197822169e-16),
        ("NEST", 6.7394175185004, 4.972417415443474, 0.0001949261518701997),
        ("BAY",  6.7394, 4.9724, 5.650116197822169e-16)
    ]
)
def test_alignment_method(method, x, y, z, digichem_options):
    result = parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Single Point (Singlet) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options, ornt = method)

    # Check axes are in the correct order.
    assert result.atoms.X_length >= result.atoms.Y_length
    assert result.atoms.Y_length >= result.atoms.Z_length

    # Check we have the roughly expected value
    assert result.atoms.X_length == pytest.approx(x, rel=1e-5)
    assert result.atoms.Y_length == pytest.approx(y, rel=1e-5)
    assert result.atoms.Z_length == pytest.approx(z, rel=1e-5)


def test_fixed_alignment(digichem_options):
    result = parse_calculation(Path(data_directory(), "Naphthalene/Gaussian 16 Single Point (Singlet) PBE1PBE (GD3BJ) Toluene 6-31G(d,p).tar.gz"), options = digichem_options, ornt = 'FIX')

    # Uniquely, our axes are not in the correct order.
    assert result.atoms.Y_length >= result.atoms.Z_length
    assert result.atoms.Z_length >= result.atoms.X_length

    # Check we have the roughly expected value
    assert result.atoms.X_length == pytest.approx(0.0, rel=1e-5)
    assert result.atoms.Y_length == pytest.approx(6.7394, rel=1e-5)
    assert result.atoms.Z_length == pytest.approx(4.9724, rel=1e-5)
