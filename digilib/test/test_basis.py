"""Tests for handling basis sets from BSE."""
import pytest

from digilib.basis import BSE_basis_set

def test_basis():
    """Can we load a simple basis set?"""
    basis = BSE_basis_set({
        "STO-3G":"all"
    })
    basis.to_format()

def test_ecps():
    """Can we load a simple basis set eith ecps?"""
    basis = BSE_basis_set({
        "LANL2DZ":"all"
    })
    assert basis.has_ECPs()

@pytest.mark.parametrize("set", [
    {"STO-3G":"all"},
    {"STO-3G":"C"},
    {"STO-3G":"C,H"},
    {"STO-3G":"H","def2-SVP":"C"}
], ids = ["all", "C", "C,H", "mixed"])
def test_meta(set):
    """Can we handle basis set metadata correctly?"""
    basis = BSE_basis_set(set)
    str(basis)

@pytest.mark.parametrize("filter, number", [
    ("C", 1),
    ("C,H", 2),
    (("C", "H"), 2)
], ids = ["C", "C,H", list])
def test_filter(filter, number):
    """Can we return only parts of a basis set?"""
    basis = BSE_basis_set({
        "STO-3G":"all"
    })
    output = basis.to_format("gaussian94", filter)
    assert output.count("****") == number

def test_mixed_basis():
    """Can we combine multiple basis sets?"""
    basis = BSE_basis_set({
        "STO-3G":"H",
        "def2-SVP" : "C"
    })

    basis.to_format("gaussian94")

def test_mixed_basis_filter():
    """Can we combine multiple basis sets?"""
    basis = BSE_basis_set({
        "STO-3G":"H",
        "def2-SVP" : "C"
    })

    output = basis.to_format("gaussian94", "C,H")
    assert output.count("****") == 2

def test_wrong_basis():
    """Check a non-existent basis set is flagged"""
    basis = BSE_basis_set({
        "STO-3G":"H",
        "defs2-SVP" : "C"
    })

    with pytest.raises(KeyError):
        basis.to_format("gaussian94")