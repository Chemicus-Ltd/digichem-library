"""Tests for translation classes."""
import pytest

from digilib.translate import Solvent, Basis_set, Multiplicity, Functional,\
    Cube_grid_points


def test_solvent():
    """
    Test conversions for DCM.
    """
    # Get some different ways to represent DCM.
    solvents = (Solvent("dcm"), Solvent("DCM"), Solvent("dichloroMethane"), Solvent("Dichloro Methane"))
    
    # Check all have equivalent string representations.
    assert all([str(solvents[0]) == str(solvent) for solvent in solvents])
    
    # And the correct gaussian name.
    assert all(["DiChloroMethane" == solvent.to_gaussian() for solvent in solvents])
    
    # And epsilon
    assert all([pytest.approx(8.93) == solvent.translate("epsilon") for solvent in solvents])
    

def test_631g_str_str():
    """
    Test conversions for basis sets.
    """
    # Test both types of 6-31G**.
    basis_sets = (Basis_set("6-31g(d,p)"), Basis_set("6-31G**"))
    
    # Check both are equivalent.
    assert str(basis_sets[0]) == str(basis_sets[1])
    
    # Check Gaussian version.
    assert all(["6-31G**" == basis_set.to_gaussian() for basis_set in basis_sets])
    
    # Check Turbomole version.
    assert all(["6-31G**" == basis_set.to_turbomole() for basis_set in basis_sets])
    
def test_def2_svp_p():
    """
    Test weird gaussian name for def2-SV(P).
    """
    assert Basis_set("def2-SV(P)").to_gaussian() == "def2SVPP"

@pytest.mark.parametrize("basis_set_name", ["def2-SVP", "def2-TZV", "def2-TZVP", "def2-QZV", "def2-QZVP"])
def test_karlsruhe(basis_set_name):
    """
    Test weird gaussian names for Karlsruhe basis sets.
    """
    assert "-" not in Basis_set(basis_set_name).to_gaussian()
    
def test_ccpvdz():
    assert Basis_set("cc-PvDz").to_gaussian() == "cc-pVDZ"
    assert str(Basis_set("cc-PvDz")) == "cc-pVDZ"
    

def test_b3lyp():
    """
    Test conversions for DFT functionals.
    """
    funcs = (Functional("b3-LYp"), Functional("B3lyp"))
    
    assert all(["b3-lyp" == func.to_turbomole() for func in funcs])
    
    assert all(["B3LYP" == func.to_gaussian() for func in funcs])


def test_pbe0gaussian():
    # Turbomole functionals are case sensitive. This is particularly annoying because functional
    # names are all lower-case, except for the word 'Gaussian'. Check this.
    assert Functional("PBE0 gaussian").to_turbomole() == "pbe0 Gaussian"


def test_multiplicity():
    """
    Test conversion for multiplicity.
    """
    mults = (Multiplicity(1), Multiplicity("1"), Multiplicity("singlet"), Multiplicity("sInGlEt"), Multiplicity("S"))
    
    # Check conversions.
    assert all(["Singlet" == mult.to_gaussian() for mult in mults])
    
    assert all([1 == mult.to_turbomole() for mult in mults])
    
    assert all(["S" == mult.symbol for mult in mults])
    
    assert all([1 == mult.number for mult in mults])
    
    assert all(["Singlet" == mult.string for mult in mults])
    
def test_grid_points():
    """
    Test different names for cube grid points.
    """
    # Check the special 'Default' value is correct.
    default = Cube_grid_points("default")
    
    assert default.translate("gaussian")    == 0
    assert default.translate("orca")        == 100
    assert default.translate("turbomole")   == "m3"
    
    # Check the named values.
    for name, value in [("Tiny", 25), ("Small", 50), ("Medium", 100), ("Large", 200), ("Huge", 500)]:
        assert Cube_grid_points(name).translate("points") == value
        
    # Check a random value.
    assert Cube_grid_points(123).translate("points") == 123



