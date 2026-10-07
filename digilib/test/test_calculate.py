"""Tests for calculated results"""

import pytest
from pathlib import Path
import shutil

from digilib.test.util import data_directory, digichem_options
from digilib.parse.util import parse_calculation


@pytest.fixture(scope="module")
def formaldehyde_SOC_result(digichem_options):
    return parse_calculation(Path(data_directory(), "Formaldehyde/Gaussian 09 Excited States TD-DFT 4 Singlets 4 Triplets B3LYP Gas Phase TZVP.tar.gz"), options = digichem_options)

@pytest.mark.skipif(not shutil.which("rwfdump"), reason="No rwfdump available")
@pytest.mark.parametrize("result_set, singlet, triplet, SOC", [
    # These SOC values are taken from the original paper J. Chem. Theory Comput. 2017, 13, 515−524 for PySOC,
    # the program Silioc uses to calculate SOC. We are limited in terms of accuracy by the number of decimals
    # (0) that they report...
        (pytest.lazy_fixture("formaldehyde_SOC_result"), "S(1)", "T(1)", 0.0),
        (pytest.lazy_fixture("formaldehyde_SOC_result"), "S(1)", "T(2)", 45.0)
    ])
def test_soc(result_set, singlet, triplet, SOC):
    """Test calculation of spin-orbit coupling values."""
    
    # Check we have SOC between each singlet state (including the ground) and each triplet state.
    assert len(result_set.soc) == (result_set.excited_states.num_singlets +1) * result_set.excited_states.num_triplets
    
    # Check some SOC values.
    assert result_set.soc.between(singlet, triplet).wavenumbers == pytest.approx(SOC, abs=1e-0)