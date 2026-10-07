"""Test for memory representations."""

import pytest

from digilib.memory import Memory


@pytest.mark.parametrize("value, valid", [
        ["1", True],
        [1, True],
        [1.0, True],
        ["1 B", True],
        ["1 KB", True],
        ["1 MB", True],
        ["1 GB", True],
        ["1 TB", True],
        ["foo bar", False],
     ])
def test_validation(value, valid):
    try:
        mem = Memory(value)

        str(mem)
    
    except Exception:
        if valid:
            raise
        
        else:
            return
    
    if not valid:
        raise