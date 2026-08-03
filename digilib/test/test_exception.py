import pytest

from digilib.exception.uncatchable import Submission_paused, Signal_caught, Uncatchable_exception

@pytest.mark.parametrize("exception", [
    Submission_paused(None), Signal_caught(None, None)
], ids=["Submission_paused", "Signal_caught"])
def test_uncatachable(exception):
    """Can we only catch uncatachables if we really try?"""
    with pytest.raises(Uncatchable_exception):
        try:
            raise exception
        
        except Exception:
            # This should get ignored.
            pass
