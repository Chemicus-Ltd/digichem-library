import numpy
import os
from pathlib import Path

# Set numpy errors (not sure why this isn't the default...)
numpy.seterr(invalid = 'raise', divide = 'raise')

# Expand path to include mocks.
os.environ["PATH"] = str(Path(__file__).parent / "mock") + os.pathsep + os.environ["PATH"]