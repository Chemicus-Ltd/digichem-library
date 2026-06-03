"""Computational chemistry management"""

from datetime import datetime
import sys

import digichem.log
from digichem.datas import get_resource

####################
# Package metadata.#
####################

__version__ = "8.0.0-pre.3"
_v_parts = __version__.split("-")[0].split(".")
major_version = int(_v_parts[0])
minor_version = int(_v_parts[1])
revision = int(_v_parts[2])
if "-pre" in __version__:
    development = True
    prerelease = int(__version__.split(".")[3])
else:
    development = False
    prerelease = None


# Those who wrote this.
__author__ = [
    "Oliver S. Lee"
]

# Program date (when we were last updated). This is changed automatically.
_last_updated_string = "03/06/2026"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")

# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True

    # We are frozen, expand PATH to include the possibly bundled oprattle.
    from pathlib import Path
    import os
    freeze_dir = Path(sys._MEIPASS)
    if freeze_dir.name == "_internal":
        freeze_dir = freeze_dir.parent
    
    os.environ['PATH'] = os.environ.get("PATH", "") + ":" + str(freeze_dir)

else:
    frozen = False

import rdkit.RDLogger
# WrapLogs() outputs rdkit logging to python's stderr (which might be redirected to an urwid widget).
# If/when rdkit is further intergrated into digichem, this call will likely be moved elsewhere. 
#rdkit.Chem.rdchem.WrapLogs()
# Sadly the behaviour of WrapLogs() is a bit bizzare, although we do get redirection to our custom widgets etc,
# logs are also still dumped to screen...
# for now, disable logging...
rdkit.RDLogger.DisableLog('rdApp.*')

# Setup the logger
digichem.log.init_logger()

# At end to avoid circular imports.
import digichem.config
