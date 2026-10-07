from pathlib import Path

import digilib
from digilib.datas import get_resource


# A config file included in the digichem install directory.
master_config_path = get_resource("data/config/digichem.yaml")

# The location of the system-wide config file (which has precedence over the master).
system_config_location = Path("/etc/digichem{}/digichem.yaml".format(digilib.major_version))

# The location of the user specific config file (which has precedence over the system).
user_config_location = Path(Path.home(), ".config/digichem{}/digichem.yaml".format(digilib.major_version))