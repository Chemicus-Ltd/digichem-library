import pytest
from pathlib import Path
import yaml
import copy

from digichem.test.util import digichem_options, data_directory
from digichem.config import Auto_type, get_config
from digichem.config.parse import Config_parser

def test_config_save(digichem_options, tmpdir):
    """Can we save config options to file?"""
    # First, check an exception is thrown if we try and save to a non-existent dir:
    # with pytest.raises(Exception):
    #     digichem_options.save(Path(tmpdir, "foobar", "conf.yaml"))

    # Now save properly.
    save_file = Path(tmpdir, "conf.yaml")
    digichem_options.save(save_file)

    # Open the file and read it back.
    with open(save_file) as file:
        data = yaml.safe_load(file)

    # Check the data is the same.
    assert digichem_options.dump() == data

@pytest.mark.parametrize("var, typeis", [
    ["1", int],
    ["1.1", float],
    ["one", str]
], ids=["int", "float", "str"])
def test_autotype(var, typeis):
    """Does autotype work?"""
    assert type(Auto_type(var)) is typeis

def test_string_parser():
    """Can we parse config options from a string."""
    Config_parser("foo: bar").load() == {"foo": "bar"}

def test_bad_string_parser():
    """Can we not parse config options from a malformed string."""
    with pytest.raises(Exception):
        Config_parser("foo: {bar").load()

def test_del_option(digichem_options):
    """Can we reset an option to its default by deleting it?"""
    options = copy.deepcopy(digichem_options)
    # First, set a value.
    options.alignment = "FAP"

    # Delete and check.
    del options.alignment
    assert options.alignment == "NEST"

def test_config_merge():
    """Can we set additional options from a string?"""
    conf = get_config(
        clear_cache = True,
        extra_config_strings = ["alignment: FAP"]
    )

    assert conf.alignment == "FAP"
    get_config(clear_cache = True)

def test_config_file_merge(tmpdir):
    """Can we set additional options from another file?"""
    with open(Path(tmpdir, "conf.yaml"), "w") as conf_file:
        yaml.dump({'alignment': 'FAP'}, conf_file)

    conf = get_config(
        clear_cache = True,
        extra_config_files = [Path(tmpdir, "conf.yaml")]
    )

    assert conf.alignment == "FAP"

    # Reset the config object.
    get_config(clear_cache = True)