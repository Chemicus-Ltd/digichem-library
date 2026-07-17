import atexit
from contextlib import ExitStack

try:
    import importlib.resources
    importlib_resources = importlib.resources

except ImportError:
    import importlib_resources

def get_resource(name, pkg = 'digichem'):
    """
    Get a pathlib path object to a package resource.
    """
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    ref = importlib_resources.files(pkg) / name
    return file_manager.enter_context(importlib_resources.as_file(ref))