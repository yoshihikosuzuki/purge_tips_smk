################################################################################
#                         Common requirements/imports
################################################################################

from snakemake.utils import min_version
min_version("6.0")
import os
configfile: "config.yaml"
shell.executable("/bin/bash")
shell.prefix(config["shell_prefix"])

################################################################################
#                              Utility functions
################################################################################

def parse_path(path):
    """Decompose /path/to/foo.txt into `("/path/to", "foo", "txt")`."""
    dirname, filename = os.path.split(path)
    basename, ext = os.path.splitext(filename)
    if ext.startswith('.'):
        ext = ext[1:]
    else:
        assert ext == ''
    return (dirname, basename, ext)

def get_modules(names):
    return ' '.join([config["modules"][name] for name in names])
