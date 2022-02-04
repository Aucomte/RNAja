from pathlib import Path
from .run import run_local
from .install import install_local
from .usefull_function import get_version
from .edit_files import create_config
from .global_variable import *

__version__=get_version()

__doc__ = """RNAja"""

description_tools = f"""RNAja"""