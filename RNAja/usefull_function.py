from pathlib import Path
import urllib.request
from tqdm import tqdm
import os
import click

INSTALL_PATH = Path(__file__).resolve().parent

def get_version():
    """Read VERSION file to know current version
    Returns:
        version: actual version read on the VERSION file
    Examples:
        version = get_version()
        print(version)
            1.3.0
    """
    with open(INSTALL_PATH.joinpath("VERSION"), 'r') as version_file:
        return version_file.readline().strip()