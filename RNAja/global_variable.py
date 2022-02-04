import sys
from pathlib import Path

#DOCS = "https://culebront-pipeline.readthedocs.io/en/latest/"
#GIT_URL = "https://github.com/SouthGreenPlatform/CulebrONT_pipeline"

# Hack for build docs with unspecified path install
args = str(sys.argv)
if "sphinx" in args:
    RNAJA_PATH = Path("/Path/to/culebrONT_install/")
else:
    RNAJA_PATH = Path(__file__).resolve().parent
    RNAJA_SNAKEFILE = RNAJA_PATH.joinpath("snakefiles", "Snakefile")
    RNAJA_SCRIPTS = RNAJA_PATH.joinpath("scripts")
    RNAJA_CONFIG_PATH = RNAJA_PATH.joinpath("install_files", "config.yaml")
    RNAJA_MODE = RNAJA_PATH.joinpath(".mode.txt")

SINGULARITY_URL_FILES = [('https://itrop.ird.fr/RNAja_utilities/Singularity.RNAja.sif', f'{RNAJA_PATH}/containers/Singularity.RNAja.sif')]