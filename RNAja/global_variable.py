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

#SINGULARITY_URL_FILES = [('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.culebront_tools.sif',
#              f'{CULEBRONT_PATH}/containers/Singularity.culebront_tools.sif'),
#             ('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.report.sif',
#              f'{CULEBRONT_PATH}/containers/Singularity.report.sif')]