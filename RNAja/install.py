import click
from click import Abort
from pathlib import Path
from shutil import rmtree, copyfile, unpack_archive
import RNAja
import os
import re
from cookiecutter.main import cookiecutter
from RNAja.global_variable import *
from RNAja.usefull_function import command_required_option_from_option, multiprocessing_download, get_install_mode

required_options = {
    True: 'modules_dir',
    False: 'scheduler'
}

def create_bash_completion():
    bashrc_file = Path("~/.bashrc").expanduser().as_posix()
    import subprocess
    output = subprocess.run(
        ["bash", "-c", "echo ${BASH_VERSION}"], stdout=subprocess.PIPE
    )
    match = re.search(r"^(\d+)\.(\d+)\.\d+", output.stdout.decode())
    if match is not None:
        major, minor = match.groups()
    if major < "4" or (major == "4" and minor < "4"):
        click.secho(f"\n    WARNNING Shell completion is not supported for Bash versions older than 4.4.", fg="red", nl=False)
        # raise RuntimeError("Shell completion is not supported for Bash versions older than 4.4.")
    else:
        if not Path(f"{RNAJA_PATH}/RNAja-complete.sh").exists():
            build_completion = f"_CULEBRONT_COMPLETE=bash_source RNAja > {RNAJA_PATH}/RNAja-complete.sh"
            os.system(build_completion)
        with open(bashrc_file, "r") as bash_file_read:
            if not [True for line in bash_file_read if "RNAJA" in line]:
                with open(bashrc_file, "a") as bash_file_open:
                    append_bashrc = f"\n#Add autocompletion for RNAJA\n. {RNAJA_PATH}/RNAja-complete.sh"
                    bash_file_open.write(append_bashrc)
                    click.secho(f"\n    INSTALL autocompletion for RNAJA on {bashrc_file} with command {append_bashrc}\nUpdate with commande:\nsource ~/.bashrc",
                                fg="yellow")
            else:
                path_culebront = ""
                with open(bashrc_file, "r") as bash_file_read:
                    for line in bash_file_read:
                        if "RNAJA" in line:
                            path_bash = bash_file_read.readline().strip()
                            path_culebront = f"{line}{path_bash}"
                load = f"{path_bash[2:]}"
                if f"{RNAJA_PATH}/RNAja-complete.sh" != load:
                    click.secho(
                        f"\n    WARNNING autocompletion for RNAJA  already found on {bashrc_file}, with other path please fix the good:",
                        fg="red", nl=False)
                    click.secho( f"\n    Load on bashrc: {load}\n    New install:    {RNAJA_PATH}/RNAja-complete.sh", fg='bright_red')


def create_envmodules(modules_dir):
    from RNAja import MODULE_FILE
    modules_dir = Path(modules_dir)
    modules_dir.mkdir(parents=True, exist_ok=True)
    modules_dir.joinpath(f"{RNAja.__version__}").open("w").write(MODULE_FILE)
    click.edit(require_save=False, extension='', filename=modules_dir.joinpath(f"{RNAja.__version__}").as_posix())
    click.secho(f"\n    Success install module file for version {RNAja.__version__} on path '{modules_dir}'", fg="yellow")

def clean_home():
    if Path("~/.config/RNAja/").expanduser().exists():
        rmtree(Path("~/.config/RNAja/").expanduser().as_posix())

def check_and_download_singularity():
    # check if already download if true pop from list
    SINGULARITY_URL_FILES_DOWNLOAD = []
    for imgs_list in SINGULARITY_URL_FILES:
        url, path_install = imgs_list
        if not Path(path_install).exists():
            SINGULARITY_URL_FILES_DOWNLOAD.append(imgs_list)
        else:
            click.secho(f"    File: {path_install} already download, done.", fg="yellow", nl=True)
    results = multiprocessing_download(SINGULARITY_URL_FILES_DOWNLOAD)
    for r in results:
        click.secho(r, fg="blue")
    click.secho(f"\n    WARNNING please check if binding is active on your singularity configuration, see https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html", fg = "bright_red")

@click.command("install_local", short_help='Install RNAja on local computer', context_settings=dict(max_content_width=800), no_args_is_help=False)
@click.option('--bash_completion/--no-bash_completion', is_flag=True, required=False, default=True, show_default=True,
              help='Allow bash completion of RNAja commands on the bashrc file')
def install_local(bash_completion):
    """
    \b
    Run installation for local computer with download the singularity images.
    """
    # add file for installation mode:
    mode_file = RNAJA_MODE
    # rm previous install (ie @julie cluster then local)
    clean_home()
    # add path to download
    RNAJA_PATH.joinpath("envs").mkdir(exist_ok=True, parents=True)
    try:
        check_and_download_singularity()
        # export to add bash completion
        if bash_completion:
           create_bash_completion()
        click.secho(f"\n    Congratulations, you have successfully installed RNAja !!!\n\n", fg="green", bold=True)
        mode_file.open("w").write("local")
    except Exception as e:
        mode_file.unlink(missing_ok=True)
        click.secho(f"\n    ERROR : an error was detected, please check {e}", fg="red")
        exit()
