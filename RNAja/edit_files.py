import click
from shutil import copyfile
from RNAja.global_variable import *
from RNAja.usefull_function import get_install_mode

@click.command("create_config", short_help='Create config.yaml for run', no_args_is_help=True)
@click.option('--configyaml', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create config.yaml')
def create_config(configyaml):
    configyaml = Path(configyaml)
    configyaml.parent.mkdir(parents=True, exist_ok=True)
    copyfile(RNAJA_CONFIG_PATH.as_posix(), configyaml.as_posix())
    click.edit(require_save=True, extension='.yaml', filename=configyaml)
    click.secho(f"\n    Success create config file on path '{configyaml}'\n    add to command:", fg="yellow")
    mode = get_install_mode()
    click.secho(f"    RNAja {'run_cluster' if mode == 'cluster' else 'run_local'} --config {configyaml}\n\n", fg="bright_blue")

