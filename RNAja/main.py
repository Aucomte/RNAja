#!/usr/bin/env python3
import os
import click
import RNAja
from RNAja.global_variable import *
from RNAja.usefull_function import get_install_mode, check_privileges


@click.group(help=click.secho(RNAja.description_tools, fg='green', nl=False), context_settings={'help_option_names': ('-h', '--help'),"max_content_width":800},
             invoke_without_command=True, no_args_is_help=True)
@click.option('--restore', '-r', is_flag=True, required=False, default=False, show_default=True, help='Restore installation mode (need root or sudo)')
@click.version_option(RNAja.__version__, '--version', '-v')
@click.pass_context
def main(ctx, restore):
    if ctx.invoked_subcommand is None and restore and check_privileges():
        if RNAJA_MODE.exists():
            RNAJA_MODE.unlink(missing_ok=False)
            click.secho(f"\n    Remove installation mode, now run:\n    RNAja install_local or install_cluster\n\n", fg="yellow")
        else:
            click.secho(f"\n    No reset need, RNAja not install !!!!!\n    Please run: RNAja install_local or install_cluster !!!!\n\n", fg="red")
    pass


# Hack for build docs with unspecified install
args = str(sys.argv)
if "sphinx" in args:
    #main.add_command(RNAja.run_cluster)
    #main.add_command(RNAja.create_cluster_config)
    main.add_command(RNAja.create_config)
    #main.add_command(RNAja.edit_tools)
    main.add_command(RNAja.run_local)
    #main.add_command(RNAja.install_cluster)
    main.add_command(RNAja.install_local)
    #main.add_command(RNAja.test_install)
else:
    mode = get_install_mode()
    if mode == "cluster":
        #main.add_command(RNAja.test_install)
        #main.add_command(RNAja.run_cluster)
        #main.add_command(RNAja.create_cluster_config)
        main.add_command(RNAja.create_config)
        #main.add_command(RNAja.edit_tools)
    elif mode == "local":
        #main.add_command(RNAja.test_install)
        main.add_command(RNAja.run_local)
        main.add_command(RNAja.create_config)
    else:
        #main.add_command(RNAja.install_cluster)
        main.add_command(RNAja.install_local)

if __name__ == '__main__':
    main()
