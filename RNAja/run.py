import click
from shutil import copyfile
from RNAja.global_variable import *
import os


@click.command("run_local", short_help='Run workflow on local computer (use singularity mandatory)', context_settings=dict(ignore_unknown_options=True, max_content_width=800),
               no_args_is_help=True)
@click.option('--config', '-c', type=click.Path(exists=True, file_okay=True, readable=True, resolve_path=True), required=True, help='Configuration file for run culebrONT')
@click.option('--threads', '-t', type=int, required=True, help='number of threads')
@click.option('--pdf', '-p', is_flag=True, required=False, default=False, help='run snakemake with --dag, --rulegraph and --filegraph')
@click.argument('snakemake_other', nargs=-1, type=click.UNPROCESSED)
def run_local(config, threads, pdf, snakemake_other):
    """    Run snakemake command line with mandatory parameters.
    SNAKEMAKE_OTHER: append other snakemake command such '--dry-run'
    Example:
        culebrONT run_local -c config.yaml --dry-run
    """
    # get user arguments
    click.secho(f'    Config file: {config}', fg='yellow')

    cmd_snakemake_base = f"snakemake --latency-wait 60000 --cores {threads} --use-singularity --singularity-args \'--bind $HOME\' --show-failed-logs --printshellcmds -s {RNAJA_SNAKEFILE} --configfile {config}  {' '.join(snakemake_other)}"
    click.secho(f"\n    {cmd_snakemake_base}\n", fg='bright_blue')

    os.system(cmd_snakemake_base)
    if pdf:
        dag_cmd_snakemake = f"{cmd_snakemake_base} --dag | dot -Tpdf > schema_pipeline_dag.pdf"
        click.secho(f"    {dag_cmd_snakemake}\n", fg='bright_blue')
        os.system(dag_cmd_snakemake)
        rulegraph_cmd_snakemake = f"{cmd_snakemake_base} --rulegraph | dot -Tpdf > schema_pipeline_global.pdf"
        click.secho(f"    {rulegraph_cmd_snakemake}\n", fg='bright_blue')
        os.system(rulegraph_cmd_snakemake)
        filegraph_cmd_snakemake = f"{cmd_snakemake_base} --filegraph | dot -Tpdf > schema_pipeline_files.pdf"
        click.secho(f"    {filegraph_cmd_snakemake}\n", fg='bright_blue')
        os.system(filegraph_cmd_snakemake)