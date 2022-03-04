.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Run RNAja
========================

You must provide data paths, options and singulartity path in the ``config.yaml`` file

To create file juste run

.. code-block:: bash

   RNAja create_config --help
   RNAja create_config -configyaml Path/to/file

Then edit section on file according to your data.

First, indicate the data path in the configuration ``config.yaml`` file:

.. code-block:: YAML

    'DATA':
        'directories':
            'out_dir': "output"
        'files':
            'reference': "DATA/ref/allcon.fasta"
            'annotation': "DATA/ref/msu7.gtf"
            'sample_info': "DATA/sample_info.txt"


.. csv-table::

        "Input", "Description"
        "out_dir", "path / name of the output directory"
        "reference", "Only one REFERENCE genome file (fasta file)"
        "annotation", "Annotation file corresponding to the reference file (gtf or gff)"
        "sample_info", "path of the file containing information on the sample and the experimental design (csv or txt, separator = ',')"


exemple of sample_info file:

.. csv-table::

        "FileName", "SampleName", "Condition1", "Condition2"
        "/path/to/fastq/Control1", "C_1", "Control", "E1"
        "/path/to/fastq/Control2", "C_2", "Control", "E2"
        "/path/to/fastq/Control3", "C_3", "Control", "E3"
        "/path/to/fastq/Treatment1", "T_1", "Treat", "E1"
        "/path/to/fastq/Treatment2", "T_2", "Treat", "E2"
        "/path/to/fastq/Treatment3", "T_3", "Treat", "E3"


.. warning::

    For the column SampleName of the sampleInfo file, the names should be: condition_repetition separated by a _

.. ############################################################

How to run the workflow
=======================

Before run RNAja please be sure you have already modified the ``config.yaml`` file.
You can optionally also pass to snakemake more options by using non option parameter (check it https://snakemake.readthedocs.io/en/stable/executing/cli.html).

.. code-block:: bash

    # in LOCAL using maximum 8 threads
    RNAja run_local -c config.yaml --cores 8 --dry-run

    # in LOCAL using 6 threads for hisat2 from the total 8 threads
    RNAja run_local -c config.yaml ---cores 8 --set-threads hisat2_map=6

Output on RNAja
===================

The architecture of RNAja output is designed as follows:

.. code-block:: bash

    OUTPUT_RNAJA/
    ├── 1_QC
    │   ├── fastqc
    ├── 2_mapping
    ├── 3_count
    │   ├── STRINGTIE
    ├── 4_DE_analysis
    ├── 5_MULTIQC
    └── LOGS
