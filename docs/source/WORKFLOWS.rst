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

Then edit section on file according to your workflow.

First, indicate the data path in the configuration ``config.yaml`` file:

.. code-block:: YAML

    'DATA':
        'directories':
            'out_dir': "output2"
        'files':
            'reference': "DATA/ref/allcon.fasta"
            'annotation': "DATA/ref/msu7.gtf"
            'sample_info': "DATA/sample_info.txt"
            'de_comparisons_file': "DATA/treatmentsComparisons.csv"

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
