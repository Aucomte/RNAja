.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

RNAja requires |PythonVersions|, |SnakemakeVersions| and |graphviz|.

------------------------------------------------------------------------

Steps for installation
======================

RNAja and dependencies and also every tool used to create a pipeline are available through a ``Singularity images``.

Install RNAja package
----------------------------

First install RNAja python package with pip.

.. code-block:: bash

   python3 -m pip install RNAja
   RNAja --help

Then run the commande line to install on LOCAL

.. code-block:: bash

   RNAja install_local --help
   RNAja install_local

The script automatically download singularity images required and configure RNAja.

Test RNAja install (Optional but recommended) using an available dataset.
See the section :ref:`Check install` for details.

------------------------------------------------------------------------

Check install
=============

It's is time to prepare the configuration file ``config.yaml`` file to indicate the paths to your data and the parameters.

Finally run RNAja!!

.. code-block:: bash

    RNAja run_local --config config.yaml --threads 6 --dry-run

------------------------------------------------------------------------

Advance installation
====================

How to build singularity images
-------------------------------

You can build your own image using the available *.def* recipe from the ``RNAja/RNAja/envs/`` directory.

.. warning::
    Be careful, you need root rights to build singularity images

.. code-block:: bash

    cd RNAja/RNAja/envs/
    sudo singularity build Singularity.RNAja.sif Singularity.RNAja.def

------------------------------------------------------------------------


.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
   :alt: Python 3.7+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+

.. |graphviz| image:: https://img.shields.io/badge/graphviz-%3E%3D2.40.1-green
   :target: https://graphviz.org/
   :alt: graphviz 2.40.1+
