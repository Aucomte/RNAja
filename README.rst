|PythonVersions| |SnakemakeVersions| |Singularity|

.. contents:: Table of Contents
    :depth: 2

About RNAja
===============

RNAja is a pipeline written in snakemake.
It's aim is to analyse RNAseq data, perform differential expression analysis on it and output a report in Rmarkdown.
The output are compatible with the database `diffexDB <https://bioinfo-web.mpl.ird.fr/cgi-bin2/microarray/public/diffexdb.cgi>`_.

|readthedocs|

**Homepage:**

TODO


Authors
_______

* Aurore Comte (IRD)

Some parts of RNAja code and documentation were inspired or came from the pipelines below.

- Culebront (Julie Orjuela *et al.*) https://github.com/SouthGreenPlatform/culebrONT

- sRNAmake (Sebastien Cunnac *et al.*) https://github.com/Aucomte/sRNAmake

- BulkRNA (Camille Cohen) https://github.com/CamilleCohen/ProjetTuteur-_BulkRNA

Thanks and aknowledgements
==========================

Thanks to Ndomassi Tando (i-Trop IRD) by administration support.

The authors acknowledge the `IRD i-Trop HPC <https://bioinfo.ird.fr/>`_ (`South Green Platform <http://www.southgreen.fr>`_) at IRD
Montpellier for providing HPC resources that have contributed to this work.

Thanks to Alexis Dereeper for his help and the developpement of `diffexDB <https://bioinfo-web.mpl.ird.fr/cgi-bin2/microarray/public/diffexdb.cgi>`_.

License
=======

Licencied under `CeCill-C <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html>`_ and GPLv3.

Intellectual property belongs to IRD and authors.

.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
.. |readthedocs| image:: https://pbs.twimg.com/media/E5oBxcRXoAEBSp1.png
   :target: https://culebront-pipeline.readthedocs.io/en/latest/
   :width: 400px
