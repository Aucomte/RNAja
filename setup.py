#python3.8 -m pip install -e .
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
from setuptools import setup, find_packages

NAME = 'RNAja'
URL = 'https://github.com/Aucomte/RNAja'
CURRENT_PATH = Path(__file__).resolve().parent
VERSION = CURRENT_PATH.joinpath(f'RNAja','VERSION').open('r').readline().strip()

def main():
    setup(
        # Project information
        name=NAME,
        version=VERSION,
        url=URL,
        project_urls={
            'Bug Tracker': f'{URL}/issues',
            #'Documentation': f'https://{NAME}.readthedocs.io/en/latest/',
            'Source Code': URL
        },
        download_url=f'{URL}/archive/{VERSION}.tar.gz',
        author='''Aurore Comte (IRD)''',
        author_email='aurore.comte@ird.fr',
        description=CURRENT_PATH.joinpath('README.rst').open('r', encoding='utf-8').read(),
        long_description=CURRENT_PATH.joinpath('README.rst').open('r', encoding='utf-8').read(),
        long_description_content_type='text/x-rst',
        license='GPLv3',

        # docs compilation utils
        #command_options={
        #    'build_sphinx': {
        #        'project': ('setup.py', NAME),
        #        'version': ('setup.py', VERSION),
        #        'release': ('setup.py', VERSION),
        #        'source_dir': ('setup.py', CURRENT_PATH.joinpath('docs','source').as_posix()),
        #        'build_dir': ('setup.py', CURRENT_PATH.joinpath('docs','build').as_posix()),
        #    }},

        # Package information
        packages=find_packages(),
        package_data={
            '': ['*'],
        },
        include_package_data=True,
        python_requires='>=3.6',
        install_requires=[
            'snakemake',
            'docutils<0.18',
            'click>=8.0.3'
        ],
        extras_require={
            'docs': ['sphinx_copybutton',
                     'sphinx_rtd_theme',
                     'sphinx_click'],
            'dev': ['tox'],
        },
        entry_points={
            'RNAja': ['RNAja = __init__'],
            'console_scripts': [f'RNAja = RNAja.__init__:run_local']
            },

        # Pypi information
        platforms=['unix', 'linux'],
        keywords=[
            'snakemake',
            'RNAseq',
            'workflow',
            'edgeR'
        ],
        classifiers=[
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'License :: CeCILL-C Free Software License Agreement (CECILL-C)',
            'License :: Free for non-commercial use',
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Intended Audience :: End Users/Desktop',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: R',
            'Natural Language :: English',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
        options={
            'bdist_wheel': {'universal': False}
        },
        zip_safe=False,  # Don't install the lib as an .egg zipfile
    )


if __name__ == '__main__':
    main()
