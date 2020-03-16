#!/usr/bin/env python

import os

from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'gtdb_migration_tk', 'VERSION'), 'r') as f:
        return f.readline().strip()


setup(
    name='gtdb_migration_tk',
    python_requires='>=3.6',
    version=version(),
    author='Pierre-Alain Chaumeil',
    author_email='p.chaumeil@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil, Aaron Mussig, and Donovan Parks',
    maintainer_email='p.chaumeil@uq.edu.au',
    packages=['gtdb_migration_tk'],
    scripts=['bin/gtdb_migration_tk'],
    package_data={'gtdb_migration_tk': ['VERSION']},
    url='https://github.com/Ecogenomics/gtdb-migration-tk',
    description='Toolkit for updating the GTDB to the next release and test data.',
    install_requires=['requests', 'unidecode', 'biolib>=0.1.0', 'pandas', 'numpy'],
)
