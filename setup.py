#!/usr/bin/env python

import os

from setuptools import setup, find_packages


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'gtdb_migration_tk', 'VERSION'), 'r') as f:
        return f.readline().strip()


setup(
    name='gtdb_migration_tk',
    python_requires='>=3.8',
    version=version(),
    author='Pierre-Alain Chaumeil',
    author_email='p.chaumeil@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil, Aaron Mussig, and Donovan Parks',
    maintainer_email='p.chaumeil@uq.edu.au',
    packages=find_packages(include=['gtdb_migration_tk', 'gtdb_migration_tk.*']),
    scripts=['bin/gtdb_migration_tk'],
    package_data={'gtdb_migration_tk': ['VERSION']},
    url='https://github.com/Ecogenomics/gtdb-migration-tk',
    description='Toolkit for updating the GTDB to the next release and test data.',
    install_requires=[
        'requests>=2.27.1',
        'unidecode>=1.3.4',
        'pandas>=1.4.1',
        'numpy>=1.22.3',
        'sqlalchemy>=1.4.45',
        'beautifulsoup4>=4.11.1',
        'dendropy>=4.5.2',
        'tqdm>=4.63.0',
        'atpbar>=1.1.4',
        'python-dateutil>=2.8.2',
        'psycopg2-binary>=2.9.3',
    ],
)
