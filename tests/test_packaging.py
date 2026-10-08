###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""Offline unit tests for what an install of the package holds.

The conda package of a release is built by 'pip install .' (conda/meta.yaml),
which installs what setup.py says. package_data named VERSION alone, so the
files the package reads from gtdb_migration_tk/data_files -- the descriptions
update_metadata_db loads tables against, the barrnap models rna_silva searches
with -- were in no install but one made against a checkout, where they are
read from the tree. A wheel is built from a copy of the package, as pip builds
one, without the network, and its files are compared with the tree's.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PACKAGE = 'gtdb_migration_tk'


class WhatAnInstallHolds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(prefix='packaging_test.')
        source = os.path.join(cls.dir, 'source')
        os.makedirs(source)
        shutil.copy(os.path.join(REPO, 'setup.py'), source)
        shutil.copytree(os.path.join(REPO, 'bin'), os.path.join(source, 'bin'))
        shutil.copytree(os.path.join(REPO, PACKAGE), os.path.join(source, PACKAGE),
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        wheels = os.path.join(cls.dir, 'wheels')
        built = subprocess.run([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '--no-build-isolation',
                                '--no-index', '-q', '-w', wheels, source],
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        if built.returncode != 0:
            shutil.rmtree(cls.dir, True)
            raise unittest.SkipTest('no wheel could be built here: ' + built.stdout[-500:])
        wheel = [name for name in os.listdir(wheels) if name.endswith('.whl')][0]
        cls.names = set(zipfile.ZipFile(os.path.join(wheels, wheel)).namelist())

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, True)

    def test_every_data_file_of_the_tree_is_installed(self):
        data = os.path.join(REPO, PACKAGE, 'data_files')
        tree = sorted(os.path.relpath(os.path.join(root, name), REPO).replace(os.sep, '/')
                      for root, _dirs, names in os.walk(data) for name in names)

        self.assertGreater(len(tree), 20)
        self.assertEqual([name for name in tree if name not in self.names], [])

    def test_the_version_is_installed(self):
        self.assertIn(PACKAGE + '/VERSION', self.names)


if __name__ == '__main__':
    unittest.main()
