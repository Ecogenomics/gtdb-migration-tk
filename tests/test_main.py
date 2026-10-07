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

"""Offline unit tests for main() in gtdb_migration_tk/__main__.py -- the exit
status of a run.

main() catches the SystemExit of a command that ends itself, and exited 0
whatever code the command gave it: a step that failed with sys.exit(-1) was
taken for one that succeeded by whatever ran it, a shell's && or a script. The
command is stood in for; what is tested is what main() does with how it ended.
"""

import argparse
import contextlib
import io
import logging
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module


class TheExitStatus(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='main_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        for name in ('timestamp', 'no_timestamp'):
            self.addCleanup(self.drop_handlers, logging.getLogger(name))

    @staticmethod
    def drop_handlers(logger):
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()

    def run_main(self, command_ends):
        """main() over a command that ends as command_ends does.

        @return: the code main() exited with, or None where it returned.
        """

        argv = ['gtdb_migration_tk', 'set_gtdb_domain', '--db_service', 'gtdb_r237',
                '-l', os.path.join(self.dir, 'run.log'), '--silent']
        with mock.patch.object(sys, 'argv', argv), \
                mock.patch.object(main_module.OptionsParser, 'parse_options',
                                  side_effect=command_ends), \
                contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            try:
                main_module.main()
            except SystemExit as exc:
                return exc.code
        return None

    def test_a_command_that_ends_itself_with_a_code_exits_with_that_code(self):
        self.assertEqual(self.run_main(SystemExit(-1)), -1)

    def test_a_command_that_ends_itself_with_a_message_exits_1(self):
        self.assertEqual(self.run_main(SystemExit('no genomes identified')), 1)

    def test_a_command_that_ends_itself_with_no_code_or_0_exits_0(self):
        self.assertIsNone(self.run_main(SystemExit()))
        self.assertIsNone(self.run_main(SystemExit(0)))

    def test_a_code_a_command_returns_is_still_the_exit_status(self):
        # ncbi_genome_sync returns its exit status rather than raising it
        self.assertEqual(self.run_main(lambda options: 3), 3)
        self.assertIsNone(self.run_main(lambda options: 0))



class TheCpusDefault(unittest.TestCase):
    """Every command's -c/--cpus defaults to 1."""

    @staticmethod
    def commands(parser, path=()):
        """(command path, parser) of every subcommand, nested ones included."""
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                for name, sub in action.choices.items():
                    yield path + (name,), sub
                    yield from TheCpusDefault.commands(sub, path + (name,))

    def test_every_command_taking_cpus_defaults_to_1(self):
        # create_tables and list_genomes defaulted to 8 and update_genomes to 16
        taking = {}
        for path, parser in self.commands(main_module.get_main_parser()):
            for action in parser._actions:
                if '--cpus' in action.option_strings:
                    taking[' '.join(path)] = action.default

        self.assertGreater(len(taking), 20)
        for command in ('create_tables', 'list_genomes', 'update_genomes'):
            self.assertIn(command, taking)
        self.assertEqual({command: default for command, default in taking.items() if default != 1}, {})


if __name__ == '__main__':
    unittest.main()
