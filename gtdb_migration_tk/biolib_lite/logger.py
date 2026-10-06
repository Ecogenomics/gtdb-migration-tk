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

import os
import sys
import logging

try:                                         # Python 2 fallback; kept in place
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists


# What a secret given on the command line is logged as.
REDACTED = '********'


def redacted_command_line(argv, secrets):
    """The command line as the log records it, with every secret masked.

    Masked by value rather than by flag: -p is the database password of the *_db
    commands and an output prefix or an rRNA path of others, so it is the value
    the parsed options call a password that is hidden, as its own argument or
    after '=' in --password=<value>.

    Parameters
    ----------
    argv : list of str
        The arguments, sys.argv[1:].
    secrets : iterable of str
        Values never to be written; empty and None ones are ignored.

    @return: the arguments joined by spaces, each secret replaced by REDACTED.
    """

    secrets = {secret for secret in secrets if secret}
    shown = []
    for arg in argv:
        if arg in secrets:
            arg = REDACTED
        elif '=' in arg and arg.split('=', 1)[1] in secrets:
            arg = arg.split('=', 1)[0] + '=' + REDACTED
        shown.append(arg)
    return ' '.join(shown)


def logger_setup(log_dir, log_file, program_name,software_name, version, silent, secrets=()):
    """Setup loggers.

    Two logger are setup which both print to the stdout and a
    log file when the log_dir is not None. The first logger is
    named 'timestamp' and provides a timestamp with each call,
    while the other is named 'no_timestamp' and does not prepend
    any information. The attribution 'is_silent' is also added
    to each logger to indicate if the silent flag is thrown.

    Parameters
    ----------
    log_dir : str
        Output directory for log file.
    log_file : str
        Desired name of log file.
    program_name : str
        Name of program.
    software_name : str
        Name of software.
    version : str
        Program version number.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    secrets : iterable of str
        Values of the command line never to be logged, e.g. a database password.
    """

    # setup general properties of loggers
    timestamp_logger = logging.getLogger('timestamp')
    timestamp_logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    no_timestamp_logger = logging.getLogger('no_timestamp')
    no_timestamp_logger.setLevel(logging.DEBUG)

    # the handlers of an earlier call go before this one's are added. Both loggers
    # are named, so they are the SAME objects on a second call, and the handlers
    # it added are still on them: every line of the run then reached the console
    # twice. That happens whenever the first call could not open the log it was
    # given -- __main__ falls back to ./gtdb_migration_tk.log and calls again --
    # and a command whose --log named a file where a directory was wanted, which
    # is one keystroke, printed itself double from beginning to end
    for logger in (timestamp_logger, no_timestamp_logger):
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()

    # setup logging to console
    timestamp_stream_logger = logging.StreamHandler(sys.stdout)
    timestamp_stream_logger.setFormatter(log_format)
    timestamp_logger.addHandler(timestamp_stream_logger)

    no_timestamp_stream_logger = logging.StreamHandler(sys.stdout)
    no_timestamp_stream_logger.setFormatter(None)
    no_timestamp_logger.addHandler(no_timestamp_stream_logger)

    timestamp_logger.is_silent = False
    no_timestamp_stream_logger.is_silent = False
    if silent:
        timestamp_logger.is_silent = True
        timestamp_stream_logger.setLevel(logging.ERROR)
        no_timestamp_stream_logger.is_silent = True

    if log_dir:
        make_sure_path_exists(log_dir)
        timestamp_file_logger = logging.FileHandler(os.path.join(log_dir, log_file), 'a')
        timestamp_file_logger.setFormatter(log_format)
        timestamp_logger.addHandler(timestamp_file_logger)

        no_timestamp_file_logger = logging.FileHandler(os.path.join(log_dir, log_file), 'a')
        no_timestamp_file_logger.setFormatter(None)
        no_timestamp_logger.addHandler(no_timestamp_file_logger)

    timestamp_logger.info('%s v%s' % (program_name, version))
    timestamp_logger.info(software_name + ' ' + redacted_command_line(sys.argv[1:], secrets))

def log_directory():
    """The directory the run's log is written to, for files that belong beside it.

    @return: the directory of the 'timestamp' logger's log file, or the current
             directory where the run is logged to the console alone.
    """

    for handler in logging.getLogger('timestamp').handlers:
        if isinstance(handler, logging.FileHandler):
            return os.path.dirname(handler.baseFilename)
    return os.getcwd()
