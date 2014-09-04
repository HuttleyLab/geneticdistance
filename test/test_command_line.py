try:
    from mpi4py import MPI
    USING_MPI = MPI.COMM_WORLD.Get_size() > 1
except ImportError:
    USING_MPI = False

from nose.tools import assert_equal, assert_in
from os.path import realpath, abspath, dirname, join
from shutil import rmtree
from inspect import getfile, currentframe
from tempfile import gettempdir
import sys

import src
from data import get_data_dir

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

class TestCommandLine(object):
    def setup(self):
        self._output = join(gettempdir(), 'test_cmd_line_output')
        self._argv = sys.argv

    def teardown(self):
        try:
            rmtree(self._output)
        except OSError:
            pass
        sys.argv = self._argv

    def test_command_line(self):
        datadir = get_data_dir()
        logfile = join(self._output, 'nsl.log')
        cmd = 'nonstationary_lengths.py -i ' + datadir + ' -o ' + \
                self._output + ' -l ' + logfile + ' -L DEBUG -c 3 -u 20 -F seq_fit'
        sys.argv = cmd.split()
        from nonstationary_lengths import main
        assert_equal(main(), 0)
        if USING_MPI:
            MPI.COMM_WORLD.Barrier()

        log = ''
        with open(logfile) as logfile:
            log = logfile.read()
        assert_in('Done Mouse/Opossum/Human in ENSG00000111145', log)

        logfile = join(self._output, 'gs.log')
        cmd = 'g_stats.py -o ' + self._output + ' -l ' + logfile + \
                ' -N 1 -u 20 -P 1 -F seq_fit -L DEBUG'
        sys.argv = cmd.split()
        from g_stats import main
        assert_equal(main(), 0)
        if USING_MPI:
            MPI.COMM_WORLD.Barrier()

        log = ''
        with open(logfile) as logfile:
            log = logfile.read()
        assert_in('Done General and Mouse/Opossum/Human in ENSG00000111145', log)

def main():
    tester = TestCommandLine()
    tester.setup()
    try:
        tester.test_command_line()
    finally:
        tester.teardown()
    return 0

if __name__ == '__main__':
    sys.exit(main())
