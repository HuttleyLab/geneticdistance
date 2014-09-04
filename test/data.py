from gzip import GzipFile
from os.path import realpath, abspath, dirname, join
from inspect import getfile, currentframe

from cogent import Alignment

_data_dir = join(realpath(abspath(dirname(getfile(currentframe())))), 'data')

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def get_data_dir():
    return _data_dir

def get_aln(model, aln_len):
    filename = '_'.join((model, str(aln_len))) + '.fasta.gz'
    data = ''
    with GzipFile(join(_data_dir, filename)) as fastafile:
        data = fastafile.read()
    return Alignment(data=data)
