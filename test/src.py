# thank you http://stackoverflow.com/questions/279237/import-a-module-from-a-relative-path

from os.path import realpath, abspath, dirname, join
from inspect import getfile, currentframe
import sys

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

parent = realpath(dirname(abspath(dirname(getfile(currentframe())))))
src = join(parent, 'src')

if src not in sys.path:
    sys.path.insert(0, src)
