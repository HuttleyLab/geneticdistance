# thank you http://stackoverflow.com/questions/279237/import-a-module-from-a-relative-path

from os.path import realpath, abspath, dirname, join
from inspect import getfile, currentframe
import sys

parent = realpath(dirname(abspath(dirname(getfile(currentframe())))))
lib = join(parent, 'lib')

if lib not in sys.path:
    sys.path.insert(0, lib)
