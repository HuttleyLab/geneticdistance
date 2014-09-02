from __future__ import division

import lib
import masterslave

from nose.tools import assert_equal, assert_less, assert_in
from tempfile import gettempdir
from socket import gethostname
import logging
import types

from shutil import rmtree
import os
import sys

def test_farm():
    size = masterslave.size()

    def func(dummy):
        return dummy, masterslave.rank()

    def test_n(n):
        array = range(n)
        result = masterslave.farm(func, array)
        if result:
            result = list(result)
            if result:
                array, ranks = zip(*result)
                assert_equal(len(array), n)
                assert_equal(set(array), set(range(n)))
                proper_ranks = set(range(min(size-1, 1), min(n+1, size)))
                assert_equal(set(ranks), proper_ranks)

    for n in [size//2, size-1, size, size*2]:
        test_n(n)

def test_farm_gav():
    size = masterslave.size()

    def func(d1, d2):
        return d1, d2, masterslave.rank()

    def test_n(n):
        array = range(n)
        result = masterslave.farm(func, array, array)
        if result:
            result = list(result)
            if result:
                array, array1, ranks = zip(*result)
                assert_equal(array, array1)
                assert_equal(len(array), n)
                assert_equal(set(array), set(range(n)))
                proper_ranks = set(range(min(size-1, 1), min(n+1, size)))
                assert_equal(set(ranks), proper_ranks)

    for n in [size//2, size-1, size, size*2]:
        test_n(n)

def test_map():
    size = masterslave.size()

    def func(dummy):
        return dummy, masterslave.rank()

    def test_n(n):
        array = range(n)
        result = masterslave.map(func, array)
        if result:
            result = list(result)
            if result:
                output, ranks = zip(*result)
                assert_equal(list(output), array)
                proper_ranks = set(range(min(size-1, 1), min(n+1, size)))
                assert_equal(set(ranks), proper_ranks)
    
    for n in [size//2, size-1, size, size*2]:
        test_n(n)

def test_map_gav():
    size = masterslave.size()

    def func(d1, d2):
        return d1, d2, masterslave.rank()

    def test_n(n):
        array = range(n)
        result = masterslave.map(func, array, array)
        if result:
            result = list(result)
            if result:
                output, output1, ranks = zip(*result)
                assert_equal(output, output1)
                assert_equal(list(output), array)
                proper_ranks = set(range(min(size-1, 1), min(n+1, size)))
                assert_equal(set(ranks), proper_ranks)
    
    for n in [size//2, size-1, size, size*2]:
        test_n(n)

def test_imap():
    size = masterslave.size()

    def func(dummy):
        return dummy, masterslave.rank()

    def test_n(n):
        array = range(n)
        result = masterslave.imap(func, array)
        if size != 1:
            assert_equal(type(result), types.GeneratorType)
        result = list(result)
        if result:
            output, ranks = zip(*result)
            assert_equal(list(output), array)
            proper_ranks = set(range(min(size-1, 1), min(n+1, size)))
            assert_equal(set(ranks), proper_ranks)

    for n in [size//2, size-1, size, size*2]:
        test_n(n)

class TestIO(object):
    def setup(self):
        self.tempdir = os.path.join(gettempdir(), 'test_io')
        masterslave.checkmakedirs(self.tempdir)
        self.tempfilename = os.path.join(self.tempdir, 'tempfile')

    def teardown(self):
        try:
            rmtree(self.tempdir)
        except OSError:
            pass

    def test_io(self):
        tempfilename = self.tempfilename
        tempfile = masterslave.open(tempfilename, mode='w')
        tempfile.write(str(masterslave.rank())+'\n')
        tempfile.flush()
        contents = ''
        with open(tempfilename) as stillopen:
            contents = stillopen.read()
        assert_in(str(masterslave.rank())+'\n', contents)
        tempfile.close()
        with open(tempfilename) as nowclosed:
            ranks = [int(s.strip()) for s in nowclosed]
        assert_equal(set(ranks), set(range(masterslave.size())))
        tempfile = masterslave.open(tempfilename, 'a')
        tempfile.write(str(masterslave.rank())+'\n')
        tempfile.close()
        with open(tempfilename) as nowclosed:
            ranks = [int(s.strip()) for s in nowclosed]
        assert_equal(len(ranks), masterslave.size()*2)
        assert_equal(set(ranks), set(range(masterslave.size())))
        if masterslave.size() != 1:
            masterslave._comm.Barrier()

class TestMPIFileHandler(object):
    def setup(self):
        self.tempdir = os.path.join(gettempdir(), 'test_MPIFileHandler')
        masterslave.checkmakedirs(self.tempdir)
        self.tempfilename = os.path.join(self.tempdir, 'tempfile')

    def teardown(self):
        try:
            rmtree(self.tempdir)
        except OSError:
            pass

    def test_MPIFileHandler(self):
        tempfilename = self.tempfilename
        handler = masterslave.MPIFileHandler(tempfilename, 'w')
        host = gethostname()
        formatter = logging.Formatter(host+':%(message)s')
        handler.setFormatter(formatter)
        logging.root.addHandler(handler)
        logging.warn('TEST')
        logging.getLogger().removeHandler(handler)
        handler.flush()
        handler.close()
        lines = {}
        with open(tempfilename) as tempfile:
            lines = set(tempfile.readlines())
        assert_equal(lines, set([host+':TEST\n']))
        if masterslave.size() != 1:
            masterslave._comm.Barrier()

def main():
    test_farm()
    test_map()
    test_imap()
    tester = TestIO()
    tester.setup()
    try:
        tester.test_io()
    finally:
        tester.teardown()
    tester = TestMPIFileHandler()
    tester.setup()
    try:
        tester.test_MPIFileHandler()
    finally:
        tester.teardown()

if __name__ == '__main__':
    sys.exit(main())
