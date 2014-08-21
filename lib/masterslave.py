try:
    from mpi4py import MPI
    USING_MPI = MPI.COMM_WORLD.Get_size() > 1
except ImportError:
    USING_MPI = False

from itertools import izip
from logging import StreamHandler
import sys
import os

__version__ = '0.0.9-dev'

_MASTER = 0
_WORKTAG = 1
_DIETAG = 2

def am_master():
    return _rank == _MASTER

def rank():
    return _rank

def size():
    return _size

def checkmakedirs(dirname):
    if _rank == _MASTER and not os.path.exists(dirname):
        os.makedirs(dirname)
    if USING_MPI:
        _comm.Barrier()

def exit(code):
    if USING_MPI:
        _comm.Abort(code)
    sys.exit(code)

def _open(filename, mode='a'):
    """Create an MPI safe file object for use with logging.StreamHandler"""
    return _file(filename, mode=mode)

class _file(object):
    """Only intended to fulfill the requirements of logging.StreamHandler"""
    def __init__(self, filename, mode='a'):
        if mode not in 'aw':
            raise ValueError(mode+' is not a valid mode for this log file')
        filemode = MPI.MODE_CREATE | MPI.MODE_WRONLY
        if 'a' in mode:
            filemode |= MPI.MODE_APPEND
        self._am_open = True
        self._stream = MPI.File.Open(_comm, filename, filemode)
        if 'a' not in mode: # MPI.File.Open never truncates
            self._stream.Set_size(0)
            _comm.Barrier()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def write(self, b):
        self._stream.Write_shared(b.encode('utf-8')) # unicode crashes mpi4py

    def flush(self):
        self._stream.Sync()

    def close(self):
        _comm.Barrier() # make sure everyone agrees
        self._stream.Close()

def _farm(function, *args):
    if _rank == _MASTER:
        return _master(i for i in izip(*args))
    else:
        _slave(lambda arg: function(*arg))

def _imap(function, *args):
    def tracker(envelope):
        order, arg = envelope
        return order, function(*arg)
    if _rank == _MASTER:
        cache = {}
        upto = 0
        for order, result in _master(enumerate(izip(*args))):
            cache[order] = result
            while cache and min(cache) == upto:
                yield cache[upto]
                del cache[upto]
                upto += 1
    else:
        _slave(tracker)
            
def _map(function, *args):
    return list(_imap(function, *args))

if USING_MPI:
    _comm = MPI.COMM_WORLD
    _rank = _comm.Get_rank()
    _size = _comm.Get_size()
    farm = _farm
    map = _map
    imap = _imap
    open = _open
    file = _file
else:
    _rank = 0
    _size = 1
    map = map
    from itertools import imap
    imap = imap
    farm = imap
    open = open
    file = file

class MPIFileHandler(StreamHandler):
    def __init__(self, filename, filemode='a'):
        self._file = open(filename, filemode)
        super(MPIFileHandler, self).__init__(self._file)

    def close(self):
        super(MPIFileHandler, self).close()
        self._file.close()

def _master(sequence):
    slave = 0
    for slave, item in enumerate(sequence, _WORKTAG):
        _comm.send(item, dest=slave, tag=_WORKTAG)
        if slave == _size-1:
            outstanding = slave
            break
    else:
        outstanding = slave
        for slave in range(slave+1, _size):
            _comm.send(None, dest=slave, tag=_DIETAG)

    status = MPI.Status()
    for item in sequence:
        result = _comm.recv(source=MPI.ANY_SOURCE, tag=_WORKTAG, status=status)
        slave = status.Get_source()
        _comm.send(item, dest=slave, tag=_WORKTAG)
        yield result

    for i in range(outstanding):
        result = _comm.recv(source=MPI.ANY_SOURCE, tag=_WORKTAG, status=status)
        slave = status.Get_source()
        _comm.send(None, dest=slave, tag=_DIETAG)
        yield result

def _slave(function):
    status = MPI.Status()
    while True:
        item = _comm.recv(source=_MASTER, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag() != _WORKTAG:
            return
        try:
            result = function(item)
        except:
            import traceback, sys
            sys.stderr.write('Uncaught exception:\n'+traceback.format_exc())
            _comm.Abort(1)
        _comm.send(result, dest=_MASTER, tag=_WORKTAG)

def main():
    return 0

if __name__ == '__main__':
    sys.exit(main())
