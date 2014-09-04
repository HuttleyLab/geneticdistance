'''
    A topological sort algorithm based on a depth-first sort.

    According to Wikipedia it is described here:

    Cormen, Thomas H.; Leiserson, Charles E.; Rivest, Ronald L.; 
    Stein, Clifford (2001), "Section 22.4: Topological sort", 
    Introduction to Algorithms (2nd ed.), MIT Press and McGraw-Hill, 
    pp. 549-552, ISBN 0-262-03293-7.
'''

from collections import deque

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def sort(graph):
    L = deque()
    marked = set([])
    temporarily_marked = set([])

    def visit(n):
        assert n not in temporarily_marked, 'graph is not a DAG'
        if n not in marked:
            temporarily_marked.add(n)
            for m in graph[n]:
                visit(m)
            temporarily_marked.remove(n)
            marked.add(n)
            L.appendleft(n)

    for n in graph:
        visit(n)
    return L
