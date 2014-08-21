from __future__ import division

from scipy.stats import chisqprob
from nest import inflate_likelihood_function, _update

from gzip import GzipFile
from collections import defaultdict, Counter

import logging
import traceback
import json
import glob
import sys
import csv
import os

def parse_stats(from_directory, model_pos, stats=None, **kw):
    gene = None
    if stats is None:
        stats = defaultdict(lambda: defaultdict(list))
    for filename in glob.iglob(os.path.join(from_directory, '*.json*')):
        try:
            if filename.endswith('_g.json.gz') or filename.endswith('_g.json'):
                continue
            if filename.endswith('.json.gz'):
                gene = os.path.basename(filename[:-8])
                with GzipFile(filename) as jsonfile:
                    flat_stats = json.load(jsonfile)
            elif filename.endswith('.json'):
                gene = os.path.basename(filename[:-5])
                with open(filename) as jsonfile:
                    flat_stats = json.load(jsonfile)
            else:
                continue
            
            gene_stats = defaultdict(lambda: defaultdict(list))
            for row in flat_stats:
                if 'log_file' in row:
                    continue

                triad = frozenset(row['tip_names'])
                js = {frozenset(eval(k)): v for k, v in row['js'].items()}
                row['js'] = js
                row['gene'] = gene
                row['from_directory'] = from_directory
                gene_stats[triad][gene].append(row)
            
            for triad in gene_stats:
                for gene in gene_stats[triad]:
                    try:
                        stats[triad][gene].append(
                                gene_stats[triad][gene][model_pos])
                    except IndexError:
                        pass # clean them up later

            gene = None
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            logging.warning((' Skipping ' + gene + ':\n' if gene else ' ') +
                    traceback.format_exc())
    return stats

def clean_and_flatten(stats, checks):
    lived = Counter()
    died = Counter()
    length = max(len(_v) for k, v in stats.items() for _k, _v in v.items())
    def check(rows):
        if len(rows) == length:
            for check in checks:
                if not check(rows[0]):
                    return False
            return True
        return False
    for t, ts in stats.items():
        for v in ts.values():
            died[v[0]['from_directory']] += 1
        stats[t] = filter(check, ts.values())
        for v in stats[t]:
            lived[v[0]['from_directory']] += 1
            died[v[0]['from_directory']] -= 1
        if len(stats[t]) == 0:
            del stats[t]
    for fd in lived:
        print 'Got', lived[fd], 'from', fd
        print 'Lost', died[fd], 'from', fd

def lrt_p(lhn, lha, dfn, dfa):
    LR = 2*(lha - lhn)
    df = dfa - dfn
    return chisqprob(LR, df)

class Check(object):
    def __init__(self, name):
        self._c = Counter()
        self._n = name

    def __call__(self, row):
        result = self.check(row)
        if not result:
            self._c[row['from_directory']] += 1
        return result

    def __repr__(self):
        r = self._n + ':\n'
        fallen = []
        for k, v in self._c.items():
            fallen.append('    Lost ' + str(v) + ' from ' + k)
        r += '    Lost none' if len(fallen) == 0 else '\n'.join(fallen)
        return r

class UpperGCheck(Check):
    def __init__(self, level):
        super(UpperGCheck, self).__init__('G-Statistic Upper')
        self._l = level

    def check(self, row):
        return gs_p(row['gs_p']) <= self._l

class GCheck(Check):
    def __init__(self, level):
        super(GCheck, self).__init__('G-Statistic')
        self._l = level

    def check(self, row):
        return gs_p(row['gs_p']) > self._l

class LengthCheck(Check):
    def __init__(self, length):
        super(LengthCheck, self).__init__('Alignment Length')
        self._l = length

    def check(self, row):
        return row['aln_length'] >= self._l

class DLCCheck(Check):
    def __init__(self):
        super(DLCCheck, self).__init__('DLC')
    
    def check(self, row):
        return inflate_likelihood_function(row).isDLC()

class UniqueCheck(Check):
    def __init__(self):
        super(UniqueCheck, self).__init__('Generator Unique')

    def check(self, row):
        try:
            armu = inflate_likelihood_function(row).areRateMatricesUnique()
        except (ArithmeticError, NotImplementedError):
            logging.info(traceback.format_exc())
            return False
        return armu

def build_checks(g_sig_level=None, unique_check=False, DLC_check=False,
        aln_length=None, upper_g_sig_level=None, **kw):
    checks = []
    if aln_length:
        checks.append(LengthCheck(aln_length))
    if upper_g_sig_level:
        checks.append(UpperGCheck(upper_g_sig_level))
    if g_sig_level:
        checks.append(GCheck(g_sig_level))
    if DLC_check:
        checks.append(DLCCheck())
    if unique_check:
        checks.append(UniqueCheck())
    return checks

def gs_p(gs):
    return 1 - (gs[0] + 101 - gs[1])/101

def get_results():
    import argparse
    parser = argparse.ArgumentParser(
            description='Analyse nonstaionary_length output')
    parser.add_argument('-d', '--data_directory', 
            help='Directory containing original alignments')
    parser.add_argument('-f', '--input_directory_file', required=True,
            type = os.path.expanduser,
            help='File containing list of input directories and model positions')
    parser.add_argument('-o', '--output_file',
            help='File into which results will be placed. Defaults to stdout')
    parser.add_argument('-O', '--outliers', type=float, default=0.,
            help='For a linear regression, remove this proportion of samples')
    parser.add_argument('-l', '--log_file', help='Defaults to stderr')
    parser.add_argument('-g', '--g_sig_level', type=float,
            help='Reject results with G-stat p-value less than this')
    parser.add_argument('-G', '--upper_g_sig_level', type=float,
            help='Reject results with G-stat p-value greater than this')
    parser.add_argument('-a', '--aln_length', type=int,
            help='Reject results with G-stat p-value less than this')
    parser.add_argument('-p', '--pair', type=lambda x: frozenset(eval(x)),
            help='Restrict analysis to this pair of taxa')
    parser.add_argument('-T', '--triad', type=lambda x: frozenset(eval(x)),
            help='Restrict analysis to this triad')
    parser.add_argument('-m', '--model',
            help='Restrict analysis to this model')
    parser.add_argument('-t', '--third_position', action='store_true',
            help='Only use the third codon position in analysis')
    parser.add_argument('-D', '--DLC_check', action='store_true',
            help='Reject results with any P matrices that do not satisfy DLC')
    parser.add_argument('-u', '--unique_check', action='store_true',
            help='Only keep results with Q matrices that map uniquely to P')
    parser.add_argument('-r', '--time_tree', type=os.path.expanduser,
            help='File containing time tree information')
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(message)s',
            filename=args.log_file)
    
    with open(args.input_directory_file) as in_file:
        reader = csv.reader(in_file, delimiter='\t')
        input_dirs = [(os.path.expanduser(d), int(n)) for d, n in reader]

    checks = build_checks(**vars(args))
    stats = None
    for d in input_dirs:
        stats = parse_stats(*d, stats=stats, **vars(args))
    clean_and_flatten(stats, checks)
    for check in checks:
        print check
    
    return stats, args

def main():
    stats, args = get_results()

    print len(stats)
    print vars(args)

     
if __name__ == '__main__':
    main()

