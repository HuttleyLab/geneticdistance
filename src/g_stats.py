# disable cogent.utils.parallel to stop it from conflicting with masterslave
import os
os.environ['DONT_USE_MPI'] = '1'

import masterslave
import nest

from collections import defaultdict

import cogent

from itertools import combinations
from gzip import GzipFile
from functools import partial
from traceback import format_exc
from socket import gethostname
import random
import time
import logging
import json
import numpy as np
import glob
import sys
import os

__version__ = '0.0.9-dev'

_versions = {
        'g_stats'     : __version__,
        'cogent'      : cogent.__version__,
        'nest'        : nest.__version__,
        'masterslave' : masterslave.__version__
        }

def setup_logging(log_level=None, log_file=None, **kw):
    try:
        if log_file:
            log_dir = os.path.dirname(log_file)
            masterslave.checkmakedirs(log_dir)
            handler = masterslave.MPIFileHandler(log_file)
        else:
            handler = logging.StreamHandler()
        log_level = getattr(logging, log_level.upper())
        handler.setLevel(log_level)
        hostpid = ''
        if masterslave.USING_MPI:
            hostpid = gethostname()+':'+str(os.getpid())+':'
        formatter = logging.Formatter('%(asctime)s:'+hostpid+
                '%(levelname)s:%(message)s')
        handler.setFormatter(formatter)
        logging.root.addHandler(handler)
        logging.root.setLevel(log_level)
    except:
        sys.stderr.write(' Unable to set up logging:\n'+format_exc())
        masterslave.exit(1)

def istats(output_directory=None, mode='s', **kw):
    try:
        filenames = glob.iglob(os.path.join(output_directory, '*.json'))
    except:
        logging.critical(' Unable to open ' + output_directory)
        masterslave.exit(1)
    for filename in filenames:
        try:
            if filename[-7:] == '_g.json':
                continue
            gene = os.path.basename(filename[:-5])

            with open(filename) as f_file:
                f_stats = json.load(f_file)

            g_stats = []
            g_filepath = os.path.join(output_directory, gene + '_g.json')
            if os.path.exists(g_filepath):
                if mode == 's':
                    logging.info(' Skipping ' + gene + ': output exists')
                    continue
                elif mode == 'a':
                    with open(g_filepath) as g_file:
                        g_stats = json.load(g_file)
            
            last_triad = None
            stats = []
            for f_row in f_stats:
                if 'tip_names' not in f_row:
                    continue
                triad = frozenset(f_row['tip_names'])
                if last_triad != triad:
                    last_triad = triad
                    stats.append([])
                stats[-1].append(f_row)

            if g_stats == []:
                g_stats = [{} for i in range(len(stats))]

            if len(g_stats) != len(stats):
                logging.error(' Skipping ' + filename + ':\n' + ' found ' +
                        str(len(stats)) + ' triad(s) and ' + str(len(g_stats)) + 
                        ' bootstrap records(s). Numbers should match.')
                continue

            for f_stats, g_stats in zip(stats, g_stats):
                yield gene, f_stats, g_stats

        except:
            logging.warning(' Skipping '+filename+':\n'+format_exc())

def param_bootstrap(stats, num_reps=None, model_pos=None, fitter=None, **kw):
    gene, f_stats, g_stats = stats
    try:
        f_row = f_stats[model_pos]
    except IndexError:
        logging.error(' Skipping ' + '/'.join(f_stats[0]['tip_names']) + ' in '
                + gene + ': position ' + model_pos + ' invalid')
        return
    model = f_row['name']

    if model_pos in g_stats:
        g_row = g_stats[model_pos]
    else:
        g_row = {'name'     : model,
                'tip_names' : f_row['tip_names'],
                'gs_samples': [],
                'll_samples': [],
                'en_samples': []}
        g_stats[model_pos] = g_row
    gs_samples = g_row['gs_samples']
    ll_samples = g_row['ll_samples']
    en_samples = g_row['en_samples']
    if 'state' in g_row:
        random.setstate(eval(g_row['state']))

    lf = nest.inflate_likelihood_function(f_row)
    aln_length = f_row['aln_length']
    start = time.time()
    for i in 10*range(num_reps):
        if len(gs_samples) >= num_reps:
            break
        try:
            aln = lf.simulateAlignment(aln_length, random_series=random)
            lfs = fitter(aln, lf.tree, return_lfs=model, **kw)
            fitted_lf = lfs[model_pos]
            ll_samples.append(fitted_lf.getLogLikelihood())
            gs_samples.append(fitted_lf.getGStatistic())
            if 'Q' in fitted_lf.defn_for:
                en_samples.append(nest.get_expected_no_subs(fitted_lf))
        except:
            logging.warning(' Missed a G stat for ' + model + ' and ' +
                    '/'.join(f_row['tip_names']) + ' in ' + gene + ':\n' +
                    format_exc())
    else:
        logging.error(' Failed to compile sufficient bootstrap repetitions for '
                + model + ' and ' + '/'.join(f_row['tip_names']) + ' in ' + gene)
    g_row['state'] = repr(random.getstate())
    f_row['gs_p'] = (sum(1 for g in gs_samples if g < f_row['gs']),
        len(gs_samples) + 1)
    f_row['ll_p'] = (sum(1 for l in ll_samples if l < f_row['ll']),
        len(ll_samples) + 1)
    logging.info(' Done ' + model + ' and ' + '/'.join(f_row['tip_names']) + 
            ' in ' + gene + ' in ' + str(time.time() - start) + ' secs')

    return gene, f_stats, g_stats

def write_files(gene, f_stats, g_stats, output_directory, log_file):
    try:
        f_stats.append({'log_file': log_file})
        filepath = os.path.join(output_directory, gene+'.json')
        with open(filepath) as infile:
            old_stats = json.load(infile)
            for row in old_stats:
                if 'log_file' in row:
                    f_stats.append(row)
        with open(filepath, 'w') as outfile:
            json.dump(f_stats, outfile)
        filepath = os.path.join(output_directory, gene+'_g.json')
        with open(filepath, 'w') as outfile:
            json.dump(g_stats, outfile)
    except:
        if gene:
            logging.error(' Problem collecting output for '+gene+':\n'+
                    format_exc())
        else:
            logging.critical(' Null gene in write_files()')
            masterslave.exit(1)

def output_results(results, log_file=None, output_directory=None, **kw):
    if log_file:
        log_file = os.path.relpath(log_file, output_directory)
    else:
        log_file = 'stderr'
    
    f_stats = []
    g_stats = []
    last_gene = None
    for result in results:
        if not result:
            continue
        gene, f_s, g_s = result
        if gene != last_gene and f_stats:
            write_files(last_gene, f_stats, g_stats, output_directory, log_file)
            f_stats = []
            g_stats = []
        last_gene = gene
        f_stats += f_s
        g_stats.append(g_s)
           
    if f_stats:
        write_files(gene, f_stats, g_stats, output_directory, log_file)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Estimate expected number of \
substitutions under stationary and non-stationary models.')
    parser.add_argument('-o', '--output_directory', required=True,
            type=os.path.expanduser,
            help='Directory in which results will be stored')
    parser.add_argument('-l', '--log_file', type=os.path.expanduser,
            help='Defaults to stderr')
    parser.add_argument('-L', '--log_level', default='INFO', 
            help='Set logging level')
    parser.add_argument('-m', '--mode', default='s', choices='sra',
            help='Write mode (skip, replace, append)')
    parser.add_argument('-R', '--num_range', type=lambda x: slice(*eval(x)),
            default=None, help='Range of files to work on')
    parser.add_argument('-N', '--num_reps', type=int, default=100,
            help='Number of G stat samples')
    parser.add_argument('-u', '--param_limit', type=float, default=20,
            help='Upper limit for length and substitution parameters')
    parser.add_argument('-P', '--model_pos', default=-1, type=int,
            help='Model for which to produce G statistics')
    parser.add_argument('-F', '--fitter', type=lambda x: getattr(nest, x),
            default='seq_fit', help='Nest function to use for fitting')
    parser.add_argument('-O', '--outgroup', type=str,
            help='Outgroup to use for clock test fitting')

    kw_args = vars(parser.parse_args())
    setup_logging(**kw_args)
    if masterslave.rank() == 0:
        logging.info(kw_args)
        logging.info(_versions)
    processor = partial(param_bootstrap, **kw_args)
    results = masterslave.imap(processor, istats(**kw_args))
    if results:
        output_results(results, **kw_args)
    return 0

if __name__ == '__main__':
    main()

