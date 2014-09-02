# disable cogent.utils.parallel to stop it from conflicting with masterslave
import os
os.environ['DONT_USE_MPI'] = '1'
from warnings import filterwarnings
filterwarnings('ignore', 'Not using MPI', UserWarning)

import lib
import masterslave
import nest
import general_ben

from cogent import LoadSeqs, LoadTree
from cogent import DNA
from cogent.evolve.substitution_model import General, GeneralStationary
from cogent.evolve.models import GTR
import cogent

from itertools import combinations
from gzip import GzipFile
from socket import gethostname
from functools import partial
from traceback import format_exc
import logging
import json
import glob
import re
import sys
import os

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.9-dev'

_versions = {
        'nonstationary_length' : __version__,
        'cogent'               : cogent.__version__,
        'nest'                 : nest.__version__,
        'masterslave'          : masterslave.__version__,
        'general_ben'          : general_ben.__version__
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

def itriads(input_directory=None, codon_position=-1, output_directory=None,
        force_recalculation=False, triad_file=None, num_range=None, **kw):
    try:
        filenames = glob.glob(os.path.join(input_directory, '*.fasta*'))
        if num_range:
            filenames = filenames[num_range]
    except:
        logging.critical(' Unable to open input directory:\n'+format_exc())
        masterslave.exit(1)

    try:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory) # only executed by master
    except:
        logging.critical(' Unable to open output directory:\n'+format_exc())
        masterslave.exit(1)

    for filename in filenames:
        try:
            if filename[-6:] == '.fasta':
                gene = os.path.basename(filename[:-6])
                with open(filename) as fastafile:
                    fastadata = fastafile.read()
            elif filename[-9:] == '.fasta.gz':
                gene = os.path.basename(filename[:-9])
                with GzipFile(filename) as fastafile:
                    fastadata = fastafile.read()
            else:
                continue

            out_filepath = os.path.join(output_directory, gene + '.json')
            if not force_recalculation and os.path.exists(out_filepath):
                logging.info(' Skipping '+gene+': output exists')
                continue

            sequences = LoadSeqs(data=fastadata)
            if codon_position > 0:
                c = codon_position
                indices = [(i, i+1) for i in range(c-1, len(sequences), 3)]
                pos3 = sequences.addFeature('pos3', 'pos3', indices)
                sequences = pos3.getSlice()
            if triad_file:
                with open(triad_file) as f:
                    triads = [l.split() for l in f]
            else:
                triads = combinations(sequences.getSeqNames(), 3)

            num_triads = 0
            for triad in triads:
                for taxon in triad:
                    if taxon not in sequences.Names:
                        logging.info(' Skipping ' + '/'.join(triad) + ' in ' +
                                gene + ': ' + taxon + ' is missing')
                        break
                else:
                    num_triads += 1
                    sa = sequences.takeSeqs(triad)
                    sa = sa.filtered(lambda x: set(''.join(x)) <= set(DNA))
                    if len(sa) == 0:
                        logging.info(' Skipping ' + '/'.join(triad) + ' in ' + 
                                gene + ': filtered alignment has length zero')
                        continue
                    yield gene, sa

            if num_triads == 0:
                logging.info(' Skipping ' + gene + ': found no valid triads')

        except:
            logging.warning(' Skipping ' + filename + ':\n'+format_exc())

def process_triad((gene, sa), fitter=None, **kw):
    try: 
        st = LoadTree(tip_names=sa.Names)
        stats = fitter(sa, st, **kw)
        logging.debug(' Done '+'/'.join(sa.Names)+' in '+gene)
        return gene, stats
    
    except:
        logging.warning(' Skipping '+'/'.join(sa.Names)+' in '+gene+':\n'+
                format_exc())

def output_results(results, log_file=None, output_directory=None,
        pretty_print=None, **kw):
    if log_file:
        log_file = os.path.relpath(log_file, output_directory)
    else:
        log_file = 'stderr'
    gene_stats = []
    last_gene = None
    for result in results:
        try:
            if not result:
                continue
            gene, stats = result
            if gene != last_gene and gene_stats:
                gene_stats.append({'log_file': log_file})
                outfile_path = os.path.join(output_directory,last_gene+'.json')
                with open(outfile_path, 'w') as outfile:
                    pp = {'sort_keys':True, 'indent':1} if pretty_print else {}
                    json.dump(gene_stats, outfile, **pp)
                gene_stats = []
            last_gene = gene
            gene_stats += stats
        except:
            logging.error(' Problem collecting output for '+last_gene+':\n'+
                    format_exc())
    if gene_stats:
        try:
            gene_stats.append({'log_file': log_file})
            outfile_path = os.path.join(output_directory, last_gene+'.json')
            with open(outfile_path, 'w') as outfile:
                pp = {'sort_keys':True, 'indent':1} if pretty_print else {}
                json.dump(gene_stats, outfile, **pp)
        except:
            logging.error(' Problem collecting output for '+last_gene+':\n'+
                    format_exc())

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Estimate expected number of \
substitutions under stationary and non-stationary models.')
    parser.add_argument('-i', '--input_directory', required=True, 
            type=os.path.expanduser,
            help='Directory containing .fasta or .fasta.gz alignments')
    parser.add_argument('-o', '--output_directory', required=True,
            type=os.path.expanduser,
            help='Directory in which results will be stored')
    parser.add_argument('-l', '--log_file', type=os.path.expanduser, 
            help='Defaults to stderr')
    parser.add_argument('-L', '--log_level', default='INFO', type=str,
            help='Set logging level')
    parser.add_argument('-f', '--force_recalculation', action='store_true',
            help='Force recalculation of existing results')
    parser.add_argument('-c', '--codon_position', type=int, default=-1,
            help='Only use this codon position in analysis')
    parser.add_argument('-g', '--global_optimisation', action='store_true',
            help='Use global optimisation when fitting models')
    parser.add_argument('-R', '--num_range', type=lambda x: slice(*eval(x)), 
            default=None, help='Range of files to work on')
    parser.add_argument('-T', '--triad_file', type=os.path.expanduser,
            help='Tab delimited list of triads to be processed')
    parser.add_argument('-u', '--param_limit', type=float, default=20,
            help='Upper limit for substitution parameters')
    parser.add_argument('-p', '--pretty_print', action='store_true',
            help='Pretty print json output')
    parser.add_argument('-F', '--fitter', type=lambda x: getattr(nest, x),
            default='seq_fit', help='Nest function to use for fitting')
    parser.add_argument('-O', '--outgroup', type=str,
            help='Outgroup to use for clock test fitting')
    
    kw_args = vars(parser.parse_args())
    setup_logging(**kw_args)
    if masterslave.rank() == 0:
        logging.info(kw_args)
        logging.info(_versions)
    results = masterslave.imap(partial(process_triad, **kw_args), itriads(**kw_args))
    if results:
        output_results(results, **kw_args)
    return 0

if __name__ == '__main__':
    sys.exit(main())
