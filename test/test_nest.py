from __future__ import division

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from nose.tools import assert_equal, assert_less
from numpy import array, sqrt, allclose
import sys

from cogent.util.misc import ConstraintError
from cogent.evolve.substitution_model import General, GeneralStationary
from cogent.evolve.models import GTR

import lib
import nest
from data import get_aln

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.22-dev'

def test_expected_no_subs():
    """expected_no_subs should return -int_0^t p0 exp(Qt) dt diag(Q)"""
    result = 0.7295333
    Q = array([[0.5, 0.2, 0.1, 0.2]]*4)
    for i in range(4):
        Q[i, i] = 0.
        Q[i, i] = -sum(Q[i,])
    p0 = array([0.2, 0.3, 0.3, 0.2])
    t = 1.
    assert_almost_equal(result, nest.expected_no_subs(p0, Q, t))

def test_get_expected_no_subs():
    """expected_no_subs should return dictionary of ENS by edge"""
    GS_lf = nest.inflate_likelihood_function(_GeneralStationary)
    EN = nest.get_expected_no_subs(GS_lf)
    for name in GS_lf.tree.getTipNames():
        assert_almost_equal(EN[name], GS_lf.getParamValue('length', name))

def test_inflate_deflate_likelihood_function():
    """deflate/inflate_likelihood_function are reciprocal maps"""
    lf = nest.inflate_likelihood_function(_GTRplusGamma)
    aln = get_aln('GTRplusGamma', _GTRplusGamma['aln_length'])
    lf.setAlignment(aln)

    down = nest.deflate_likelihood_function(lf)
    down_up = nest.inflate_likelihood_function(down)
    down_up.setAlignment(aln)
    down_up_down = nest.deflate_likelihood_function(down_up)
    
    assert_equal(down, down_up_down)

def test_deflate_likelihood_function():
    """deflate_likelihood_function produces internally consistent statistics"""
    lf = nest.inflate_likelihood_function(_General)
    aln = get_aln('General', _General['aln_length'])
    lf.setAlignment(aln)
    EN = nest.deflate_likelihood_function(lf)['EN']
    assert_equal(EN, nest.get_expected_no_subs(lf))

def test_seq_fit():
    """seq_fit should fit nested GTR and General models"""
    for model in 'GTR', 'General':
        pre_lf = nest.inflate_likelihood_function(eval('_'+model))
        prefit = nest.get_expected_no_subs(pre_lf)
        
        aln = get_aln(model, 100000)
        lfs = nest.seq_fit(aln, pre_lf.tree, param_limit=20, return_lfs=model)
        if model == 'General':
            assert_less(lfs[0].getLogLikelihood(),
                    lfs[1].getLogLikelihood())
        lf = lfs[-1]
        postfit = nest.get_expected_no_subs(lf)
        for taxon in prefit:
            assert_almost_equal(postfit[taxon], prefit[taxon], decimal=2) 

def test_hetero_fit():
    """hetero_fit should fit GTR plus Gamma models"""
    pre_lf = nest.inflate_likelihood_function(_GTRplusGamma)
    prefit = nest.get_expected_no_subs(pre_lf)
    aln = get_aln('GTRplusGamma', 100000)
    lfs = nest.hetero_fit(aln, pre_lf.tree, param_limit=20, return_lfs=True)
    postfit = nest.get_expected_no_subs(lfs[-1])
    for taxon in prefit:
        assert_almost_equal(postfit[taxon], prefit[taxon], decimal=2) 
    
def test_hetero_clock_fit():
    """hetero_clock_fit should fit a molecular clock constrained GTR plus Gamma
    model nested in a GTR plus Gamma model"""
    pre_lf = nest.inflate_likelihood_function(_GTRplusGammaClockTest)
    prefit = nest.get_expected_no_subs(pre_lf)
    aln = get_aln('GTRplusGammaClockTest', 100000)
    lfs = nest.hetero_clock_fit(aln, pre_lf.tree, outgroup='Opossum',
            param_limit=20, return_lfs=True)
    lf_equal_length, lf = lfs
    assert_less(lf_equal_length.getLogLikelihood(), lf.getLogLikelihood())
    postfit = nest.get_expected_no_subs(lf)
    postfit_equal_length = nest.get_expected_no_subs(lf_equal_length)
    for taxon in prefit:
        assert_almost_equal(postfit[taxon], prefit[taxon], decimal=2) 
        assert_almost_equal(postfit_equal_length[taxon], prefit[taxon], decimal=2) 

def test_clock_fit():
    """clock_fit should fit nested GTR, General, and GeneralBen models,
    some with equal branch lengths"""
    for modelname in ('GTRClockTest', 'GeneralBen'):
        model = eval('_' + modelname)
        pre_lf = nest.inflate_likelihood_function(model)
        prefit = nest.get_expected_no_subs(pre_lf)
        aln = get_aln(modelname, 100000)
        lfs = nest.clock_fit(aln, pre_lf.tree, outgroup='Opossum', param_limit=20, 
                return_lfs='GTR' if modelname.startswith('GTR') else 'General')
        lf_equal_length, lf = lfs[:2] if modelname[:3] == 'GTR' else lfs[2:]
        assert_less(lf_equal_length.getLogLikelihood(), lf.getLogLikelihood())
        if modelname == 'GeneralBen':
            assert_less(lfs[0].getLogLikelihood(),
                        lf_equal_length.getLogLikelihood())
        postfit = nest.get_expected_no_subs(lf)
        postfit_equal_length = nest.get_expected_no_subs(lf_equal_length)
        for taxon in prefit:
            assert_almost_equal(postfit[taxon], prefit[taxon], decimal=2) 
            assert_almost_equal(postfit_equal_length[taxon], prefit[taxon], 
                    decimal=2) 

def test_populate_parameters():
    """populate_parameters should set up a nested likelihood function"""
    lf_GTR = nest.inflate_likelihood_function(_GTR)
    lf_General = nest.inflate_likelihood_function(_General)
    for edge in lf_GTR.tree.getTipNames():
        assert not allclose(lf_GTR.getPsubForEdge(edge),
                lf_General.getPsubForEdge(edge)), 'models started close'
    nest.populate_parameters(lf_General, lf_GTR)
    for edge in lf_GTR.tree.getTipNames():
        assert_array_almost_equal(lf_GTR.getPsubForEdge(edge),
                lf_General.getPsubForEdge(edge))

def generate_alignments():
    from gzip import GzipFile
    from data import get_data_dir
    from os.path import join
    alns = [('GTRplusGamma', _GTRplusGamma['aln_length']),
            ('General', _General['aln_length']),
            ('GTR', 100000), ('General', 100000),
            ('GTRplusGamma', 100000), ('GTRplusGammaClockTest', 100000),
            ('GTRClockTest', 100000), ('GeneralBen', 100000)]
    alns = [('GTRClockTest', 100000), ('GeneralBen', 100000)]
    for model, aln_len in alns:
        lf = nest.inflate_likelihood_function(eval('_' + model))
        aln = lf.simulateAlignment(aln_len)
        filename = '_'.join((model, str(aln_len))) + '.fasta.gz'
        with GzipFile(join(get_data_dir(), filename), 'w') as aln_file:
            aln_file.write(aln.toFasta())
    return 0

def main():
#    generate_alignments()
#    return 0

    test_populate_parameters()
    print 'done test_populate_parameters'
    test_deflate_likelihood_function()
    print 'done test_deflate_likelihood_function'
    test_inflate_deflate_likelihood_function()
    print 'done test_inflate_deflate_likelihood_function'
    test_get_expected_no_subs()
    print 'done test_get_expected_no_subs'
    test_expected_no_subs()
    print 'done test_expected_no_subs'
    test_seq_fit()
    print 'done test_seq_fit'
    test_clock_fit()
    print 'done test_clock_fit'
    test_hetero_fit()
    print 'done test_hetero_fit'
    test_hetero_clock_fit()
    print 'done test_hetero_clock_fit'

_GTR = {'EN': {'Human': 0.16115955204676646,
  'Mouse': 0.20707912415193264,
  'Opossum': 0.5705441379355761},
 'aln_length': 1031,
 'df': 11,
 'gs': 88.62960090332447,
 'gs_p': [400, 401],
 'js': {"('Mouse', 'Human')": 6.897788687165729e-05,
  "('Opossum', 'Human')": 0.0002508522059996299,
  "('Opossum', 'Mouse')": 0.00034500250182167846,
  "('Opossum', 'Mouse', 'Human')": 0.00029555333898034775},
 'll': -3254.2589369487173,
 'll_p': [258, 401],
 'mprobs': {'A': 0.1548831923848898,
  'C': 0.2942615106104341,
  'G': 0.32981023466256343,
  'T': 0.22104506234211266},
 'name': 'GTR',
 'params': {'A/C': {'Human': {'bin0': 2.3393306975927097},
   'Mouse': {'bin0': 2.3393306975927097},
   'Opossum': {'bin0': 2.3393306975927097}},
  'A/G': {'Human': {'bin0': 9.344059564212516},
   'Mouse': {'bin0': 9.344059564212516},
   'Opossum': {'bin0': 9.344059564212516}},
  'A/T': {'Human': {'bin0': 3.2032765719404717},
   'Mouse': {'bin0': 3.2032765719404717},
   'Opossum': {'bin0': 3.2032765719404717}},
  'C/G': {'Human': {'bin0': 0.8533134041203625},
   'Mouse': {'bin0': 0.8533134041203625},
   'Opossum': {'bin0': 0.8533134041203625}},
  'C/T': {'Human': {'bin0': 7.9388705060058715},
   'Mouse': {'bin0': 7.9388705060058715},
   'Opossum': {'bin0': 7.9388705060058715}},
  'length': {'Human': 0.16115955204676646,
   'Mouse': 0.20707912415193264,
   'Opossum': 0.5705441379355761}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}

_GeneralStationary = {'EN': {'Human': 0.17607207857188917,
  'Mouse': 0.19116991136366723,
  'Opossum': 0.5666588424785894},
 'aln_length': 1031,
 'df': 30,
 'gs': 67.18125031647061,
 'js': {"('Mouse', 'Human')": 6.897788687165729e-05,
  "('Opossum', 'Human')": 0.0002508522059996299,
  "('Opossum', 'Mouse')": 0.00034500250182167846,
  "('Opossum', 'Mouse', 'Human')": 0.00029555333898034775},
 'll': -3243.433751359517,
 'mprobs': {'A': 0.15566482092906736,
  'C': 0.2938116602212761,
  'G': 0.32971926086041814,
  'T': 0.2208042579892385},
 'name': 'GeneralStationary',
 'params': {'A/C': {'Human': {'bin0': 1.0000002592239173e-06},
   'Mouse': {'bin0': 0.8979246962472756},
   'Opossum': {'bin0': 0.22189637671268791}},
  'A/T': {'Human': {'bin0': 0.20081813074737317},
   'Mouse': {'bin0': 0.04522781553687011},
   'Opossum': {'bin0': 0.7366662970979752}},
  'C/A': {'Human': {'bin0': 1.0000000061914567e-06},
   'Mouse': {'bin0': 0.45611798642647716},
   'Opossum': {'bin0': 0.27570514601674384}},
  'C/G': {'Human': {'bin0': 0.06302750384413647},
   'Mouse': {'bin0': 0.12977939240502762},
   'Opossum': {'bin0': 0.13193396653952372}},
  'C/T': {'Human': {'bin0': 0.5552423331167203},
   'Mouse': {'bin0': 1.7385870037984923},
   'Opossum': {'bin0': 0.7608842850832422}},
  'T/A': {'Human': {'bin0': 0.20242616918841963},
   'Mouse': {'bin0': 0.3994091387668378},
   'Opossum': {'bin0': 0.5481866444320843}},
  'T/C': {'Human': {'bin0': 0.5094164523004237},
   'Mouse': {'bin0': 1.379764061138782},
   'Opossum': {'bin0': 0.8887865419442014}},
  'T/G': {'Human': {'bin0': 0.11103386998525465},
   'Mouse': {'bin0': 0.2742576468769069},
   'Opossum': {'bin0': 0.05764699023138122}},
  'length': {'Human': 0.17607207857188917,
   'Mouse': 0.1911699113636672,
   'Opossum': 0.5666588424785894}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}

_General = {'EN': {'Human': 0.17394311948884797,
  'Mouse': 0.18337633546820462,
  'Opossum': 0.5277313281367219},
 'aln_length': 1031,
 'df': 39,
 'gs': 38.91419249155015,
 'gs_p': [350, 401],
 'js': {"('Mouse', 'Human')": 6.897788687165729e-05,
  "('Opossum', 'Human')": 0.0002508522059996299,
  "('Opossum', 'Mouse')": 0.00034500250182167846,
  "('Opossum', 'Mouse', 'Human')": 0.00029555333898034775},
 'll': -3229.623539877141,
 'll_p': [208, 401],
 'mprobs': {'A': 0.1123444708566534,
  'C': 0.2813768679430518,
  'G': 0.370030286770501,
  'T': 0.2362483744297938},
 'name': 'General',
 'params': {'A/C': {'Human': {'bin0': 1.0000000328426524e-06},
   'Mouse': {'bin0': 0.3925393588891999},
   'Opossum': {'bin0': 0.1478553375135336}},
  'A/G': {'Human': {'bin0': 0.17893527744728566},
   'Mouse': {'bin0': 0.01570457062959878},
   'Opossum': {'bin0': 0.47498498782711635}},
  'A/T': {'Human': {'bin0': 0.03862827058765341},
   'Mouse': {'bin0': 0.020403184760213428},
   'Opossum': {'bin0': 0.4409610661755024}},
  'C/A': {'Human': {'bin0': 0.014350770864713126},
   'Mouse': {'bin0': 0.24855021600952792},
   'Opossum': {'bin0': 0.1983117026588957}},
  'C/G': {'Human': {'bin0': 0.03647011721194428},
   'Mouse': {'bin0': 0.06297026512348824},
   'Opossum': {'bin0': 0.07117691546340325}},
  'C/T': {'Human': {'bin0': 0.31956378066059915},
   'Mouse': {'bin0': 0.8473736112274767},
   'Opossum': {'bin0': 0.34392815285161416}},
  'G/C': {'Human': {'bin0': 0.04021915329707346},
   'Mouse': {'bin0': 0.076577088869986},
   'Opossum': {'bin0': 0.07944341799359622}},
  'G/T': {'Human': {'bin0': 0.043063939489069114},
   'Mouse': {'bin0': 0.04952836073964198},
   'Opossum': {'bin0': 0.049593659340448584}},
  'T/A': {'Human': {'bin0': 0.15032451543501174},
   'Mouse': {'bin0': 0.251887229659776},
   'Opossum': {'bin0': 0.517848807166159}},
  'T/C': {'Human': {'bin0': 0.23874438553864197},
   'Mouse': {'bin0': 0.5212130058653592},
   'Opossum': {'bin0': 0.7411793015168008}},
  'T/G': {'Human': {'bin0': 0.052123648102814564},
   'Mouse': {'bin0': 0.10980957763178682},
   'Opossum': {'bin0': 0.049862645271916126}},
  'length': {'Human': 0.17636419803484066,
   'Mouse': 0.1841994732886849,
   'Opossum': 0.527834100897092}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}

_GeneralBen = {'EN': {'Human': 0.16787723594561116,
  'Mouse': 0.16787723594561105,
  'Opossum': 0.5559953925462553},
 'aln_length': 1031,
 'dependencies': {'A/C': [{'edges': ['Mouse']},
   {'edges': ['Opossum']},
   {'edges': ['Human']}],
  'A/G': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'A/T': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'C/A': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'C/G': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'C/T': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'G/C': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'G/T': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'T/A': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'T/C': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'T/G': [{'edges': ['Mouse']}, {'edges': ['Opossum']}, {'edges': ['Human']}],
  'length': [{'edges': ['Mouse', 'Human']}, {'edges': ['Opossum']}]},
 'df': 38,
 'gs': 34.81057786032392,
 'js': {"('Mouse', 'Human')": 0.00047733890135637225,
  "('Opossum', 'Human')": 0.001786606133242996,
  "('Opossum', 'Mouse')": 0.0030992419482862577,
  "('Opossum', 'Mouse', 'Human')": 0.0023905001310877694},
 'll': -3281.684054230924,
 'mprobs': {'A': 0.12175323021721553,
  'C': 0.2873114084342382,
  'G': 0.356104146936461,
  'T': 0.2348312144120853},
 'name': 'GeneralBen',
 'params': {'A/C': {'Human': {'bin0': 0.01739536754166563},
   'Mouse': {'bin0': 0.31745089230190804},
   'Opossum': {'bin0': 0.031160939029535328}},
  'A/G': {'Human': {'bin0': 0.1589263536565361},
   'Mouse': {'bin0': 1.0000005129275605e-06},
   'Opossum': {'bin0': 0.6696459364015507}},
  'A/T': {'Human': {'bin0': 0.08262312582993601},
   'Mouse': {'bin0': 0.08780938814878483},
   'Opossum': {'bin0': 0.7067280076115665}},
  'C/A': {'Human': {'bin0': 0.019613342496517496},
   'Mouse': {'bin0': 0.3207076944885724},
   'Opossum': {'bin0': 0.34817557223609363}},
  'C/G': {'Human': {'bin0': 0.028648944395845095},
   'Mouse': {'bin0': 0.06648480555681702},
   'Opossum': {'bin0': 0.09847467575121126}},
  'C/T': {'Human': {'bin0': 0.36969234213962504},
   'Mouse': {'bin0': 0.8896999113146077},
   'Opossum': {'bin0': 0.30236440806287757}},
  'G/C': {'Human': {'bin0': 0.02482482251587386},
   'Mouse': {'bin0': 0.1478740509395621},
   'Opossum': {'bin0': 0.13928635020135502}},
  'G/T': {'Human': {'bin0': 0.045949784308164696},
   'Mouse': {'bin0': 0.03664525483387962},
   'Opossum': {'bin0': 0.08670575951485372}},
  'T/A': {'Human': {'bin0': 0.25398048829732006},
   'Mouse': {'bin0': 0.306993614058762},
   'Opossum': {'bin0': 0.430888826580508}},
  'T/C': {'Human': {'bin0': 0.31591561978289817},
   'Mouse': {'bin0': 0.32715765459617424},
   'Opossum': {'bin0': 0.7685261161283372}},
  'T/G': {'Human': {'bin0': 0.06826071583310062},
   'Mouse': {'bin0': 0.10931678856385355},
   'Opossum': {'bin0': 0.03554895005731518}},
  'length': {'Human': 0.16787723594561108,
   'Mouse': 0.16787723594561108,
   'Opossum': 0.5559953925462553}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}

_GTRClockTest = {'EN': {'Human': 0.18620418050573892,
  'Mouse': 0.18620418050573892,
  'Opossum': 0.5024388154721869},
 'dependencies': {'A/C': [{'edges': ['Mouse', 'Opossum', 'Human']}],
  'A/G': [{'edges': ['Mouse', 'Opossum', 'Human']}],
  'A/T': [{'edges': ['Mouse', 'Opossum', 'Human']}],
  'C/G': [{'edges': ['Mouse', 'Opossum', 'Human']}],
  'C/T': [{'edges': ['Mouse', 'Opossum', 'Human']}],
  'length': [{'edges': ['Mouse', 'Human']}, {'edges': ['Opossum']}]},
 'df': 10,
 'mprobs': {'A': 0.15949328415361136,
  'C': 0.30686718265087093,
  'G': 0.3219823362062739,
  'T': 0.2116571969892438},
 'name': 'GTR',
 'params': {'A/C': {'Human': {'bin0': 1.0160711625703538},
   'Mouse': {'bin0': 1.0160711625703538},
   'Opossum': {'bin0': 1.0160711625703538}},
  'A/G': {'Human': {'bin0': 8.01180076472396},
   'Mouse': {'bin0': 8.01180076472396},
   'Opossum': {'bin0': 8.01180076472396}},
  'A/T': {'Human': {'bin0': 4.944114919997411},
   'Mouse': {'bin0': 4.944114919997411},
   'Opossum': {'bin0': 4.944114919997411}},
  'C/G': {'Human': {'bin0': 1.6936045022784383},
   'Mouse': {'bin0': 1.6936045022784383},
   'Opossum': {'bin0': 1.6936045022784383}},
  'C/T': {'Human': {'bin0': 8.575093670136177},
   'Mouse': {'bin0': 8.575093670136177},
   'Opossum': {'bin0': 8.575093670136177}},
  'length': {'Human': 0.1862041805057389,
   'Mouse': 0.1862041805057389,
   'Opossum': 0.5024388154721869}},
 'tip_names': ['Mouse', 'Opossum', 'Human'],
 'with_rate': False}

_GTRplusGamma = {'EN': {'Human': 0.17417499333872027,
  'Mouse': 0.23382018203154617,
  'Opossum': 0.57285966730339655},
 'aln_length': 1000,
 'dependencies': {'A/C': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'A/G': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'A/T': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'C/G': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'C/T': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'bprobs': [{'loci': ['locus0']}],
  'length': [{'edges': ['Mouse']},
   {'edges': ['Human']},
   {'edges': ['Opossum']}],
  'mprobs': [{'edges': ['root', 'Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'rate_shape': [{}]},
 'df': 15,
 'gs': 131.92929398977799,
 'js': {"('Mouse', 'Human')": 0.00031578008667842994,
  "('Opossum', 'Human')": 0.002957388508205261,
  "('Opossum', 'Mouse')": 0.002089472266039394,
  "('Opossum', 'Mouse', 'Human')": 0.0023909587890782458},
 'll': -3261.8051461517257,
 'mprobs': {'A': 0.15887312725601699,
  'C': 0.29836773363959851,
  'G': 0.32556410224792659,
  'T': 0.21719503685645786},
 'name': 'GTR',
 'params': {'A/C': {'Human': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218},
   'Mouse': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218},
   'Opossum': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218}},
  'A/G': {'Human': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624},
   'Mouse': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624},
   'Opossum': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624}},
  'A/T': {'Human': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905},
   'Mouse': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905},
   'Opossum': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905}},
  'C/G': {'Human': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712},
   'Mouse': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712},
   'Opossum': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712}},
  'C/T': {'Human': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478},
   'Mouse': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478},
   'Opossum': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478}},
  'length': {'Human': 0.17417499333872027,
   'Mouse': 0.23382018203154617,
   'Opossum': 0.57285966730339655},
  'rate_shape': 1.0},
 'tip_names': ['Human', 'Mouse', 'Opossum'],
 'tree': '(Human,Mouse,Opossum)root;',
 'with_rate': True}

_GTRplusGammaClockTest = {'EN': {'Human': 0.23382018203154617,
  'Mouse': 0.23382018203154617,
  'Opossum': 0.57285966730339655},
 'aln_length': 1000,
 'dependencies': {'A/C': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'A/G': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'A/T': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'C/G': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'C/T': [{'bins': ['bin1', 'bin0', 'bin3', 'bin2'],
    'edges': ['Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'bprobs': [{'loci': ['locus0']}],
  'length': [{'edges': ['Mouse']},
   {'edges': ['Human']},
   {'edges': ['Opossum']}],
  'mprobs': [{'edges': ['root', 'Opossum', 'Mouse', 'Human'],
    'loci': ['locus0']}],
  'rate_shape': [{}]},
 'df': 15,
 'gs': 131.92929398977799,
 'js': {"('Mouse', 'Human')": 0.00031578008667842994,
  "('Opossum', 'Human')": 0.002957388508205261,
  "('Opossum', 'Mouse')": 0.002089472266039394,
  "('Opossum', 'Mouse', 'Human')": 0.0023909587890782458},
 'll': -3261.8051461517257,
 'mprobs': {'A': 0.15887312725601699,
  'C': 0.29836773363959851,
  'G': 0.32556410224792659,
  'T': 0.21719503685645786},
 'name': 'GTR',
 'params': {'A/C': {'Human': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218},
   'Mouse': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218},
   'Opossum': {'bin0': 2.7163265190242218,
    'bin1': 2.7163265190242218,
    'bin2': 2.7163265190242218,
    'bin3': 2.7163265190242218}},
  'A/G': {'Human': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624},
   'Mouse': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624},
   'Opossum': {'bin0': 8.3106755354521624,
    'bin1': 8.3106755354521624,
    'bin2': 8.3106755354521624,
    'bin3': 8.3106755354521624}},
  'A/T': {'Human': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905},
   'Mouse': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905},
   'Opossum': {'bin0': 2.641225566697905,
    'bin1': 2.641225566697905,
    'bin2': 2.641225566697905,
    'bin3': 2.641225566697905}},
  'C/G': {'Human': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712},
   'Mouse': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712},
   'Opossum': {'bin0': 0.84442704955960712,
    'bin1': 0.84442704955960712,
    'bin2': 0.84442704955960712,
    'bin3': 0.84442704955960712}},
  'C/T': {'Human': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478},
   'Mouse': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478},
   'Opossum': {'bin0': 7.8368242144958478,
    'bin1': 7.8368242144958478,
    'bin2': 7.8368242144958478,
    'bin3': 7.8368242144958478}},
  'length': {'Human': 0.23382018203154617,
   'Mouse': 0.23382018203154617,
   'Opossum': 0.57285966730339655},
  'rate_shape': 1.0},
 'tip_names': ['Human', 'Mouse', 'Opossum'],
 'tree': '(Human,Mouse,Opossum)root;',
 'with_rate': True}

if __name__ == '__main__':
    sys.exit(main())
