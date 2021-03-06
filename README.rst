#####################################################
Instructions for reproducing results of Kaehler et al
#####################################################

:authors: Benjamin D. Kaehler, Von Bing Yap, Rongli Zhang, and Gavin A. Huttley

**************************
Contents of the Repository
**************************

1. `src` contains two command line scripts that perform all model fitting and goodness-of-fit testing.
2. `figures` contains several scripts that process the output from the executables in `src` and produce figures and statistics.
3. `test` contains unit tests.
4. `lib` contains local library scripts.

************
Dependencies
************

Required
========

- Python_ 2.7.8
- `The latest PyCogent from github <https://github.com/pycogent/pycogent/archive/master.zip>`_
- numpy_ 1.8.2

.. _numpy: https://pypi.python.org/pypi/numpy/1.8.2
.. _python: http://https://www.python.org

Optional
========

- mpi4py_ 1.3.1
- nose_ 1.3.1
- rpy2_ 2.4.3
- R_ 3.1.1
    - ggplot2 1.0.0
    - stats 3.1.1
    - quantreg 5.05

.. _mpi4py: http://mpi4py.scipy.org
.. _nose: https://pypi.python.org/pypi/nose/1.3.4
.. _scipy: http://www.scipy.org
.. _rpy2: https://pypi.python.org/pypi/rpy2
.. _R: http://cran.r-project.org

********************
Installing this code
********************

All code is written to run in place.

Checking this code installed correctly
======================================

The unit tests will confirm that the model fitting and goodness-of-fit scripts are present with all the necessary requirements. From the top directory run::

    $ nosetests

or::

    $ mpiexec -n 3 nosetests

to confirm that the MPI functionality is working properly.

********************************************************
Reproduce the figures (plots and numbers) from the paper
********************************************************

1. Download the multiple sequence alignments from Dryad_.
2. Move the microbial multiple sequence alignments into the directory ``~/Data/greengenes/uniform/``
3. Move the mammal multiple sequence alignments into the directory ``~/Data/release-68/exons/aligned/``
4. Move the mitochondrial multiple sequence alignments into the directory ``~/Data/release-68/mtDNA/aligned/``
5. Change your working directory to src and run the following commands. (Launch theses scripts using MPI with a few hundred processes, or they will take a very, very long time to complete.)::

    $ python nonstationary_lengths.py -i ~/Data/release-68/exons/aligned -o ~/revisions/release-68/exons/aligned/general -l ~/revisions/release-68/exons/aligned/general/nsl.log -c 3 -u 20 -F seq_fit 
    $ python nonstationary_lengths.py -i ~/Data/release-68/exons/aligned -o ~/revisions/release-68/exons/aligned/gtrplusgamma -l ~/revisions/release-68/exons/aligned/gtrplusgamma/nsl.log -c 3 -u 20 -F hetero_fit
    $ python nonstationary_lengths.py -i ~/Data/release-68/exons/aligned -o ~/revisions/release-68/exons/aligned/clock/general -l ~/revisions/release-68/exons/aligned/clock/general/nsl.log -c 3 -u 20 -F clock_fit -O Opossum 
    $ python nonstationary_lengths.py -i ~/Data/release-68/exons/aligned -o ~/revisions/release-68/exons/aligned/clock/gtrplusgamma -l ~/revisions/release-68/exons/aligned/clock/gtrplusgamma/nsl.log -c 3 -u 20 -F hetero_clock_fit -O Opossum
    $ python nonstationary_lengths.py -i ~/Data/release-68/mtDNA/aligned -o ~/revisions/release-68/mtDNA/aligned/general -l ~/revisions/release-68/mtDNA/aligned/general/nsl.log -u 20 -F seq_fit 
    $ python nonstationary_lengths.py -i ~/Data/release-68/mtDNA/aligned -o ~/revisions/release-68/mtDNA/aligned/gtrplusgamma -l ~/revisions/release-68/mtDNA/aligned/gtrplusgamma/nsl.log -u 20 -F hetero_fit
    $ python nonstationary_lengths.py -i ~/Data/greengenes/uniform -o ~/revisions/greengenes/uniform/general -l ~/revisions/greengenes/uniform/general/nsl.log -u 20 -F seq_fit 
    $ python nonstationary_lengths.py -i ~/Data/greengenes/uniform/aligned -o ~/revisions/greengenes/uniform/gtrplusgamma -l ~/revisions/greengenes/uniform/gtrplusgamma/nsl.log -u 20 -F hetero_fit

    $ python g_stats.py -o ~/revisions/release-68/exons/aligned/general -l ~/revisions/release-68/exons/aligned/general/gs.log -N 100 -u 20 -P 1 -F seq_fit
    $ python g_stats.py -o ~/revisions/release-68/exons/aligned/general -l ~/revisions/release-68/exons/aligned/general/gs.log -N 100 -u 20 -P 0 -F seq_fit
    $ python g_stats.py -o ~/revisions/release-68/exons/aligned/gtrplusgamma -l ~/revisions/release-68/exons/aligned/gtrplusgamma/gs.log -N 100 -u 20 -P 0 -F hetero_fit
    $ python g_stats.py -o ~/revisions/release-68/exons/aligned/clock/general -l ~/revisions/release-68/exons/aligned/clock/general/gs.log -N 100 -u 20 -P 3 -F clock_fit -O Opossum 
    $ python g_stats.py -o ~/revisions/release-68/mtDNA/aligned/general -l ~/revisions/release-68/mtDNA/aligned/general/gs.log -N 100 -u 20 -P 1 -F seq_fit
    $ python g_stats.py -o ~/revisions/release-68/mtDNA/aligned/general -l ~/revisions/release-68/mtDNA/aligned/general/gs.log -N 100 -u 20 -P 0 -F seq_fit
    $ python g_stats.py -o ~/revisions/release-68/mtDNA/aligned/gtrplusgamma -l ~/revisions/release-68/mtDNA/aligned/gtrplusgamma/gs.log -N 100 -u 20 -P 0 -F hetero_fit
    $ python g_stats.py -o ~/revisions/greengenes/uniform/general -l ~/revisions/greengenes/uniform/general/gs.log -N 100 -u 20 -P 1 -F seq_fit
    $ python g_stats.py -o ~/revisions/greengenes/uniform/general -l ~/revisions/greengenes/uniform/general/gs.log -N 100 -u 20 -P 0 -F seq_fit
    $ python g_stats.py -o ~/revisions/greengenes/uniform/gtrplusgamma -l ~/revisions/greengenes/uniform/gtrplusgamma/gs.log -N 100 -u 20 -P 0 -F hetero_fit

6. Change your working directory to figures and run generate.sh. The plots will be saved in the figures directory as pdf files and the numbers will be output to screen.

.. _Dryad: http://doi:10.5061/dryad.g7g0n
