G=0.05
A=500
echo $G $A
set -e
set -x
twopi -Tpdf triad_example.dot > figure0.pdf
dot -Tpdf seq_fit.dot > figure1a.pdf
dot -Tpdf clock_fit.dot > figure1b.pdf
dot -Tpdf clock_fit_gamma.dot > figure1c.pdf
python overall_nt.py -f input_files/exons_combined -D -u -a $A -g $G -p '"Mouse", "Human"'
python human_gt_mouse.py -f input_files/clock_test -D -u -a $A -g $G
python stdev.py -f input_files/stdev -D -u -a $A -g $G -p '"Human", "Mouse"'
python akaike.py -f input_files/g_tests -D -u -a $A -g $G
python bar.py -f input_files/g_tests -D -u -a $A -o figure2.pdf
python jsd_plots_facets.py -f input_files/exons_combined -D -u -a $A -g $G -o figure3a.pdf
python jsd_plots_facets.py -f input_files/16s_combined -D -u -a $A -g $G -o figure3b.pdf
python gcb_facets.py  -g $G -a $A -D -u -f input_files/exons_combined -d ~/Data/release-68/exons/aligned -t -o figure4.pdf
python jsdlength.py -g $G -a $A -D -u -f input_files/exons -p '"Mouse", "Human"'  -o figure5.pdf
python paralinear.py -f input_files/paralinear -d ~/Data/release-68/exons/aligned -D -u -a $A -g $G -p '"Mouse", "Human"' -t -o figure6.pdf
python clock_test.py -f input_files/clock_test -D -u -a $A -g $G -o figure7.pdf
python n_t_ratios.py -f input_files/exons_gtr -D -u -a $A -g $G -p '"Mouse", "Human"' -o figure8a.pdf
python n_t_ratios.py -f input_files/exons_gtrplusgamma -D -u -a $A -g $G -p '"Mouse", "Human"' -o figure8b.pdf

python stdev.py -f input_files/stdev -D -u -a $A -p '"Human", "Mouse"' -G $G
python akaike.py -f input_files/g_tests -D -u -a $A -G $G
python jsd_plots_facets.py -f input_files/exons_combined -D -u -a $A -o online_figure3a.pdf -G $G
python jsd_plots_facets.py -f input_files/16s_combined -D -u -a $A -o online_figure3b.pdf -G $G
python gcb_facets.py  -a $A -D -u -f input_files/exons_combined -d ~/Data/release-68/exons/aligned -t -o online_figure4.pdf -G $G
python jsdlength.py -a $A -D -u -f input_files/exons -p '"Mouse", "Human"'  -o online_figure5.pdf -G $G
python paralinear.py -f input_files/paralinear -d ~/Data/release-68/exons/aligned -D -u -a $A -p '"Mouse", "Human"' -t -o online_figure6.pdf -G $G
python clock_test.py -f input_files/clock_test -D -u -a $A -o online_figure7.pdf -G $G
python n_t_ratios.py -f input_files/exons_gtr -D -u -a $A -p '"Mouse", "Human"' -o online_figure8a.pdf -G $G
python n_t_ratios.py -f input_files/exons_gtrplusgamma -D -u -a $A -p '"Mouse", "Human"' -o online_figure8b.pdf -G $G
