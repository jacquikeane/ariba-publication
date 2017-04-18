#!/usr/bin/env bash
set -euv


root_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
script_dir=$root_dir/Scripts

# _________________________ Run time and memory _______________________________
rm -rf $root_dir/Time_and_memory
mkdir $root_dir/Time_and_memory
$script_dir/get_time_and_memory.py $root_dir
cd $root_dir/Time_and_memory
R CMD BATCH $script_dir/plot_resources.R


# ________________________________ E. faecium _________________________________
cd $root_dir/E_faecium
rm -f summary.*
$script_dir/e_faecium_summary.py $PWD summary
$script_dir/e_faecium_summary.py --lenient $PWD summary.lenient
rm -f Rplots.pdf

# _________________________________ S. sonnei _________________________________
cd $root_dir/S_sonnei
rm -f summary.*
$script_dir/s_sonnei_summary.py $PWD Ref/srst2.fa Ref/ariba_db/ summary > summary.debug
$script_dir/s_sonnei_summary.py --lenient $PWD Ref/srst2.fa Ref/ariba_db/ summary.lenient > summary.lenient.debug
R CMD BATCH $script_dir/s_sonnei_depth_plots.R
rm -f Rplots.pdf


# _________________________________ N. gonorrhoeae ____________________________
cd  $root_dir/N_gonorrhoeae
rm -f summary.*

ls ARIBA_reports/* | awk -F/ '{print $0"\t"$2}' | sed 's/\.tsv$//' > summary.fofn

#Create summary for Phandango containing all variants for 23S and mtrR
ariba summary --row_filter n --cluster_cols assembled,pct_id,known_var,novel_var --only_clusters 23S,mtrR --v_groups --novel_variants --no_tree --fofn summary.fofn summary.Allvariants

#Create summary of known 23S and mtrR mutations and include assembled so that interrupted mtrR can be identified
ariba summary --row_filter n --cluster_cols assembled,known_var --only_clusters 23S,mtrR --v_groups --no_tree --fofn summary.fofn summary.AZMknowngroups

#Micplot of known AZM mutations including interrupted without combinations
ariba micplot $root_dir/N_gonorrhoeae/Ref/Ngo_ARIBAdb/ --no_combinations --interrupted --hlines 0.25,2 --dot_y_text_size 10 --point_size 0 --number_of_colours 3 Azithromycin mic_data.tsv summary.AZMknowngroups.csv summary.supp_fig.AZMknowngroups_nocombos

#MICplot of known AZM mutations including interrupted
ariba micplot $root_dir/N_gonorrhoeae/Ref/Ngo_ARIBAdb/ --interrupted --hlines 0.25,2 --dot_y_text_size 10 --point_size 0 --number_of_colours 3 Azithromycin mic_data.tsv summary.AZMknowngroups.csv summary.figure_4_AZMknowngroups

#MICplot of known AZM mutations including interrupted
ariba micplot $root_dir/N_gonorrhoeae/Ref/Ngo_ARIBAdb/ --use_hets exclude --interrupted --hlines 0.25,2 --dot_y_text_size 10 --point_size 0 --number_of_colours 3 Azithromycin mic_data.tsv summary.AZMknowngroups.csv summary.supp_fig.AZMknowngroups_nohets

# Figure 5
ariba summary --no_tree --col_filter y --row_filter n --cluster_cols assembled,ref_seq,pct_id,known_var,novel_var --known_variants --v_groups summary.fig5 ARIBA_reports/*.tsv

R CMD BATCH $script_dir/n_gonorrhoeae_figure5.R


# ____________________ copying figure pdf files and making xlsx _______________
cd $root_dir
rm -f Figures/figure_{2,3,4,5}.pdf Supplementary_LaTeX/pics-autogen
mkdir Supplementary_LaTeX/pics-autogen

cp E_faecium/summary.depth_plot.all.pdf Figures/figure_2.pdf
for x in B H R S X; do cp E_faecium/summary.depth_plot.Van$x.pdf Supplementary_LaTeX/pics-autogen/e_faecium.depth_plot.Van$x.pdf; done
cp E_faecium/summary.lenient.depth_plot.all.pdf Supplementary_LaTeX/pics-autogen/e_faecium.lenient.depth_plot.all.pdf
for x in B H R S X; do cp E_faecium/summary.lenient.depth_plot.Van$x.pdf Supplementary_LaTeX/pics-autogen/e_faecium.lenient.depth_plot.Van$x.pdf; done

cp S_sonnei/summary.s_sonnei.called_depth.0-150.pdf Supplementary_LaTeX/pics-autogen/s_sonnei.called_depth.0-150.pdf
cp S_sonnei/summary.s_sonnei.called_depth.pdf Supplementary_LaTeX/pics-autogen//s_sonnei.called_depth.pdf
cairosvg -o Figures/figure_3.pdf $root_dir/S_sonnei/summary.upset.svg
cairosvg -o Supplementary_LaTeX/pics-autogen/s_sonnei.lenient.upset.pdf $root_dir/S_sonnei/summary.lenient.upset.svg


cp N_gonorrhoeae/summary.figure_4_AZMknowngroups.pdf Figures/figure_4.pdf
cp N_gonorrhoeae/summary.supp_fig.AZMknowngroups_nocombos.pdf Supplementary_LaTeX/pics-autogen/n_gono_AZMknowngroups_nocombos.pdf
cp N_gonorrhoeae/summary.supp_fig.AZMknowngroups_nohets.pdf Supplementary_LaTeX/pics-autogen/n_gono_AZMknowngroups_nohets.pdf
cp N_gonorrhoeae/summary.Azithromycin_23S.23S.2597T_copy_number.pdf Figures/figure_5.pdf

cp Time_and_memory/*pdf Supplementary_LaTeX/pics-autogen/

./Scripts/make_xlsx.py


# ____________________________ making supplementary pdf _______________________
cd $root_dir/Supplementary_LaTeX/
pdflatex main.tex
pdflatex main.tex

#__ _____________________________ end _________________________________________
