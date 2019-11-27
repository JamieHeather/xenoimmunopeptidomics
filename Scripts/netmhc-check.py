# -*- coding: utf-8 -*-

"""
netmhc-check.py

Validate the results of MHCflurry, i.e. whether mice do have fewer A2/more NPB peptides than cultured samples

Ran netMHC as follows (basically netMHC-4, default settings, relevant alleles):

$ Data/peptide/
for f in [CM]*txt
do
	base=$(echo $f | cut -d\- -f1)
	echo $base
	netMHC -a HLA-A0201,HLA-B0702,HLA-C0702 -p -f $f | grep 0\ \ \ \ HLA > netMHC/$base.netmhc.txt
done

"""

from __future__ import division
import functions as fxn
import os
import collections as coll
import seaborn as sns
import matplotlib.pyplot as plt


__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


if __name__ == '__main__':

    fxn.check_scripts_dir()
    sns.set(font="Arial", font_scale=1.5)

    # Sort directories, get data
    plot_dir = fxn.plot_dir + fxn.get_date() + '-netmhc-check/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    prediction_data_dir = '../Data/peptide/netMHC/'

    all_input_files = os.listdir(prediction_data_dir)
    all_input_files.sort()

    predictions = []
    names = []

    # Scroll through all netMHC produced files
    for f in all_input_files:
        nam = f.split('-')[0]
        print nam
        names.append(nam)

        # Pull out the necessary info, store in a temporary dictionary
        tmp_dict = coll.defaultdict(fxn.nest)
        with open(prediction_data_dir + f, 'rU') as in_file:
            for line in in_file:
                bits = [x for x in line.rstrip().split(' ') if x]
                allele = bits[1]
                peptide = bits[2]
                rank = float(bits[13])
                tmp_dict[peptide][allele] = rank

        # Then for each peptide, figure out best matching allele, and (if it's predicted to bind) notch one up for it
        allele_count = coll.Counter()
        for p in tmp_dict:
            best_allele = min(tmp_dict[p], key=tmp_dict[p].get)
            best_val = tmp_dict[p][best_allele]
            if best_val <= 2:
                if best_allele == 'HLA-A0201':
                    allele_count['HLA-A*\n02:01'] += 1
                elif best_allele == 'HLA-B0702':
                    allele_count['HLA-B*\n07:02'] += 1
                elif best_allele == 'HLA-C0702':
                    allele_count['HLA-C*\n07:02'] += 1
            else:
                allele_count['NPB'] += 1  # No predicted binder

        # Write out to nested list, to be turned into long proportions df, to be plotted as in peptide-analysis.py
        type_key = {'C': 'Cultured', 'M': 'Mouse'}
        for allele in allele_count:
            predictions.append([type_key[nam[0]], allele, allele_count[allele]/sum(allele_count.values())])

    props = fxn.list_to_df(predictions, ['Type', 'Allele', 'Proportion'], False)
    hla_alleles = ['HLA-A*\n02:01', 'HLA-B*\n07:02', 'HLA-C*\n07:02', 'NPB']

    # HLA allele contribution barplot
    fig = plt.figure(figsize=(4.5, 5))
    ax = fig.add_subplot(111)
    sns.barplot(data=props, x='Allele', y='Proportion', hue='Type', order=hla_alleles)

    sns.swarmplot(data=props, x='Allele', y='Proportion', hue='Type', order=hla_alleles,
                  dodge=True, edgecolor='black', alpha=.5, linewidth=1)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[2:], labels[2:], loc=1, borderaxespad=0.)  # Only plot the barplot legend items
    plt.savefig(plot_dir + 'hla-allele-contributions-netmhc.png', dpi=300, bbox_inches='tight')
    plt.close()
