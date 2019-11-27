# -*- coding: utf-8 -*-

"""
other-lines.py

Analyse the immunopeptidomes of the non-JY cell lines

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
    plot_dir = fxn.plot_dir + fxn.get_date() + '-other-cell-lines/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    pep_data_dir = '../Data/peptide/'
    phospho_data_dir = '../Data/phospho/'

    nam_key = {'C': 'Cultured', 'M': 'Mouse'}
    hla_types = {'Colo205': ['HLA-A0101', 'HLA-A0201', 'HLA-B0702', 'HLA-B0801', 'HLA-C0702'],
                 'MDA-MB-436': ['HLA-A0101', 'HLA-B0801', 'HLA-C0702']}
    # NB: Colo205 and MDA-MB-436 actually have HLA-C*07:01, not :02, but MHCflurry doesn't cover 01

    # First do the non-phospho immunopeptidomes
    all_pep_files = [x for x in os.listdir(pep_data_dir) if '205' in x or '436' in x]
    peptides = coll.defaultdict(fxn.nest)

    fxn.initialise_prediction()

    for fl in all_pep_files:
        print fl
        growth_type = fl.split('.')[0].split('-')[0]
        cell_line = '-'.join(fl.split('.')[0].split('-')[1:])
        with open(pep_data_dir + fl, 'rU') as in_file:
            for line in in_file:
                peptides[cell_line][growth_type].append(line.rstrip())

    # Ensure these are unique lists
    lengths_raw = coll.defaultdict()

    for cell_line in peptides:
        totals = []
        lengths_raw[cell_line] = coll.defaultdict(fxn.nest_counter)
        peptide_lens_props = []

        for growth_type in peptides[cell_line]:
            peptides[cell_line][growth_type] = list(set(peptides[cell_line][growth_type]))
            totals.append([cell_line, nam_key[growth_type], len(peptides[cell_line][growth_type])])

            # Relative length distribution calculation
            for peptide in peptides[cell_line][growth_type]:
                lengths_raw[cell_line][growth_type][len(peptide)] += 1

            for i in range(8, 16):
                proportion = lengths_raw[cell_line][growth_type][i] / sum(lengths_raw[cell_line][growth_type].values())
                peptide_lens_props.append([cell_line, i, proportion, nam_key[growth_type]])

        # Relative length distribution plotting
        peptide_lens_props = fxn.list_to_df(peptide_lens_props, ['Cell Line', 'Length', 'Proportion', 'Type'], False)

        fig = plt.figure(figsize=(3, 5))
        ax = fig.add_subplot(111)
        sns.barplot(x='Length', y='Proportion', hue='Type', data=peptide_lens_props)
        plt.xlabel("Length (amino acids)")
        plt.savefig(plot_dir + cell_line + '-length-bars.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Total overlap Venn (unsampled)
        fxn.save_venn2(plot_dir + cell_line, peptides[cell_line],
                       'C', 'M', 'Cultured', 'Mouse')

        # Total number of peptides
        totals = fxn.list_to_df(totals, ['Cell Line', 'Growth Type', 'Number Peptides'], False)
        fig = plt.figure(figsize=(4, 5))
        ax = fig.add_subplot(111)
        sns.barplot(data=totals.loc[totals['Cell Line'] == cell_line], x='Growth Type', y='Number Peptides')
        save_name = plot_dir + 'total-peptides.png'
        plt.xlabel("")
        plt.savefig(plot_dir + cell_line + '-number-peptides.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Predict HLA binding affinities/alleles
    hla_predictions = coll.defaultdict()
    for cell_line in peptides:
        # Get all unique peptides for a given cell line into one list
        hla_predictions[cell_line] = []
        all_line_peptides = coll.Counter()
        for growth_type in peptides[cell_line]:
            for peptide in peptides[cell_line][growth_type]:
                all_line_peptides[peptide] += 1
        for pep in all_line_peptides:
            affin = 1e6
            allele = ''
            for cell_allele in hla_types[cell_line]:
                prediction = fxn.hla_prediction(pep, cell_allele)[0]
                if prediction < affin:
                    affin = prediction
                    allele = cell_allele
            if affin <= 50:
                assignment = '+++'
            elif affin <= 500:
                assignment = '++'
            elif affin <= 1000:
                assignment = '+'
            else:
                assignment = '-'
                allele = 'NPB'  # No predicted binder
            hla_predictions[cell_line].append([str(pep), allele, affin, assignment])
        hla_predictions[cell_line] = fxn.list_to_df(hla_predictions[cell_line],
                                                    ['Peptide', 'Allele', 'Affinity', 'Assignment'], True)

        # Having calculated the predicted HLA bindings for each of the unique peptides, count across the growth types
        hla_prop_dat = []
        for growth_type in peptides[cell_line]:
            hla_rows = hla_predictions[cell_line].loc[peptides[cell_line][growth_type]]
            for allele in hla_types[cell_line] + ['NPB']:
                if allele == 'NPB':
                    print_name = allele
                else:
                    print_name = fxn.tidy_hla_name(allele)
                subset = hla_rows.loc[hla_rows['Allele'] == allele]
                hla_prop_dat.append([print_name, len(subset) / len(hla_rows), nam_key[growth_type]])

        hla_prop_dat = fxn.list_to_df(hla_prop_dat, ['Allele', 'Proportion', 'Type'], False)

        # HLA allele contribution barplot
        fig = plt.figure(figsize=(6.5, 5))
        ax = fig.add_subplot(111)
        sns.barplot(x='Allele', y='Proportion', hue='Type', data=hla_prop_dat)
        plt.savefig(plot_dir + cell_line + '-hla-allele-contributions.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Finally analyse non-JY breast cancer line phospho-immunopeptidomes
    all_phospho_files = os.listdir(phospho_data_dir)
    line_key = {'MDAMB436': 'MDA-MB-436', 'Colo205': 'Colo205'}
    p_peptides = coll.defaultdict(fxn.nest)
    for cell_line in ['MDAMB436', 'Colo205']:
        tmp_C = [x for x in all_phospho_files if x.startswith('C-') and cell_line in x][0]
        tmp_M = [x for x in all_phospho_files if x.startswith('M-') and cell_line in x][0]
        for growth_type in ['C', 'M']:
            with open(phospho_data_dir + vars()['tmp_' + growth_type]) as in_file:
                for line in in_file:
                    p_peptides[line_key[cell_line]][nam_key[growth_type]].append(line.rstrip())

        fxn.save_venn2(plot_dir + 'phospho-' + line_key[cell_line], p_peptides[line_key[cell_line]],
                       'Cultured', 'Mouse', 'Cultured', 'Mouse')
