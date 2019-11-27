# -*- coding: utf-8 -*-

"""
phospho-analysis.py

Analyse the JY phospho-immunopeptidomes

"""

from __future__ import division
import functions as fxn
import os
import collections as coll
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot
from scipy import stats

__version__ = '0.8.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def val_df_to_sharedness(value_df, headers):
    """
    NB only for use with the phospho data
    :param value_df: A pandas dataframe with samples for columns, peps/prots for rows, with vals = affinity/abundance
    :param headers: headers for final table
    :return: a longform dataframe for plotting
    """
    vals = []
    # Just in cultured
    for thing in [x for x in list(value_df.index) if value_df.loc[x][0] > 0 and sum(value_df.loc[x][1:]) == 0]:
        vals.append(['C-only', thing] + [x for x in value_df.loc[thing] if x > 0])

    # Just in one single mouse
    for thing in [x for x in list(value_df.index) if len([y for y in list(value_df.loc[x]) if y > 0]) == 1
                and value_df.loc[x][0] == 0]:
        vals.append(['1m-only', thing] + [x for x in value_df.loc[thing] if x > 0])

    # In cultured and X #s mice mouse only
    for i in range(2, 6):
        for thing in [x for x in list(value_df.index) if len([y for y in list(value_df.loc[x]) if y > 0]) == i
          and value_df.loc[x][0] > 0]:
            for val in [x for x in value_df.loc[thing] if x > 0]:
                vals.append(['C+' + str(i-1) + 'm', thing, val])

    return fxn.list_to_df(vals, headers, False)


if __name__ == '__main__':

    fxn.check_scripts_dir()
    sns.set(font="Arial", font_scale=1.5)

    # Sort directories, get data
    plot_dir = fxn.plot_dir + fxn.get_date() + '-phospho-analysis/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    pep_data_dir = '../Data/phospho/'
    prot_data_dir = '../Data/phospho-prot/'

    all_files = os.listdir(pep_data_dir)

    Mouse = [x for x in all_files if x.startswith('M') and '436' not in x and 'MCF7' not in x]
    Cultured = [x for x in all_files if x.startswith('C') and '436' not in x and 'MCF7' not in x]

    iterations = 100
    sample_to = 50

    pep_dat = coll.defaultdict(fxn.nest)
    prot_dat = coll.defaultdict(fxn.nest)
    experiments = coll.defaultdict(list)
    peptide_dict = coll.defaultdict(fxn.nest_counter)
    protein_dict = coll.defaultdict(fxn.nest_counter)
    type_key = {'C': 'Cultured', 'M': 'Mouse'}

    for growth_type in ['Mouse', 'Cultured']:

        for f in vars()[growth_type]:

            if not f.endswith('txt') and not f.endswith('tsv'):
                continue

            nam = f.split('-')[0]
            experiments[growth_type].append(nam)
            pep_dat[nam] = fxn.read_data_in(pep_data_dir + f)
            prot_dat[nam] = fxn.read_data_in(prot_data_dir + f)

            # Need to compile peptide/protein data into frames
            for peptide in pep_dat[nam]:
                peptide_dict[peptide][nam] = pep_dat[nam][peptide]
            for protein in prot_dat[nam]:
                protein_dict[protein][nam] = prot_dat[nam][protein]

    samples = pep_dat.keys()
    samples.sort()

    peptide_list = []
    for peptide in peptide_dict:
        # peptide_list.append([peptide, peptide_dict[peptide]['Hunt'], peptide_dict[peptide]['Cultured'],
        peptide_list.append([peptide, peptide_dict[peptide]['C'],
                             peptide_dict[peptide]['M01'], peptide_dict[peptide]['M04'],
                             peptide_dict[peptide]['M07'], peptide_dict[peptide]['M10']])

    peptides = fxn.list_to_df(peptide_list, ['Peptide', 'C', 'M01', 'M04', 'M07', 'M10'], True)

    protein_list = []
    for protein in protein_dict:
        protein_list.append([protein, protein_dict[protein]['C'],
                             protein_dict[protein]['M01'], protein_dict[protein]['M04'],
                             protein_dict[protein]['M07'], protein_dict[protein]['M10']])

    proteins = fxn.list_to_df(protein_list, ['Protein', 'C', 'M01', 'M04', 'M07', 'M10'], True)

    # Plot total number of peptides detected
    totals = peptides[peptides > 0].count()
    totals = pd.DataFrame({'Sample': list(totals.index), 'Total peptides': list(totals)})
    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    g = sns.catplot(x='Sample', y='Total peptides', data=totals, color=[fxn.blue] + [fxn.orange] * 4,
                    palette=[u'#4c72b0'] + [u'#dd8452'] * 4, kind='bar',
                    order=samples, aspect=0.6)
    g.set_xticklabels(rotation=85)
    plt.xlabel("")
    plt.ylabel("Unique phosphopeptides")
    plt.savefig(plot_dir + 'number-phosphopeptides.png', dpi=300, bbox_inches='tight')

    # Plot total abundance of peptides detected (should be the same as total number for peptides)
    total_abundance = pd.DataFrame({'Sample': samples, 'Total peptide abundance': peptides.sum(axis=0)})
    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    g = sns.catplot(x='Sample', y='Total peptide abundance', data=total_abundance,
                    color='lightslategray', kind='bar', order=samples, aspect=0.6,
                    palette=[u'#4c72b0'] + [u'#dd8452'] * 4)
    g.set_xticklabels(rotation=85)
    plt.xlabel("")
    plt.ylabel("Total phosphopeptide copies/cell")
    plt.savefig(plot_dir + 'total-abundance-phosphopeptides.png', dpi=300, bbox_inches='tight')

    # Loop through datasets 100 times, subsampling to a fixed number (< smallest dataset) and calculate Jaccard/R^2
    sampled_jacc = coll.defaultdict(fxn.nest)
    whole_plot_arr = np.zeros(shape=(len(samples), len(samples)))
    for i in range(iterations):

        print '\t' + str(i) + '\t-----------------------------------------------------'

        subbed_dat = {}
        for dataset in pep_dat:
            print '\t' + dataset
            subbed_dat[dataset] = fxn.subsample_to_number(pep_dat, dataset, sample_to)

        for x in range(len(samples)):
            for y in range(len(samples)):
                sampled_jacc[x][y].append(fxn.jaccard(subbed_dat[samples[x]], subbed_dat[samples[y]]))
                # Generate Venn plots of one examples of subsamples
                if i == 0:
                    fxn.save_venn2(plot_dir + 'sub', subbed_dat, samples[x], samples[y], samples[x], samples[y])

    sub_plot_arr = np.zeros(shape=(len(samples), len(samples)))
    for x in range(len(samples)):
        for y in range(len(samples)):
            sub_plot_arr[x, y] = np.mean(sampled_jacc[x][y])

            # Also calculate Jaccards of whole unsampled peptidomes
            whole_plot_arr[x, y] = fxn.save_venn2(plot_dir + 'whole', pep_dat,
                                                  samples[x], samples[y], samples[x], samples[y])

            if x != y:
                # And plot actual linear relationship
                f, ax = plt.subplots(figsize=(7, 7))
                sns.jointplot(samples[x], samples[y], data=peptides, kind='reg', stat_func=stats.pearsonr)
                plt.savefig(plot_dir + 'sub-Regplot-' + samples[x] + '-' + samples[y] + '.png',
                            dpi=300, bbox_inches='tight')
                plt.close()

    # Plot
    for arr in ['sub_', 'whole_']:
        temp_arr = vars()[arr + 'plot_arr']
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        p = ax.pcolor(temp_arr, cmap='gnuplot', vmin=0, vmax=np.max(temp_arr))
        ax.set_xticks([x + .5 for x in range(len(samples))])
        ax.set_xticklabels(samples, rotation=85)
        ax.set_yticks([x + .5 for x in range(len(samples))])
        ax.set_yticklabels(samples)  # , fontsize=4)
        plt.xlim(0, len(samples))
        plt.ylim(0, len(samples))
        cbar = plt.colorbar(p)
        cbar.ax.get_yaxis().labelpad = 10
        cbar.ax.set_ylabel(u'← less overlap    Jaccard index    more overlap →', rotation=90)
        plt.savefig(plot_dir + arr[:-1] + '-Jaccards-' + str(iterations) +
                    'repeats-' + str(sample_to) + 'sampled.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Generate an upsetplot showing the co-occurrence of peptides across all the sequences
    peptide_sets = {}
    for sample in pep_dat:
        peptide_sets[sample] = set(pep_dat[sample].keys())

    peptide_cooccurence = fxn.generate_upset_sets(peptide_sets)
    fig = plt.figure(figsize=(13, 6))
    ax = fig.add_subplot(111)
    upsetplot.plot(peptide_cooccurence, sort_by='cardinality', sort_sets_by='cardinality')
    plt.savefig(plot_dir + 'peptide-upset-plot.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Extract peptide lengths
    peptide_lens = []
    peptide_lens_props = []
    for sample in pep_dat:
        if sample == 'C':
            nam = sample
        else:
            nam = sample[:1]
        for i in range(8, 15):
            prop_len = len([x for x in pep_dat[sample] if len(x) == i]) / sum(pep_dat[sample].values())
            peptide_lens_props.append([sample, i, prop_len, nam])
        for pep in pep_dat[sample]:
            peptide_lens.append([sample, pep, len(pep), nam])

    # Convert to dataframe, then proportions
    peptide_lens = fxn.list_to_df(peptide_lens, ['Sample', 'Peptide', 'Length', 'Type'], False)
    peptide_lens_props = fxn.list_to_df(peptide_lens_props, ['Sample', 'Length', 'Proportion', 'Type'], False)
    prop = fxn.df_to_prop(peptide_lens, 'Length', 'Sample', 'Sample')

    # And plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    g = sns.FacetGrid(prop, row='Sample', size=3, aspect=1.2)
    g.map(sns.barplot, 'Length', 'Proportion')
    plt.savefig(plot_dir + 'length-bars.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(3, 5))
    ax = fig.add_subplot(111)
    sns.violinplot(x='Length', y='Proportion', data=peptide_lens_props, hue='Type', cut=0)
    sns.stripplot(x='Length', y='Proportion', hue='Type', data=peptide_lens_props, size=4,
                  dodge=True, edgecolor='black', alpha=.7, linewidth=.6, jitter=1.00000000000000005)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[:2], ['Cultured', 'Mouse'], loc=1, borderaxespad=0.)  # Only plot the violinplot legend items
    plt.savefig(plot_dir + 'length-violins.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.close('all')

    # Plot to see whether protein levels are conserved between the techs
    # Generate an upsetplot showing the co-occurrence of proteins across all the sequences
    protein_sets = {}
    for sample in prot_dat:
        protein_sets[sample] = set(prot_dat[sample].keys())

    protein_cooccurence = fxn.generate_upset_sets(protein_sets)
    fig = plt.figure(figsize=(13, 6))
    ax = fig.add_subplot(111)
    upsetplot.plot(protein_cooccurence, sort_by='cardinality', sort_sets_by='cardinality')
    plt.savefig(plot_dir + 'protein-upset-plot.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Print abundance of peptides/proteins based on how many samples they appear in
    col_order = ['C-only', '1m-only', 'C+1m', 'C+2m', 'C+3m', 'C+4m']
    short_col_order = ['C-only', '1m-only', 'C+4m']
    col_divisor = {'C-only': 1, '1m-only': 1, 'C+1m': 2, 'C+2m': 3, 'C+3m': 4, 'C+4m': 5}
    for p in ['peptides', 'proteins']:
        abundances = val_df_to_sharedness(vars()[p], ['Samples-in', p, 'Cumulative copies/cell'])

        fig = plt.figure(figsize=(8.5, 5))
        ax = fig.add_subplot(111)
        sns.boxplot(x='Samples-in', y='Cumulative copies/cell', data=abundances,
                    order=col_order, fliersize=0, color='lightslategray')
        for x in range(len(col_order)):
            num_to_write = int(len(abundances.loc[abundances['Samples-in'] == col_order[x]])/col_divisor[col_order[x]])
            ax.text(x, 10, num_to_write, horizontalalignment='center', size='x-small', color='black', weight='semibold')

        plt.ylim(0, 10)
        plt.xlabel("")
        plt.savefig(plot_dir + p + '-sharedness-abundances.png', dpi=300, bbox_inches='tight')
        plt.close()

        # And a small one
        fig = plt.figure(figsize=(3, 5))
        ax = fig.add_subplot(111)
        sns.boxplot(x='Samples-in', y='Cumulative copies/cell', data=abundances,
                    order=short_col_order, fliersize=0, palette=fxn.c_m_mix_cols)
        for x in range(len(short_col_order)):
            num_to_write = int(len(abundances.loc[abundances['Samples-in']
                                                  == short_col_order[x]])/col_divisor[short_col_order[x]])
            ax.text(x, 10, num_to_write, horizontalalignment='center', size='x-small', color='black', weight='semibold')

        plt.ylim(0, 10)
        plt.xlabel("")
        plt.savefig(plot_dir + p + '-sharedness-abundances-small.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Do stats on last group, i.e. protein comparisons
        just_c = list(abundances[abundances['Samples-in'] == 'C-only']['Cumulative copies/cell'])
        just_m = list(abundances[abundances['Samples-in'] == '1m-only']['Cumulative copies/cell'])
        both = list(abundances[abundances['Samples-in'] == 'C+4m']['Cumulative copies/cell'])

        stats.mannwhitneyu(just_c, just_m)
        # MannwhitneyuResult(statistic=271.0, pvalue=0.15670206311162221)
        stats.mannwhitneyu(just_c, both)
        # MannwhitneyuResult(statistic=1602.0, pvalue=0.0002979200322780869)
        stats.mannwhitneyu(just_m, both)
        # MannwhitneyuResult(statistic=1229.0, pvalue=7.028887963706661e-06)

    # Need to initialise MHCflurry to enable MHC prediction
    fxn.initialise_prediction()

    # Then predict the binding affinities to the best matched allele
    hla_predictions = []
    for pep in list(peptides.index):
        affin = 1e6
        allele = ''
        for jy_allele in ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']:
            prediction = fxn.hla_prediction(pep, jy_allele)[0]
            if prediction < affin:
                affin = prediction
                allele = jy_allele
        if affin <= 50:
            assignment = '+++'
        elif affin <= 500:
            assignment = '++'
        elif affin <= 1000:
            assignment = '+'
        else:
            assignment = '-'
            allele = 'NPB'  # No predicted binder
        hla_predictions.append([pep, allele, affin, assignment])

    # Go through each sample's peptides and record their HLA details
    hla_predictions = fxn.list_to_df(hla_predictions, ['Peptide', 'Allele', 'Affinity', 'Assignment'], True)
    sample_hla_predictions = []

    for sample in samples:
        for pep in pep_dat[sample]:
            sample_hla_predictions.append([sample] + list(hla_predictions.loc[pep]))

    # Convert to df, get proportions, then plot
    hla_alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702', 'NPB']
    sample_hla_predictions = fxn.list_to_df(sample_hla_predictions,
                                            ['Sample', 'Allele', 'Affinity', 'Assignment'], False)
    hla_prop = fxn.df_to_prop(sample_hla_predictions, 'Allele', 'Sample', 'Allele')

    sample_hla_props = []
    for sample in samples:
        counts = coll.Counter(list(sample_hla_predictions[sample_hla_predictions.Sample == sample]['Allele']))
        for hla in hla_alleles:
            if hla == 'HLA-A0201':
                print_hla = 'HLA-A*\n02:01'
            elif hla == 'HLA-B0702':
                print_hla = 'HLA-B*\n07:02'
            elif hla == 'HLA-C0702':
                print_hla = 'HLA-C*\n07:02'
            else:
                print_hla = 'NPB'
            sample_hla_props.append([sample, print_hla, counts[hla] / sum(counts.values()), type_key[sample[0]]])

    sample_hla_props = fxn.list_to_df(sample_hla_props, ['Sample', 'Allele', 'Proportion', 'Type'], False)
    hla_prop = fxn.df_to_prop(sample_hla_predictions, 'Allele', 'Sample', 'Allele')
    hla_prop = hla_prop.sort_values(by='Sample')

    # HLA allele contribution barplot
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    sns.barplot(x='Allele', y='Proportion', hue='Sample', data=hla_prop)
    plt.savefig(plot_dir + 'hla-allele-contributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    # HLA allele contribution barplot+swarmplot
    fig = plt.figure(figsize=(4.5, 5))
    ax = fig.add_subplot(111)
    sns.barplot(x='Allele', y='Proportion', hue='Type', data=sample_hla_props)
    sns.swarmplot(x='Allele', y='Proportion', hue='Type', data=sample_hla_props,
                  dodge=True, edgecolor='black', alpha=.5, linewidth=1)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[2:], labels[2:], loc=1, borderaxespad=0.)  # Only plot the barplot legend items
    # plt.setp(ax.get_xticklabels(), rotation=85)
    plt.savefig(plot_dir + 'hla-allele-contributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    # HLA affinity distributions by sample
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Sample', y='Affinity', data=sample_hla_predictions, color='lightslategray', fliersize=0)
    ax.set_yscale('log')
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    plt.ylim(1e1, 1e3)
    plt.savefig(plot_dir + 'hla-affinity-distribution.png', dpi=300, bbox_inches='tight')

    # Calculate whether there's any significant differences...
    affinities = coll.defaultdict(list)
    for sample in samples:
        affinities[sample] = list(sample_hla_predictions.loc[(sample_hla_predictions.Sample == sample), 'Affinity'])

    for sample1 in samples:
        for sample2 in samples:
            print sample1, sample2
            print stats.mannwhitneyu(affinities[sample1], affinities[sample2])
            # ... there aren't!

    # Compare affinity across different levels of peptide sharing
    sharedness = []
    for pep in list(peptides.index):
        number_samples = len([x for x in list(peptides.loc[pep]) if x > 0])
        sharedness.append([number_samples, hla_predictions.loc[pep, 'Affinity']])

    sharedness = fxn.list_to_df(sharedness, ['Number-Samples', 'Affinity'], False)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Number-Samples', y='Affinity', data=sharedness, fliersize=0, color='lightslategray')
    ax.set_yscale('log')
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    plt.ylim(1e1, 1e3)
    for x in range(5):
        ax.text(x, 1e3, len(sharedness.loc[sharedness['Number-Samples'] == x+1]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    plt.savefig(plot_dir + 'hla-sharedness-affinity.png', dpi=300, bbox_inches='tight')

    # Calculate whether there's any significant differences...
    affinities = coll.defaultdict(list)
    for x in range(1, 6):
        affinities[x] = list(sharedness.loc[sharedness['Number-Samples'] == x, 'Affinity'])

    for x in range(1, 6):
        for y in range(1, 6):
            mwu = stats.mannwhitneyu(affinities[x], affinities[y])
            if mwu[1] < 0.05:
                print x, y, mwu
                # No significant differences

    # Need to plot affinity based on sharedness with (or not with) cultured sample
    cultured_sharedness = []
    affinities = coll.defaultdict(list)

    # Just in cultured
    for pep in [x for x in list(peptides.index) if peptides.loc[x][0] > 0 and sum(peptides.loc[x][1:]) == 0]:
        cultured_sharedness.append(['C-only', hla_predictions.loc[pep, 'Affinity']])
        affinities['C-only'].append(hla_predictions.loc[pep, 'Affinity'])

    # Just in one single mouse
    for pep in [x for x in list(peptides.index) if len([y for y in list(peptides.loc[x])
                                                        if y > 0]) == 1 and peptides.loc[x][0] == 0]:
        cultured_sharedness.append(['1m-only', hla_predictions.loc[pep, 'Affinity']])
        affinities['1m-only'].append(hla_predictions.loc[pep, 'Affinity'])

    # In cultured and X #s mice mouse only
    for i in range(2, 6):
        for pep in [x for x in list(peptides.index) if len([y for y in list(peptides.loc[x])
                                                            if y > 0]) == i and peptides.loc[x][0] > 0]:
            cultured_sharedness.append(['C+' + str(i-1) + 'm', hla_predictions.loc[pep, 'Affinity']])
            affinities['C+' + str(i-1) + 'm'].append(hla_predictions.loc[pep, 'Affinity'])

    cultured_sharedness = fxn.list_to_df(cultured_sharedness, ['Samples-in', 'Affinity'], False)

    fig = plt.figure(figsize=(8.5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Affinity', data=cultured_sharedness, order=col_order,
                fliersize=0, color='lightslategray')
    ax.set_yscale('log')
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    for x in range(len(col_order)):
        ax.text(x, 1e3, len(cultured_sharedness.loc[cultured_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    plt.ylim(1e0, 1e4)
    plt.ylabel("MHCflurry predicted affinity (nM)")
    plt.savefig(plot_dir + 'hla-cultured-sharedness-affinity.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(3, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Affinity', data=cultured_sharedness, order=short_col_order,
                fliersize=0, palette=fxn.c_m_mix_cols)
    ax.set_yscale('log')
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    for x in range(len(short_col_order)):
        ax.text(x, 2e2, len(cultured_sharedness.loc[cultured_sharedness['Samples-in'] == short_col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    plt.ylim(1e0, 2e2)
    plt.ylabel("MHCflurry predicted affinity (nM)")
    plt.xlabel("")
    plt.savefig(plot_dir + 'hla-cultured-sharedness-affinity-small.png', dpi=300, bbox_inches='tight')
    plt.close()

    stats.mannwhitneyu(affinities['C-only'], affinities['1m-only'])
    # MannwhitneyuResult(statistic=270.0, pvalue=0.042483406569257226)
    stats.mannwhitneyu(affinities['C-only'], affinities['C+4m'])
    # MannwhitneyuResult(statistic=622.0, pvalue=0.06919762183621828)
    stats.mannwhitneyu(affinities['1m-only'], affinities['C+4m'])
    # MannwhitneyuResult(statistic=523.0, pvalue=0.18852470241650093)

    # See whether highly shared peptides derive from highly abundant RNA/proteins
    transcriptome, proteome = fxn.get_jy_ome_data()
    transcriptome_corr, proteome_corr = fxn.get_jy_ome_correlations(transcriptome, proteome, proteins)

    for d in list(proteins):

        # Transcriptomes
        tmp_df = transcriptome_corr[transcriptome_corr.Sample == d]
        fig = plt.figure(figsize=(5, 5))
        sns.jointplot(data=tmp_df, x='Peptidome Abundance', y='Transcriptome Abundance',
                      kind="reg", stat_func=stats.pearsonr)
        plt.xlabel("Immunopeptidome copies/cell")
        plt.ylabel("RNA abundance log2(TPM)")
        plt.savefig(plot_dir + 'transcriptome-correlation-' + d + '.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Proteomes
        tmp_df = proteome_corr[proteome_corr.Sample == d]
        fig = plt.figure(figsize=(5, 5))
        sns.jointplot(data=tmp_df, x='Peptidome Abundance', y='Proteome Abundance',
                      kind="reg", stat_func=stats.pearsonr)
        plt.xlabel("Immunopeptidome copies/cell")
        plt.ylabel("Protein abundance log2(intensity)")
        plt.savefig(plot_dir + 'proteome-correlation-' + d + '.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Need to plot affinity based on sharedness with (or not with) cultured sample
    prot_sharedness = []
    prot_abundances = coll.defaultdict(list)
    rna_sharedness = []
    rna_abundances = coll.defaultdict(list)

    # Just in cultured
    for prot in [x for x in list(proteins.index) if proteins.loc[x][0] > 0 and sum(proteins.loc[x][1:]) == 0]:
        if prot in proteome:
            prot_sharedness.append(['C-only', proteome[prot]])
            prot_abundances['C-only'].append(proteome[prot])
        if prot in transcriptome:
            rna_sharedness.append(['C-only', transcriptome[prot]])
            rna_abundances['C-only'].append(transcriptome[prot])

    # Just in one single mouse
    for prot in [x for x in list(proteins.index) if
                 len([y for y in list(proteins.loc[x]) if y > 0]) == 1 and proteins.loc[x][0] == 0]:
        if prot in proteome:
            prot_sharedness.append(['1m-only', proteome[prot]])
            prot_abundances['1m-only'].append(proteome[prot])
        if prot in transcriptome:
            rna_sharedness.append(['1m-only', transcriptome[prot]])
            rna_abundances['1m-only'].append(transcriptome[prot])

    # In cultured and X #s mice mouse only
    for i in range(2, 6):
        for prot in [x for x in list(proteins.index) if len([y for y in list(proteins.loc[x])
                                                             if y > 0]) == i and proteins.loc[x][0] > 0]:
            if prot in proteome:
                prot_sharedness.append(['C+' + str(i - 1) + 'm', proteome[prot]])
                prot_abundances['C+' + str(i - 1) + 'm'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['C+' + str(i - 1) + 'm', transcriptome[prot]])
                rna_abundances['C+' + str(i - 1) + 'm'].append(transcriptome[prot])

    prot_sharedness = fxn.list_to_df(prot_sharedness, ['Samples-in', 'Proteome Abundance'], False)
    rna_sharedness = fxn.list_to_df(rna_sharedness, ['Samples-in', 'Transcript Abundance'], False)

    stats.mannwhitneyu(prot_abundances['C-only'], prot_abundances['1m-only'])
    stats.mannwhitneyu(prot_abundances['C-only'], prot_abundances['C+4m'])
    stats.mannwhitneyu(prot_abundances['1m-only'], prot_abundances['1m-only'])
    # None significant

    stats.mannwhitneyu(rna_abundances['C-only'], rna_abundances['1m-only'])
    stats.mannwhitneyu(rna_abundances['C-only'], rna_abundances['C+4m'])
    stats.mannwhitneyu(rna_abundances['1m-only'], rna_abundances['1m-only'])
    # None significant

    fig = plt.figure(figsize=(8.5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Transcript Abundance', data=rna_sharedness, order=col_order,
                fliersize=0, color='lightslategray')
    for x in range(len(col_order)):
        ax.text(x, 9, len(rna_sharedness.loc[rna_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='RNA abundance log2(TPM)')
    plt.ylim(1, 9)
    plt.savefig(plot_dir + 'transcriptome-abundance-by-sharedness.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(8.5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Proteome Abundance', data=prot_sharedness, order=col_order,
                fliersize=0, color='lightslategray')
    for x in range(len(col_order)):
        ax.text(x, 36, len(prot_sharedness.loc[prot_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='Protein abundance log2(intensity)')
    plt.ylim(22, 36)
    plt.savefig(plot_dir + 'proteome-abundance-by-sharedness.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Transcript Abundance', data=rna_sharedness, order=short_col_order,
                fliersize=0, color='lightslategray')
    for x in range(len(short_col_order)):
        ax.text(x, 9, len(rna_sharedness.loc[rna_sharedness['Samples-in'] == short_col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='RNA abundance log2(TPM)')
    plt.ylim(1, 9)
    plt.savefig(plot_dir + 'transcriptome-abundance-by-sharedness-small.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Proteome Abundance', data=prot_sharedness, order=short_col_order,
                fliersize=0, color='lightslategray')
    for x in range(len(short_col_order)):
        ax.text(x, 36, len(prot_sharedness.loc[prot_sharedness['Samples-in'] == short_col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='Protein abundance log2(intensity)')
    plt.ylim(22, 36)
    plt.savefig(plot_dir + 'proteome-abundance-by-sharedness-small.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Output lists of all proteins that are featured in ALL mice samples or ALL cultured samples, for ORA
    ora = coll.defaultdict(list)

    for p in list(proteins.index):
        vals = list(proteins.loc[p])
        if vals[0] > 0:
            ora['in_cultured'].append(p)
        if len([x for x in vals[1:] if x > 0]) > 0:
            ora['in_mice'].append(p)

    fxn.save_venn2(plot_dir + 'ORA-lists', ora, 'in_mice', 'in_cultured', 'Mice', 'Cultured')
    ora_dir = '../Data/ORA/'
    for growth_type in ['mice', 'cultured']:
        out_file_path = ora_dir + 'proteins-in' + growth_type + '-phosphos.txt'
        with open(out_file_path, 'w') as out_file:
            for p in ora['in_' + growth_type]:
                out_file.write(p + '\n')

    # Then run these through WebGestalt's ORA process; can then unpack/read in the results and plot the correlation
    ora_results_dir = ora_dir + 'Results/'
    ora_go_p = coll.defaultdict(fxn.nest_counter)
    for growth_type in ['mice', 'cultured']:
        in_dat_dir = ora_results_dir + 'phospho-' + growth_type + '/'
        enrichment_file_path = [x for x in os.listdir(in_dat_dir) if x.startswith('enrichment')][0]
        with open(in_dat_dir + enrichment_file_path, 'rU') as in_file:
            line_count = 0
            for line in in_file:
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    headers = bits
                else:
                    ora_go_p[growth_type][bits[0]] = float(bits[7])
                line_count += 1

    fxn.save_venn2(plot_dir + 'ORA-GO-lists', ora_go_p, 'mice', 'cultured', 'Mice', 'Cultured')
    all_go = list(set(ora_go_p['mice']) | set(ora_go_p['cultured']))

    go_df = []
    for g in all_go:
        for growth_type in ['mice', 'cultured']:
            # In order to get sensible logged plotting, make non-present GO term Pvals = 1,
            # and '0' vals rounded down to order of magnitude just below lowest
            if g not in ora_go_p[growth_type]:
                ora_go_p[growth_type][g] = 1
            elif ora_go_p[growth_type][g] == 0:
                ora_go_p[growth_type][g] = 1e-16
        go_df.append([g, np.log10(ora_go_p['cultured'][g]), np.log10(ora_go_p['mice'][g])])

    go_headers = ['GO', 'log10 cultured sample GO term P val', 'log10 mice sample GO term P val']
    go_df = fxn.list_to_df(go_df, go_headers, True)

    fig = plt.figure(figsize=(5.5, 5))
    ax = fig.add_subplot(111)
    sns.jointplot(data=go_df, x='log10 cultured sample GO term P val',
                  y='log10 mice sample GO term P val', kind='reg', stat_func=stats.pearsonr)
    plt.savefig(plot_dir + 'ORA-1minus-correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
