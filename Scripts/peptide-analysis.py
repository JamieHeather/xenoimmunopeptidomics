# -*- coding: utf-8 -*-

"""
peptide-analysis.py

Analyse the standard JY immunopeptidomes

"""

from __future__ import division

import functions as fxn
import os
import collections as coll
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

__version__ = '0.4.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


if __name__ == '__main__':

    fxn.check_scripts_dir()
    sns.set(font="Arial", font_scale=1.5)

    # Sort directories, get data
    plot_dir = fxn.plot_dir + fxn.get_date() + '-peptide-analysis/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    pep_data_dir = '../Data/peptide/'
    prot_data_dir = '../Data/peptide-prot/'

    all_files = os.listdir(pep_data_dir)

    # Only take the JY samples
    Mouse = [x for x in all_files if x.startswith('M') and '205' not in x and '436' not in x]
    Cultured = [x for x in all_files if x.startswith('C0') and '205' not in x and '436' not in x]

    iterations = 100
    sample_to = 2500

    pep_dat = coll.defaultdict(fxn.nest)
    prot_dat = coll.defaultdict(fxn.nest)
    experiments = coll.defaultdict(list)
    peptide_dict = coll.defaultdict(fxn.nest_counter)
    protein_dict = coll.defaultdict(fxn.nest_counter)

    for growth_type in ['Mouse', 'Cultured']:

        for f in vars()[growth_type]:

            if not f.endswith('txt') and not f.endswith('tsv'):
                continue

            nam = f.split('.')[0].split('-')[0]
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
        peptide_list.append([peptide] + [peptide_dict[peptide][x] for x in samples])

    peptides = fxn.list_to_df(peptide_list, ['Peptide'] + samples, True)

    protein_list = []
    for protein in protein_dict:
        protein_list.append([protein.split('\t')[0]] + [protein_dict[protein][x] for x in samples])

    proteins = fxn.list_to_df(protein_list, ['Protein'] + samples, True)

    # Plot total number of peptides detected
    totals = peptides[peptides > 0].count()
    totals = pd.DataFrame({'Sample': list(totals.index), 'Total peptides': list(totals)})
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    g = sns.catplot(x='Sample', y='Total peptides', data=totals, color='lightslategray', kind='bar', order=samples,
                    palette=[u'#4c72b0'] * 3 + [u'#dd8452'] * 10)
    g.set_xticklabels(rotation=85)

    # Define and plot group averages
    c_y = np.mean(totals.loc[totals.Sample.str.startswith('C')]['Total peptides'])
    m_y = np.mean(totals.loc[totals.Sample.str.startswith('M')]['Total peptides'])
    plt.plot([-0.5, 2.5], [c_y, c_y], c=fxn.blue)
    plt.plot([2.5, 12.5], [m_y, m_y], c=fxn.orange)

    plt.xlabel("")
    plt.savefig(plot_dir + 'number-peptides.png', dpi=300, bbox_inches='tight')

    # Loop through datasets 100 times, subsampling to a fixed number (< smallest dataset) and calculate Jaccard/R^2
    sampled_jacc = coll.defaultdict(fxn.nest)
    sampled_jacc_by_type = coll.defaultdict(list)
    whole_plot_arr = np.zeros(shape=(len(samples), len(samples)))
    for i in range(iterations):

        # print '\t' + str(i) + '\t-----------------------------------------------------'

        subbed_dat = {}
        for dataset in samples:
            # print '\t' + dataset
            subbed_dat[dataset] = fxn.subsample_to_number(pep_dat, dataset, sample_to)

        done = []
        for x in range(len(samples)):
            for y in range(len(samples)):
                sampled_jacc[x][y].append(fxn.jaccard(subbed_dat[samples[x]], subbed_dat[samples[y]]))

                # Generate Venn plots of one examples of subsamples
                if i == 0:
                    fxn.save_venn2(plot_dir + 'sub', subbed_dat, samples[x], samples[y], samples[x], samples[y])

                # Make lists of unique Jaccards in simple lists for non-mapped plotting
                to_pair = [samples[x], samples[y]]
                to_pair.sort()
                group = '-'.join([p[:-2] for p in to_pair])
                pair = samples[x] + '-' + samples[y]

                if pair not in done and x != y:
                    sampled_jacc_by_type[group].append(sampled_jacc[x][y][-1])

                done = done + [pair, samples[y] + '-' + samples[x]]

    sub_plot_arr = np.zeros(shape=(len(samples), len(samples)))
    for x in range(len(samples)):
        for y in range(len(samples)):
            sub_plot_arr[x, y] = np.mean(sampled_jacc[x][y])

            # Also calculate Jaccards of whole unsampled peptidomes
            whole_plot_arr[x, y] = fxn.save_venn2(plot_dir + 'whole', pep_dat,
                                                  samples[x], samples[y], samples[x], samples[y])

    # Plot
    for arr in ['sub_', 'whole_']:
        temp_arr = vars()[arr + 'plot_arr']
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        p = ax.pcolor(temp_arr, cmap='gnuplot', vmin=0, vmax=np.max(temp_arr))
        ax.set_xticks([x + .5 for x in range(len(samples))])
        ax.set_xticklabels(samples, rotation=85)
        ax.set_yticks([x + .5 for x in range(len(samples))])
        ax.set_yticklabels(samples)  # fontsize=4)
        plt.xlim(0, len(samples))
        plt.ylim(0, len(samples))
        cbar = plt.colorbar(p)
        cbar.ax.get_yaxis().labelpad = 10
        cbar.ax.set_ylabel(u'← less overlap         Jaccard index         more overlap →', rotation=90)
        plt.savefig(plot_dir + arr[:-1] + '-Jaccards-' + str(iterations) +
                    'repeats-' + str(sample_to) + 'sampled.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Extract peptide lengths
    type_key = {'C': 'Cultured', 'M': 'Mouse', 'B': 'Both'}
    peptide_lens = []
    peptide_lens_props = []
    for sample in pep_dat:
        for i in range(8, 16):
            prop_len = len([x for x in pep_dat[sample] if len(x) == i]) / sum(pep_dat[sample].values())
            peptide_lens_props.append([sample, i, prop_len, type_key[sample[:-2]]])
        for pep in pep_dat[sample]:
            peptide_lens.append([sample, pep, len(pep), type_key[sample[:-2]]])

    # Convert to dataframe, then proportions
    peptide_lens = fxn.list_to_df(peptide_lens, ['Sample', 'Peptide', 'Length', 'Type'], False)
    peptide_lens_props = fxn.list_to_df(peptide_lens_props, ['Sample', 'Length', 'Proportion', 'Type'], False)
    prop = fxn.df_to_prop(peptide_lens, 'Length', 'Sample', 'Sample')

    # And plot
    fig = plt.figure(figsize=(3, 5))
    ax = fig.add_subplot(111)
    sns.barplot(x='Length', y='Proportion', hue='Type', data=peptide_lens_props)
    sns.swarmplot(x='Length', y='Proportion', hue='Type', data=peptide_lens_props, size=4,
                   dodge=True, edgecolor='black', alpha=0, linewidth=.6)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[2:], labels[2:], loc=1, borderaxespad=0.)  # Only plot the barplot legend items
    plt.xlabel("Length (amino acids)")
    plt.savefig(plot_dir + 'length-bars.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.close('all')

    # Then predict the binding affinities to the best matched allele (if it hasn't been calculated before)
    hla_files = [x for x in all_files if x.endswith('-peptide-hla-predictions.csv')]
    if not hla_files:
        # Need to initialise MHCflurry to enable MHC prediction
        fxn.initialise_prediction()
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

        # And save that to save time on future runs
        hla_predictions.to_csv(pep_data_dir + fxn.get_date() + '-peptide-hla-predictions.csv')

    else:
        hla_files.sort()
        hla_predictions = pd.read_csv(pep_data_dir + hla_files[-1])
        hla_predictions = hla_predictions.set_index('Peptide', drop=True)

    # Turn those into sample specific values
    sample_hla_predictions = []
    for sample in samples:
        for pep in pep_dat[sample]:
            sample_hla_predictions.append([sample] + list(hla_predictions.loc[pep]) + [type_key[sample[0]]])

    # Convert to df, get proportions, then plot
    hla_alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702', 'NPB']
    sample_hla_predictions = fxn.list_to_df(sample_hla_predictions,
                                            ['Sample', 'Allele', 'Affinity', 'Assignment', 'Type'], False)

    sample_hla_props = []
    for sample in samples:
        counts = coll.Counter(list(sample_hla_predictions[sample_hla_predictions.Sample == sample]['Allele']))
        for hla in hla_alleles:
            if hla == 'NPB':
                print_hla = hla
            else:
                print_hla = fxn.tidy_hla_name(hla)
            sample_hla_props.append([sample, print_hla, counts[hla]/sum(counts.values()), type_key[sample[0]]])

    sample_hla_props = fxn.list_to_df(sample_hla_props, ['Sample', 'Allele', 'Proportion', 'Type'], False)
    hla_prop = fxn.df_to_prop(sample_hla_predictions, 'Allele', 'Sample', 'Allele')
    hla_prop = hla_prop.sort_values(by='Sample')

    # HLA allele contribution barplot
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

    for a in list(set(sample_hla_props['Allele'])):
        allele_props = list(sample_hla_props.loc[(sample_hla_props.Allele == a)]['Proportion'])
        # print a, allele_props
        # print a, stats.mannwhitneyu(allele_props[:3], allele_props[3:])
        # # HLA-A*02:01 MannwhitneyuResult(statistic=0.0, pvalue=0.007124039836173583)
        # # HLA-B*07:02 MannwhitneyuResult(statistic=5.0, pvalue=0.05415969036500204)
        # # HLA-C*07:02 MannwhitneyuResult(statistic=5.0, pvalue=0.05415969036500204)
        # # NPB         MannwhitneyuResult(statistic=0.0, pvalue=0.007124039836173583)

    # HLA affinity distributions by sample
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Sample', y='Affinity', data=sample_hla_predictions, color='lightslategray', fliersize=0)
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    plt.ylim(2e0, 1e3)
    plt.yscale('log')
    plt.setp(ax.get_xticklabels(), rotation=85)
    plt.savefig(plot_dir + 'hla-affinity-distribution.png', dpi=300, bbox_inches='tight')

    # HLA affinity distributions by type
    fig = plt.figure(figsize=(2.5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Type', y='Affinity', data=sample_hla_predictions, color='lightslategray', fliersize=0)
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    plt.ylim(1e0, 1e3)
    plt.yscale('log')
    plt.savefig(plot_dir + 'hla-affinity-type-distribution.png', dpi=300, bbox_inches='tight')

    # Calculate whether there's any significant differences...
    affinities = coll.defaultdict(list)
    for sample in samples:
        affinities[sample] = list(sample_hla_predictions.loc[(sample_hla_predictions.Sample == sample), 'Affinity'])

    # for sample1 in samples:
    #     for sample2 in samples:
    #         print sample1, sample2
    #         print stats.mannwhitneyu(affinities[sample1], affinities[sample2])
    #         # ... there are some

    for i in range(8, 16):
        length_props = list(prop.loc[(prop.Length == i)]['Proportion'])
        # print i, stats.mannwhitneyu(length_props[:3], length_props[3:])
        # # 8    MannwhitneyuResult(statistic=0.0, pvalue=0.007124039836173583)
        # # 9    MannwhitneyuResult(statistic=8.0, pvalue=0.13594935540982417)
        # # 10   MannwhitneyuResult(statistic=11.0, pvalue=0.27705656503472276)
        # # 11   MannwhitneyuResult(statistic=8.0, pvalue=0.13594935540982417)
        # # 12   MannwhitneyuResult(statistic=8.0, pvalue=0.13594935540982417)
        # # 13   MannwhitneyuResult(statistic=6.0, pvalue=0.07539278251398564)
        # # 14   MannwhitneyuResult(statistic=5.0, pvalue=0.05415969036500204)
        # # 15   MannwhitneyuResult(statistic=0.0, pvalue=0.007124039836173583)

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
    plt.ylim(1e0, 1e5)
    for x in range(max(sharedness['Number-Samples'])):
        ax.text(x, 1e5, len(sharedness.loc[sharedness['Number-Samples'] == x+1]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    plt.savefig(plot_dir + 'hla-sharedness-affinity.png', dpi=300, bbox_inches='tight')

    # Calculate whether there's any significant differences...
    affinities = coll.defaultdict(list)
    for x in range(1, 6):
        affinities[x] = list(sharedness.loc[sharedness['Number-Samples'] == x, 'Affinity'])

    # for x in range(1, 14):
    #     for y in range(1, 14):
    #         mwu = stats.mannwhitneyu(affinities[x], affinities[y])
    #         if mwu[1] < 0.05:
    #             print x, y, mwu
    #             # Very many, as you'd expect - those found in just one sample are massively weaker affinity!

    # Need to plot affinity based on sharedness with (or not with) cultured sample
    type_sharedness = []
    affinities = coll.defaultdict(list)

    for pep in list(peptides.index):
        pep_vals = peptides.loc[pep]
        mouse_vals = list(pep_vals[3:])
        cultured_vals = list(pep_vals[:3])
        affinity = hla_predictions.loc[pep, 'Affinity']

        # Just in a single mouse or cultured sample
        if sum(mouse_vals) == 1 and sum(cultured_vals) == 0:
            type_sharedness.append(['1m-only', affinity])
            affinities['1m-only'].append(affinity)
        elif sum(mouse_vals) == 0 and sum(cultured_vals) == 1:
            type_sharedness.append(['1c-only', affinity])
            affinities['1c-only'].append(affinity)

        # Then those that are in 3 (or more) of one group, and absent in the other group
        elif sum(mouse_vals) >= 3 and sum(cultured_vals) == 0:
            type_sharedness.append(['3+m-0c', affinity])
            affinities['3+m-0c'].append(affinity)
        elif sum(mouse_vals) == 0 and sum(cultured_vals) == 3:
            type_sharedness.append(['3c-0m', affinity])
            affinities['3c-0m'].append(affinity)

        # Then those are in all samples
        elif sum(pep_vals) == len(pep_vals):
            type_sharedness.append(['All', affinity])
            affinities['All'].append(affinity)

    type_sharedness = fxn.list_to_df(type_sharedness, ['Samples-in', 'Affinity'], False)
    col_order = ['1c-only', '1m-only', '3c-0m', '3+m-0c', 'All']

    fig = plt.figure(figsize=(5.5, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Affinity', data=type_sharedness,
                fliersize=0, order=col_order, palette=fxn.c_m_mix_cols[0:2] * 2 + [fxn.c_m_mix_cols[2]])

    ax.set_yscale('log')
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    for x in range(len(col_order)):
        ax.text(x, 1e5, len(type_sharedness.loc[type_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    plt.ylim(1e0, 1e5)
    plt.xlabel("")
    plt.ylabel("MHCflurry predicted affinity (nM)")
    plt.savefig(plot_dir + 'hla-type-sharedness-affinity.png', dpi=300, bbox_inches='tight')
    plt.close()

    stats.mannwhitneyu(affinities['1c-only'], affinities['1m-only'])
    # MannwhitneyuResult(statistic=777084.0, pvalue=1.3073826188568937e-44)
    stats.mannwhitneyu(affinities['3+m-0c'], affinities['1c-only'])
    # MannwhitneyuResult(statistic=525258.5, pvalue=8.250331135935228e-31)
    stats.mannwhitneyu(affinities['1m-only'], affinities['3c-0m'])
    # MannwhitneyuResult(statistic=158230.0, pvalue=9.640479306316736e-66)
    stats.mannwhitneyu(affinities['3+m-0c'], affinities['3c-0m'])
    # MannwhitneyuResult(statistic=285626.0, pvalue=0.2778609006684889)     Only non signif!
    stats.mannwhitneyu(affinities['All'], affinities['3+m-0c'])
    # MannwhitneyuResult(statistic=1015015.0, pvalue=9.688772530997733e-15)
    stats.mannwhitneyu(affinities['All'], affinities['3c-0m'])
    # MannwhitneyuResult(statistic=98025.0, pvalue=0.00010784299757799932)

    # Output lists of group-specific peptides if no predicted binder
    C_only = [x for x in peptides.index if sum(list(peptides.loc[x])[:3]) > 0
              and sum(list(peptides.loc[x])[3:]) == 0 and hla_predictions.loc[x]['Allele'] == 'NPB']
    M_only = [x for x in peptides.index if sum(list(peptides.loc[x])[:3]) == 0
              and sum(list(peptides.loc[x])[3:]) > 0 and hla_predictions.loc[x]['Allele'] == 'NPB']
    B_only = [x for x in peptides.index if sum(list(peptides.loc[x])[:3]) > 0
              and sum(list(peptides.loc[x])[3:]) > 0 and hla_predictions.loc[x]['Allele'] == 'NPB']

    non_binder_dir = '../Data/NonPredictedBinders/'
    if not os.path.exists(non_binder_dir):
        os.mkdir(non_binder_dir)
    for t in ['C', 'M', 'B']:
        with open(non_binder_dir + type_key[t] + '-MHCflurry-predicted-non-binders.txt', 'w') as out_file:
            vars()[t + '_only'].sort()
            for p in vars()[t + '_only']:
                out_file.write(p + '\n')

    # See whether highly shared peptides derive from highly abundant RNA/proteins
    transcriptome, proteome = fxn.get_jy_ome_data()
    transcriptome_corr, proteome_corr = fxn.get_jy_ome_correlations(transcriptome, proteome, proteins)

    for d in list(proteins):

        # Transcriptomes
        tmp_df = transcriptome_corr[transcriptome_corr.Sample == d]
        fig = plt.figure(figsize=(4, 5))
        sns.jointplot(list(tmp_df['Peptidome Abundance']), list(tmp_df['Transcriptome Abundance']),
                      kind='reg', stat_func=stats.pearsonr)
        plt.xlabel("Immunopeptidome presence")
        plt.ylabel("RNA abundance log2(TPM)")
        plt.xlim(-.25, 1.25)
        plt.savefig(plot_dir + 'transcriptome-correlation-' + d + '.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Proteomes
        tmp_df = proteome_corr[proteome_corr.Sample == d]
        fig = plt.figure(figsize=(4, 5))
        sns.jointplot(list(tmp_df['Peptidome Abundance']), list(tmp_df['Proteome Abundance']),
                      kind='reg', stat_func=stats.pearsonr)
        plt.xlabel("Immunopeptidome presence")
        plt.ylabel("Protein abundance log2(intensity)")
        plt.xlim(-.25, 1.25)
        plt.savefig(plot_dir + 'proteome-correlation-' + d + '.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Need to plot affinity based on sharedness with (or not with) cultured sample
    prot_sharedness = []
    prot_abundances = coll.defaultdict(list)
    rna_sharedness = []
    rna_abundances = coll.defaultdict(list)

    for prot in list(proteins.index):

        prot_vals = proteins.loc[prot]
        mouse_vals = list(prot_vals[3:])
        cultured_vals = list(prot_vals[:3])

        # Just in a single mouse or cultured sample
        if sum(mouse_vals) == 1 and sum(cultured_vals) == 0:
            if prot in proteome:
                prot_sharedness.append(['1m-only', proteome[prot]])
                prot_abundances['1m-only'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['1m-only', transcriptome[prot]])
                rna_abundances['1m-only'].append(transcriptome[prot])

        elif sum(mouse_vals) == 0 and sum(cultured_vals) == 1:
            if prot in proteome:
                prot_sharedness.append(['1c-only', proteome[prot]])
                prot_abundances['1c-only'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['1c-only', transcriptome[prot]])
                rna_abundances['1c-only'].append(transcriptome[prot])

        # Then those that are in 3 (or more) of one group, and absent in the other group
        elif sum(mouse_vals) >= 3 and sum(cultured_vals) == 0:
            if prot in proteome:
                prot_sharedness.append(['3+m-0c', proteome[prot]])
                prot_abundances['3+m-0c'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['3+m-0c', transcriptome[prot]])
                rna_abundances['3+m-0c'].append(transcriptome[prot])

        elif sum(mouse_vals) == 0 and sum(cultured_vals) == 3:
            if prot in proteome:
                prot_sharedness.append(['3c-0m', proteome[prot]])
                prot_abundances['3c-0m'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['3c-0m', transcriptome[prot]])
                rna_abundances['3c-0m'].append(transcriptome[prot])

        # Then those are in all samples
        elif sum(prot_vals) == len(prot_vals):
            if prot in proteome:
                prot_sharedness.append(['All', proteome[prot]])
                prot_abundances['All'].append(proteome[prot])
            if prot in transcriptome:
                rna_sharedness.append(['All', transcriptome[prot]])
                rna_abundances['All'].append(transcriptome[prot])

    prot_sharedness = fxn.list_to_df(prot_sharedness, ['Samples-in', 'Proteome Abundance'], False)
    rna_sharedness = fxn.list_to_df(rna_sharedness, ['Samples-in', 'Transcript Abundance'], False)
    col_order = ['1c-only', '1m-only', '3c-0m', '3+m-0c', 'All']

    stats.mannwhitneyu(prot_abundances['1m-only'], prot_abundances['1c-only'])
    # MannwhitneyuResult(statistic=4554.0, pvalue=0.48700010373996278)
    stats.mannwhitneyu(prot_abundances['3+m-0c'], prot_abundances['3c-0m'])
    # MannwhitneyuResult(statistic=7731.5, pvalue=0.16108702752604553)
    stats.mannwhitneyu(prot_abundances['3+m-0c'], prot_abundances['All'])
    # MannwhitneyuResult(statistic=97293.0, pvalue=4.5065510948557757e-10)
    stats.mannwhitneyu(prot_abundances['3c-0m'], prot_abundances['All'])
    # MannwhitneyuResult(statistic=10078.5, pvalue=0.065018001099927808)
    stats.mannwhitneyu(prot_abundances['1m-only'], prot_abundances['All'])
    # MannwhitneyuResult(statistic=33759.5, pvalue=2.8845626328782545e-05)
    stats.mannwhitneyu(prot_abundances['1c-only'], prot_abundances['All'])
    # MannwhitneyuResult(statistic=14324.5, pvalue=0.0014614951305713328)

    stats.mannwhitneyu(rna_abundances['1m-only'], rna_abundances['1c-only'])
    # MannwhitneyuResult(statistic=28644.0, pvalue=0.1939781088731688)
    stats.mannwhitneyu(rna_abundances['3+m-0c'], rna_abundances['3c-0m'])
    # MannwhitneyuResult(statistic=25491.0, pvalue=0.0020202471186718065)
    stats.mannwhitneyu(rna_abundances['3+m-0c'], rna_abundances['All'])
    # MannwhitneyuResult(statistic=254593.0, pvalue=1.5494157262165817e-51)
    stats.mannwhitneyu(rna_abundances['3c-0m'], rna_abundances['All'])
    # MannwhitneyuResult(statistic=24795.0, pvalue=0.011459113063982137)
    stats.mannwhitneyu(rna_abundances['1m-only'], rna_abundances['All'])
    # MannwhitneyuResult(statistic=93436.0, pvalue=4.738827709249312e-39)
    stats.mannwhitneyu(rna_abundances['1c-only'], rna_abundances['All'])
    # MannwhitneyuResult(statistic=43263.0, pvalue=7.604260506183144e-14)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Transcript Abundance', data=rna_sharedness, order=col_order,
                fliersize=0, palette=fxn.c_m_mix_cols[0:2] * 2 + [fxn.c_m_mix_cols[2]])
    for x in range(len(col_order)):
        ax.text(x, 10, len(rna_sharedness.loc[rna_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='RNA abundance log2(TPM)')
    plt.ylim(-2, 10)
    plt.savefig(plot_dir + 'transcriptome-abundance-by-sharedness.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Proteome Abundance', data=prot_sharedness, order=col_order,
                fliersize=0, palette=fxn.c_m_mix_cols[0:2] * 2 + [fxn.c_m_mix_cols[2]])
    for x in range(len(col_order)):
        ax.text(x, 38, len(prot_sharedness.loc[prot_sharedness['Samples-in'] == col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='Samples present in', ylabel='Protein abundance log2(intensity)')
    plt.ylim(21, 38)
    plt.savefig(plot_dir + 'proteome-abundance-by-sharedness.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Let's plot small versions for concise visualisation
    short_col_order = ['1c-only', '1m-only', 'All']

    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Transcript Abundance', data=rna_sharedness, order=short_col_order,
                fliersize=0, palette=fxn.c_m_mix_cols)
    for x in range(len(short_col_order)):
        ax.text(x, 10, len(rna_sharedness.loc[rna_sharedness['Samples-in'] == short_col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='', ylabel='RNA abundance log2(TPM)')
    plt.ylim(-2, 10)
    plt.savefig(plot_dir + 'transcriptome-abundance-by-sharedness-small.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Samples-in', y='Proteome Abundance', data=prot_sharedness, order=short_col_order,
                fliersize=0, palette=fxn.c_m_mix_cols)
    for x in range(len(short_col_order)):
        ax.text(x, 38, len(prot_sharedness.loc[prot_sharedness['Samples-in'] == short_col_order[x]]),
                horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='', ylabel='Protein abundance log2(intensity)')
    plt.ylim(21, 38)
    plt.savefig(plot_dir + 'proteome-abundance-by-sharedness-small.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Output lists of all proteins that are featured in ALL mice samples or ALL cultured samples, for ORA
    ora = coll.defaultdict(list)
    for p in list(proteins.index):
        vals = list(proteins.loc[p])
        if sum(vals[:3]) == 3:
            ora['in_all_cultured'].append(p)
        if sum(vals[3:]) == 10:
            ora['in_all_mice'].append(p)

    fxn.save_venn2(plot_dir + 'ORA-lists', ora, 'in_all_cultured', 'in_all_mice', 'Cultured', 'Mice')
    ora_dir = '../Data/ORA/'
    for growth_type in ['mice', 'cultured']:
        out_file_path = ora_dir + 'proteins-inall' + growth_type + '-peptides.txt'
        with open(out_file_path, 'w') as out_file:
            for p in ora['in_all_' + growth_type]:
                out_file.write(p + '\n')

    # Then run these through WebGestalt's ORA process; can then unpack/read in the results and plot the correlation
    ora_results_dir = ora_dir + 'Results/'
    ora_go_p = coll.defaultdict(fxn.nest_counter)
    for growth_type in ['mice', 'cultured']:
        in_dat_dir = ora_results_dir + 'peptide-' + growth_type + '/'
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

    fxn.save_venn2(plot_dir + 'ORA-GO-lists', ora_go_p, 'cultured', 'mice', 'Cultured', 'Mice')
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
                  y='log10 mice sample GO term P val', kind='reg', stat_func=stats.pearsonr, color=fxn.brown)
    plt.savefig(plot_dir + 'ORA-1minus-correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
