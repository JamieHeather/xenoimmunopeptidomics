# -*- coding: utf-8 -*-

"""
published-comparison.py

"""

from __future__ import division
import os
import functions as fxn
import matplotlib.pyplot as plt
import matplotlib_venn as venn
import collections as coll
import itertools as it
import numpy as np
import seaborn as sns
from scipy import stats


__version__ = '0.3.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def get_data_source(sample_name):
    if sample_name.startswith('M') and len(sample_name) == 3:
        return 'mouse'
    elif sample_name.startswith('C0'):
        return 'culture'
    elif 'THP1' in sample_name or 'MM' in sample_name or 'C1866' in sample_name:
        return 'control'
    else:
        return 'published'


if __name__ == '__main__':

    plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})
    sns.set(font="Arial", font_scale=1.5)

    plot_dir = fxn.plot_dir + fxn.get_date() + '-publishedJY-analysis/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # Get data processed by generate-data-files.py
    peptides, pep_types = fxn.get_data('peptide')

    pep_samples = peptides.keys()
    pep_samples.sort()
    # Specify explicitly which are control lines, for ease of plotting
    control_data = [x for x in pep_samples if 'Jurkat' in x or 'C1866' in x or 'THP' in x or 'MM16' in x]
    mouse_data = [x for x in pep_samples if x.startswith('M')]
    culture_data = [x for x in pep_samples if x.startswith('C0')]
    published_data = [x for x in pep_samples if x not in control_data and x not in mouse_data and x not in culture_data]
    pep_samples = culture_data + mouse_data + published_data + control_data

    # Define shorter IDs for each dataset which can be used for plotting/saving labels
    ids = {'BS2015': 'Bassani',
           'Bourdetsky2014': 'Bourdetsky', 'Hassan2013': 'Hassan',
           'M01': 'M01', 'M02': 'M02', 'M04': 'M04',
           'M05': 'M05', 'M06': 'M06', 'M07': 'M07', 'M08': 'M08',
           'M09': 'M09', 'M10': 'M10', 'M11': 'M11',
           'C01': 'C01', 'C02': 'C02', 'C03': 'C03',
           'BSSysMHCTHP1v1': 'THP1',
           'Caron2015Jurkat': 'Jurkat', 'Caron2015JY': 'Caron',
           'TernetteSysMHCC1866': 'C1866', 'BSSysMHCMM16': 'MM16'}

    pep_ids = [ids[x] for x in pep_samples]

    # For a given number of repeats, go through each dataset and sample to a given number,
    # calculate Jaccard/overlap and then plot an averaged heat map
    repeats = 100
    sample_to = 1500
    sampled_jacc = coll.defaultdict(fxn.nest)

    for i in range(repeats):

        print str(i) + '\t-----------------------------------------------------'
        # Subsample each dataset to a fixed number
        dat = {}
        for dataset in peptides:
            print '\t' + dataset
            sample = fxn.subsample_to_number(peptides, dataset, sample_to)
            dat[dataset] = sample

        for x in range(len(pep_samples)):
            for y in range(len(pep_samples)):
                sampled_jacc[x][y].append(fxn.jaccard(dat[pep_samples[x]], dat[pep_samples[y]]))

    plot_arr = np.zeros(shape=(len(pep_samples), len(pep_samples)))
    for x in range(len(pep_samples)):
        for y in range(len(pep_samples)):
            plot_arr[x, y] = np.mean(sampled_jacc[x][y])

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    p = ax.pcolor(plot_arr, cmap='gnuplot', vmin=0, vmax=1)
    ax.set_xticks([x + .5 for x in range(len(pep_samples))])
    ax.set_xticklabels(pep_ids, rotation=85)
    ax.set_yticks([x + .5 for x in range(len(pep_samples))])
    ax.set_yticklabels(pep_ids)  # , fontsize=4)
    plt.xlim(0, len(pep_samples))
    plt.ylim(0, len(pep_samples))
    cbar = plt.colorbar(p)
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel(u'← less overlap         Jaccard index         more overlap →', rotation=90)
    plt.savefig(plot_dir + 'Subsampled-Jaccards-' + str(repeats) + 'repeats-' + str(sample_to) + 'sampled.png',
                dpi=300, bbox_inches='tight')
    plt.close()

    plot_arr = np.zeros(shape=(len(pep_samples), len(pep_samples)))
    for x in range(len(pep_samples)):
        for y in range(len(pep_samples)):
            plot_arr[x, y] = fxn.jaccard(peptides[pep_samples[x]], peptides[pep_samples[y]])

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    p = ax.pcolor(plot_arr, cmap='gnuplot', vmin=0, vmax=1)
    ax.set_xticks([x + .5 for x in range(len(pep_samples))])
    ax.set_xticklabels(pep_ids, rotation=85)
    ax.set_yticks([x + .5 for x in range(len(pep_samples))])
    ax.set_yticklabels(pep_ids)  # , fontsize=4)
    plt.xlim(0, len(pep_samples))
    plt.ylim(0, len(pep_samples))
    cbar = plt.colorbar(p)
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel(u'← less overlap         Jaccard index         more overlap →', rotation=90)
    plt.savefig(plot_dir + 'WholeSampleJaccards.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plotting a violin plot of the Jaccards, so I need only the important comparisons (and each only once)
    sampled_jacc_long = []
    whole_jacc_long = []
    done_comparisons = []
    for x in range(len(pep_samples)):
        for y in range(len(pep_samples)):
            if x != y and str(x) + '-' + str(y) not in done_comparisons:
                x_type = get_data_source(pep_samples[x])
                y_type = get_data_source(pep_samples[y])
                comparison_type = [x_type, y_type]
                comparison_type.sort()
                if comparison_type != ['control', 'control']:
                    for jacc in sampled_jacc[x][y]:
                        sampled_jacc_long.append(['-\n'.join(comparison_type), jacc])
                    whole_jacc_long.append(['-\n'.join(comparison_type),
                                            fxn.jaccard(peptides[pep_samples[x]], peptides[pep_samples[y]])])
                done_comparisons.append(str(x) + '-' + str(y))
                done_comparisons.append(str(y) + '-' + str(x))

    sampled_jacc_long = fxn.list_to_df(sampled_jacc_long, ['Comparison', 'Jaccard'], False)
    whole_jacc_long = fxn.list_to_df(whole_jacc_long, ['Comparison', 'Jaccard'], False)

    col_order = ['control-\npublished',
                 'control-\nculture',
                 'control-\nmouse',
                 'published-\npublished',
                 'culture-\npublished',
                 'mouse-\npublished',
                 'culture-\nmouse',
                 'mouse-\nmouse',
                 'culture-\nculture']

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    sns.violinplot(x='Comparison', y='Jaccard', data=sampled_jacc_long, order=col_order,
                   cut=0, color='lightslategray')
    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.savefig(plot_dir + 'sub-Jaccard-violin.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    sns.violinplot(x='Comparison', y='Jaccard', data=whole_jacc_long, order=col_order,
                   cut=0, color='lightslategray')
    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.savefig(plot_dir + 'whole-Jaccard-violin.png', dpi=300, bbox_inches='tight')
    plt.close()

    # And calc stats
    p_p = list(sampled_jacc_long.loc[(sampled_jacc_long.Comparison == 'published-\npublished')]['Jaccard'])
    m_p = list(sampled_jacc_long.loc[(sampled_jacc_long.Comparison == 'mouse-\npublished')]['Jaccard'])
    c_p = list(sampled_jacc_long.loc[(sampled_jacc_long.Comparison == 'culture-\npublished')]['Jaccard'])

    stats.mannwhitneyu(p_p, m_p)
    stats.mannwhitneyu(p_p, c_p)
    stats.mannwhitneyu(m_p, c_p)

    p_p = list(whole_jacc_long.loc[(whole_jacc_long.Comparison == 'published-\npublished')]['Jaccard'])
    m_p = list(whole_jacc_long.loc[(whole_jacc_long.Comparison == 'mouse-\npublished')]['Jaccard'])
    c_p = list(whole_jacc_long.loc[(whole_jacc_long.Comparison == 'culturfe-\npublished')]['Jaccard'])

    stats.mannwhitneyu(p_p, m_p)
    stats.mannwhitneyu(p_p, c_p)
    stats.mannwhitneyu(m_p, c_p)
