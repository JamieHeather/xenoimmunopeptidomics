# -*- coding: utf-8 -*-

"""
mouse-proteome-check.py

Look to see whether non-predicted HLA binding peptides are more likely to get a match to the murine proteome

"""

from __future__ import division
import functions as fxn
import os
import gzip
from acora import AcoraBuilder
import collections as coll
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

__version__ = '0.3.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


if __name__ == '__main__':

    fxn.check_scripts_dir()
    sns.set(font="Arial", font_scale=1.5)

    # Sort directories, get data
    plot_dir = fxn.plot_dir + fxn.get_date() + '-mouse-proteome-check/'
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # Read proteome into dict
    mouse_proteome_file = [x for x in os.listdir(fxn.base_data_dir) if '_mouse.fasta' in x][0]

    mouse_proteins = coll.defaultdict()
    with gzip.open(fxn.base_data_dir + mouse_proteome_file, 'rU') as in_file:
        for protein, seq, blank in fxn.read_fa(in_file):
            mouse_proteins[protein.split(' ')[0]] = seq

    # Then scroll through non-predicted binder files, build an AC trie of all the peptides per file
    data_dir = '../Data/NonPredictedBinders/'
    matches = coll.defaultdict(fxn.nest_counter)
    all_peptides = coll.defaultdict(list)
    for f in [x for x in os.listdir(data_dir) if x.endswith('.txt')]:
        nam = f.split('-')[0]
        search_builder = AcoraBuilder()
        peptides = []

        # Build trie
        with open(data_dir + f, 'rU') as in_file:
            for line in in_file:
                search_builder.add(line.rstrip())
                peptides.append(line.rstrip())
                all_peptides[f.split('-')[0]].append(line.rstrip())
        seq_search = search_builder.build()

        # Use to search all proteins in proteome
        for protein in mouse_proteins:
            seq_check = seq_search.findall(mouse_proteins[protein])
            if seq_check:
                for s in seq_check:
                    matches[nam][s[0]] += 1

        # Then fill in the zeroes (unmatched peptides) to get denominator
        for p in peptides:
            if p not in matches[nam]:
                matches[nam][p] = 0

    # Go through the matches and count!
    vals = []
    positives = {}
    totals = {}
    col_order = ['Cultured', 'Mouse', 'Both']
    for growth_type in col_order:
        positive = len([x for x in matches[growth_type] if matches[growth_type][x] > 0])
        total = len(matches[growth_type])
        print growth_type, positive, total, positive/total
        vals.append([growth_type, positive/total])
        positives[growth_type] = positive
        totals[growth_type] = total

    vals = fxn.list_to_df(vals, ['Type', 'Proportion'], False)

    # Finally plot
    fig = plt.figure(figsize=(4.5, 5))
    ax = fig.add_subplot(111)
    sns.barplot(data=vals, x='Type', y='Proportion', palette=fxn.c_m_mix_cols, order=col_order)
    for x in range(len(col_order)):
        plot_text = str(positives[col_order[x]]) + '/\n' + str(totals[col_order[x]])
        print plot_text
        ax.text(x, .64, plot_text, horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.set(xlabel='', ylabel='Proportion murine-matching peptides')
    plt.savefig(plot_dir + 'NPB-mouse-proteome-matches.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Also want to check to see whether there's any difference in the predicted ability of these peptides
    # to bind to any of the class I alleles present in NSG mice (H2d haplotype)
    fxn.initialise_prediction()
    hla_predictions = []
    for found_in in col_order:
        for pep in all_peptides[found_in]:
            affin = 1e6
            allele = ''
            for mouse_allele in ['H-2-Dd', 'H-2-Kd', 'H-2-Ld']:
                prediction = fxn.hla_prediction(pep, mouse_allele)[0]
                if prediction < affin:
                    affin = prediction
                    allele = mouse_allele
            if affin <= 50:
                assignment = '+++'
            elif affin <= 500:
                assignment = '++'
            elif affin <= 1000:
                assignment = '+'
            else:
                assignment = '-'
                allele = 'NPB'  # No predicted binder
            hla_predictions.append([pep, allele, found_in, affin, assignment])

    hla_predictions = fxn.list_to_df(hla_predictions, ['Peptide', 'Allele', 'Sample', 'Affinity', 'Assignment'], True)

    # Plot
    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    sns.boxplot(x='Sample', y='Affinity', data=hla_predictions, palette=fxn.c_m_mix_cols, order=col_order, fliersize=0)
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    # plt.ylim(1e0, 1e3)
    plt.yscale('log')
    plt.xlabel('')
    plt.ylabel("MHCflurry predicted affinity (nM)")
    plt.savefig(plot_dir + 'h2d-affinity-type-distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Let's see which human proteins the peptides map to - different here, as we care which protein it is
    human_proteome_file = [x for x in os.listdir(fxn.base_data_dir) if '_human.fasta' in x][0]

    human_proteins = coll.defaultdict(list)
    with gzip.open(fxn.base_data_dir + human_proteome_file, 'rU') as in_file:
        for protein, seq, blank in fxn.read_fa(in_file):
            # if 'GN=' in protein:
            if not protein.startswith('tr|'):
                if 'GN=' in protein:
                    gene = [x[3:] for x in protein.split(' ') if 'GN=' in x][0]
                else:
                    gene = protein.split(' ')[0].split('|')[2].replace('_HUMAN', '')
                human_proteins[gene].append(seq)
            else:  # Ignore TrEMBL
                continue

    # Then scroll through non-predicted binder files, build an AC trie of all the peptides per file
    data_dir = '../Data/NonPredictedBinders/'
    ora_dir = data_dir + 'ORA/'
    if not os.path.exists(ora_dir):
        os.mkdir(ora_dir)

    matches = coll.defaultdict(fxn.nest)
    all_peptides = coll.defaultdict(list)
    single_hit_prots = coll.defaultdict(list)
    for f in [x for x in os.listdir(data_dir) if x.endswith('.txt')]:
        nam = f.split('-')[0]
        search_builder = AcoraBuilder()
        peptides = []
        #
        # Build trie
        with open(data_dir + f, 'rU') as in_file:
            for line in in_file:
                search_builder.add(line.rstrip())
                peptides.append(line.rstrip())
                all_peptides[f.split('-')[0]].append(line.rstrip())
        seq_search = search_builder.build()
        #
        # Use to search all proteins in proteome
        for gene in human_proteins:
            for protein in human_proteins[gene]:
                seq_check = seq_search.findall(protein)
                if seq_check:
                    for s in seq_check:
                        if gene not in matches[nam][s[0]]:  # Only count one protein per gene symbol
                            matches[nam][s[0]].append(gene)

        # Output text files of peptides which match to single proteins
        single_hits = [x for x in matches[nam] if len(matches[nam][x]) == 1]
        single_hit_prots[nam] = [matches[nam][x][0] for x in single_hits]
        with open(ora_dir + 'single-Hs-hits-' + nam + '.txt', 'w') as out_file:
            for peptide in single_hits:
                out_file.write(matches[nam][peptide][0] + '\n')

    # Then run these through WebGestalt's ORA process; can then unpack/read in the results and plot the correlation
    ora_results_dir = ora_dir + 'Results/'
    ora_go_p = coll.defaultdict(fxn.nest_counter)
    for growth_type in ['mice', 'cultured']:
        in_dat_dir = ora_results_dir + growth_type + '/'
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
    fxn.save_venn2(plot_dir + 'ORA-prot-lists', single_hit_prots, 'Cultured', 'Mouse', 'Cultured', 'Mice')

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

    # Test if the enriched proteins correlate between NPBs and the conserved/predicted to bind
    ora_dir2 = '../Data/ORA/'
    ora_results_dir2 = ora_dir2 + 'Results/'
    # ora_go_p = coll.defaultdict(fxn.nest_counter)
    for growth_type in ['mice', 'cultured']:
        in_dat_dir = ora_results_dir2 + 'peptide-' + growth_type + '/'
        enrichment_file_path = [x for x in os.listdir(in_dat_dir) if x.startswith('enrichment')][0]
        with open(in_dat_dir + enrichment_file_path, 'rU') as in_file:
            line_count = 0
            for line in in_file:
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    headers = bits
                else:
                    ora_go_p['gp-' + growth_type][bits[0]] = float(bits[7])
                line_count += 1

    fxn.save_venn2(plot_dir + 'ORA-GO-NPB', ora_go_p, 'cultured', 'gp-cultured', 'NPB (any)', 'PB (all)')
    fxn.save_venn2(plot_dir + 'ORA-GO-NPB', ora_go_p, 'mice', 'gp-mice', 'NPB (any)', 'PB (all)')

    all_go = list(set(
        ora_go_p['mice']) | set(ora_go_p['cultured']) | set(ora_go_p['gp-mice']) | set(ora_go_p['gp-cultured']))

    go_df = []
    for g in all_go:
        for growth_type in ['mice', 'cultured', 'gp-mice', 'gp-cultured']:
            # In order to get sensible logged plotting, make non-present GO term Pvals = 1,
            # and '0' vals rounded down to order of magnitude just below lowest
            if g not in ora_go_p[growth_type]:
                ora_go_p[growth_type][g] = 1
            elif ora_go_p[growth_type][g] == 0:
                ora_go_p[growth_type][g] = 1e-16
        go_df.append([g, np.log10(ora_go_p['cultured'][g]), np.log10(ora_go_p['mice'][g]),
                      np.log10(ora_go_p['gp-cultured'][g]), np.log10(ora_go_p['gp-mice'][g])])

    go_df = fxn.list_to_df(go_df, ['GO', 'log10 NPB cultured (any) GO term P val',
                                   'log10 NPB mice (any) sample GO term P val',
                                   'log10 PB cultured (all) GO term P val',
                                   'log10 PB mice (all) sample GO term P val'], True)

    done = []
    count = 0
    for x_dat in go_df[1:]:
        for y_dat in go_df[1:]:
            if x_dat != y_dat:
                pair = [x_dat, y_dat]
                pair.sort()
                if pair not in done:
                    fig = plt.figure(figsize=(5.5, 5))
                    ax = fig.add_subplot(111)
                    sns.jointplot(data=go_df, x=x_dat, y=y_dat, kind='reg', stat_func=stats.pearsonr, color=fxn.brown)
                    plt.savefig(plot_dir + 'ORA-1minus-corr-' + str(count) + '.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    count += 1
                    done.append(pair)
