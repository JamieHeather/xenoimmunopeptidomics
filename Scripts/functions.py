# -*- coding: utf-8 -*-

"""
functions.py

Functions for the Xeno-Immunopeptidomics project

"""

from __future__ import division
import json
import os
import sys
import random
import collections as coll
import matplotlib.pyplot as plt
import matplotlib_venn as venn
import datetime
import operator
import numpy as np
import pandas as pd
from time import strftime

__version__ = '0.7.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})

base_data_dir = '../Data/'
plot_dir = '../Plots/'

# Pulled colours out from sns pallete directly, via sns.color_palette()[.as_hex()]
blue = (0.2980392156862745, 0.4470588235294118, 0.6901960784313725)
orange = (0.8666666666666667, 0.5176470588235295, 0.3215686274509804)
purple = (0.5058823529411764, 0.4470588235294118, 0.7019607843137254)
brown = (0.5764705882352941, 0.47058823529411764, 0.3764705882352941)
c_m_mix_cols = [u'#4c72b0', u'#dd8452', u'#937860']  # Cultured/mouse/mixed colour hexes

proteome_key = {'UP000005640': 'human', 'UP000000589': 'mouse',
                'mouse': 'UP000000589', 'human': 'UP000005640'}


def check_scripts_dir():
    """
    Check we're in the right directory (Scripts)
    :return:
    """
    if not os.getcwd().endswith('/Scripts'):
        if 'Scripts' in os.listdir(os.getcwd()):
            os.chdir('Scripts')
        else:
            print "Check your current working directory - this is designed to be run from root or Scripts folders"
            sys.exit()


def get_date():
    """
    :return: the current date in ISO 8601 format
    """
    return strftime("%Y-%m-%d")


def get_data(data_type):
    """
    Read in data to nested dictionaries, such that there's a top level dict containing counters or lists of peptides
    :param data_type: type of ligandome data, i.e. 'peptide' or 'phospho' [-peptide]
    :return: nested dict containing each file's data type, and a separate dict detailing data types
    """
    all_files = [x for x in os.listdir(base_data_dir + data_type) if not x.startswith('.') and
                 (x.endswith('tsv') or x.endswith('txt'))]
    all_files.sort()
    print all_files

    dat = {}
    abundances = coll.defaultdict(list)

    for f in all_files:
        if '-confidential' in f:
            nam = f[:-17]
        else:
            nam = f[:-4]

        if f.endswith('.tsv'):
            spec_dat = coll.Counter()
            abundances['counts'].append(nam)

            with open(base_data_dir + data_type + '/' + f, 'rU') as in_file:
                for line in in_file:
                    bits = line.rstrip().split('\t')
                    try:
                        spec_dat[bits[0]] = int(bits[1])
                    except:
                        spec_dat[bits[0]] = float(bits[1])

        elif f.endswith('.txt'):
            spec_dat = []
            abundances['list'].append(nam)

            with open(base_data_dir + data_type + '/' + f, 'rU') as in_file:
                for line in in_file:
                    spec_dat.append(line.rstrip())

        else:
            print "Error - incorrect file type detected:", f
            sys.exit()

        dat[nam] = spec_dat

    return dat, abundances


def jaccard(list1, list2):
    """
    Return Jaccard index of two lists
    """
    intersection = list(set(list1) & set(list2))
    union = list(set(list1) | set(list2))
    return len(intersection) / len(union)


def proportional(savename, dat, name1, name2, label1, label2):
    """
    Plot a scatter plot of two lists of abundances.
    :param savename: Unique name to append to file names
    :param dat: nested dictionary containing each samples' data
    :param name1: name of first samples' data in nested dict
    :param name2: name of first samples' data in nested dict
    :param label1: label to use when plotting
    :param label2: label to use when plotting
    :return: Just plots
    """
    counter1 = dat[name1]
    counter2 = dat[name2]

    sum1 = sum(counter1.values())
    sum2 = sum(counter2.values())

    xs = []
    ys = []

    for pep in counter1:
        xs.append(counter1[pep] / sum1)
        if pep in counter2:
            ys.append(counter2[pep] / sum2)
        else:
            ys.append(0)

    for pep in [p for p in counter2.keys() if p not in counter1.keys()]:
        xs.append(0)
        ys.append(counter2[pep] / sum2)

    x_fraction = max(xs) / 100
    y_fraction = max(ys) / 100

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.scatter(xs, ys)
    plt.xlabel(label1 + ' proportions')
    plt.ylabel(label2 + ' proportions')
    plt.xlim(0 - x_fraction, max(xs) + x_fraction)
    plt.ylim(0 - y_fraction, max(ys) + y_fraction)
    plt.tight_layout()

    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d")
    outname = plot_dir + date + '-' + savename + '-Proportions-' + name1 + '-' + name2 + '.png'
    plt.savefig(outname, dpi=300)
    plt.close()
    return


def save_venn2(plot_dir, dat, name1, name2, label1, label2):
    """
    Plot a proportional 2-circle Venn diagram of two lists. Also include Jaccard index on plot
    :param plot_dir: directory to plot into
    :param dat: nested dictionary containing each samples' data
    :param name1: name of first samples' data in nested dict
    :param name2: name of first samples' data in nested dict
    :param label1: label to use when plotting
    :param label2: label to use when plotting
    :return: Jaccard index of the pairing
    """
    counter1 = dat[name1]
    counter2 = dat[name2]

    jac = round(jaccard(dat[name1], dat[name2]), 3)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    v = venn.venn2([set(counter1), set(counter2)], set_labels=(label1, label2))
    v.get_patch_by_id('10').set_color(blue)
    v.get_patch_by_id('01').set_color(orange)

    v.get_patch_by_id('10').set_edgecolor('none')
    v.get_patch_by_id('01').set_edgecolor('none')

    v.get_patch_by_id('11').set_color(brown)
    v.get_patch_by_id('11').set_edgecolor('none')

    outname = plot_dir + '-Venn-' + name1 + '-' + name2 + '.png'
    plt.title('Jaccard index: ' + str(jac))
    plt.savefig(outname, dpi=300)
    plt.close()

    return jac


def get_relations(savename, dat, name1, name2, label1, label2):
    """
    Given a data dict and a set of experiments to compare, run through each of the individual comparison functions
    :param savename: A general ID included in save names, e.g. 'peptide' or 'phospho'
    :param dat: The actual dictionary read in by get_data
    :param name1: Name of the nested dictionary for first sample
    :param name2: Name of the nested dictionary for second sample
    :param label1: Plot label for first sample
    :param label2: Plot label for second sample
    :return: Saves the plots in the relevant dir as it goes
    """

    if isinstance(dat[name1], dict) and isinstance(dat[name2], dict):
        proportional(savename, dat, name1, name2, label1, label2)

    jac = save_venn2(savename, dat, name1, name2, label1, label2)

    return jac


def subsample_to_number(dat, nam, number):
    """
    Given a dictionary or list, randomly subsample it down to a given number.
    Note that for input data with abundance values this requires an iterative process,
      ... taking larger and larger subsamples until the unique number of species hits the desired number
    :param dat:
    :param nam:
    :param number:
    :return:
    """

    # If the input data is a simiple list, just simple subsample the desired number directly
    if isinstance(dat[nam], list):
        return coll.Counter(random.sample(dat[nam], number))

    # Otherwise if it's a dictionary, begin the iterative subsampling process
    elif isinstance(dat[nam], dict):

        peptides = [x for x in dat[nam].keys()]
        probabilities = [dat[nam][x] / sum(dat[nam].values()) for x in dat[nam]]

        counter = coll.Counter()
        out_number = number

        # Subsample from the desired number, increasing the sample size until unique number of peptides desired matched
        while len(counter) <= number:

            # When subsampling at a given number produces a unique species count very far from the target, don't linger
            # However when it gets close it then has multiple attempts per subsampling depth
            if abs(len(counter) - number) > 10:
                attempts = 1
            else:
                attempts = 5

            while attempts != 0:

                sample = np.random.choice(peptides, out_number, p=probabilities)
                counter = coll.Counter(sample)
                attempts -= 1

                if len(counter) == number:
                    return counter
                elif len(counter) > number:
                    attempts += 2

            out_number += 1

        return "Error"

    else:
        print "Type error: " + nam
        sys.exit()


def nest():
    """
    Create nested defaultdicts
    """
    return coll.defaultdict(list)


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def hla_prediction(peptide, hla_allele):
    """
    Actual code that performs the mhcflurry prediction. initialise_prediction must be run first
    :param peptide: peptide sequence to be tested
    :param hla_allele: HLA allele for prediction (NB: rare alleles may not be covered, and throw an error!)
    :return: the predicted (average) Kd (nM)
    """
    return predictor.predict_to_dataframe(peptides=[peptide], allele=hla_allele)['prediction'].iloc


def initialise_prediction():
    """
    Enable mhcflurry prediction
    :return:
    """
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Prevents warnings. No, I'm not building from source
    from mhcflurry import Class1AffinityPredictor
    global predictor
    predictor = Class1AffinityPredictor.load()


def generate_upset_sets(in_data):
    """
    From/for upsetplot - see https://github.com/jvivian/UpSetPlot/commit/b6dea63b3f2b4f61f54fc450fe6368a01bb556cc
    Data loader for a dict of sets
    :param in_data:
    :return: Multi-Index Series of intersection counts
    """
    # Convert dict into OrderedDict to preserve Key/Value order
    data_dict = coll.OrderedDict(in_data)
    # Construct Index
    tf_array = [[True, False]] * len(data_dict)
    index = pd.MultiIndex.from_product(tf_array, names=data_dict.keys())
    # Curate values from each intersection group
    values = []
    for i in index:
        values.append(_intersection_counts(data_dict.values(), i))
    return pd.Series(values, index=index)


def _intersection_counts(sets, bool_tuple):
    """
    From upsetplot - see https://github.com/jvivian/UpSetPlot/commit/b6dea63b3f2b4f61f54fc450fe6368a01bb556cc
    Given list of sets and boolean tuple, return count of intersection
    :param List[sets[str]] sets:
    :param Tuple[bool] bool_tuple:
    :return: Count of intersection
    """

    # For all False case, return 0
    if True not in bool_tuple:
        return 0
    # Operator dictionary
    set_ops = {True: operator.and_, False: operator.sub}
    # For each grouping, perform set operation
    zipped = sorted(list(zip(bool_tuple, sets)), reverse=True)
    _, base = zipped[0]
    for operation, s in zipped[1:]:
        base = set_ops[operation](base, s)
    return len(base)


def list_to_df(input_list, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param input_list: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(input_list)
    df = df.rename(index=str, columns=dict(zip(range(len(headers)), headers)))
    df = df.sort_values(by=[headers[0]])
    if rename is True:
        df = df.set_index(headers[0], drop=True)
    return df


def df_to_prop(df, value, group_by, sort_by):
    """
    Convert a pandas dataframe of data to a proportions
    :param df: summary pandas dataframe with absolute measures
    :param value: column of the dataframe containing the data to convert to proportions
    :param group_by: column of the dataframe containing the groups to summarise under (e.g. 'Sample')
    :param sort_by: column to sort on
    :return: dataframe with normalized proportion values
    """
    df_prop = (df[value]
               .groupby(df[group_by])
               .value_counts(normalize=True)
               .rename('Proportion')
               .reset_index())
    df_prop = df_prop.sort_values(by=sort_by)
    return df_prop


def read_data_in(file_path):
    """
    :param file_path: Path to peptide or protein list
    :return: dictionary of same data
    """

    if file_path.endswith('.tsv'):
        dat = coll.Counter()
        with open(file_path) as in_file:
            for line in in_file:
                bits = line.rstrip().split('\t')
                dat[bits[0]] += float(bits[1])

    elif file_path.endswith('.txt'):
        dat = coll.Counter()
        with open(file_path) as in_file:
            for line in in_file:
                dat[line.rstrip().split('\t')[0]] = 1

    else:
        print "Skipped un-recognised file type:", file_path
        dat = ''

    return dat


def open_json(file_name):
    """
    :param file_name: the file name of the json archive to be opened
    :return: the nested peptide dictionary it contained
    """
    with open(file_name) as json_file:
        json_dat = json.load(json_file)

    return json_dat


def save_json(save_path, object_to_save):
    """
    :param save_path: Path to save to
    :param object_to_save: Thing to save as
    :return:
    """
    with open(save_path, 'w') as output:
        json.dump(object_to_save, output)

    return


def get_jy_ome_data():
    """
    :return: Two Counters, containing the transcriptome and proteome data from get-JY-ome-data.py
    """
    transcriptome = coll.Counter()
    proteome = coll.Counter()

    path_to_ome_data = '../Data/JY-omes/JY-omes.csv'
    with open(path_to_ome_data, 'rU') as in_file:
        line_count = 0
        for line in in_file:
            if line_count != 0:
                bits = line.rstrip().split(',')

                if bits[1]:
                    transcriptome[bits[0]] += float(bits[1])
                if bits[2]:
                    proteome[bits[0]] += float(bits[2])

            line_count += 1

    return transcriptome, proteome


def get_jy_ome_correlations(transcriptome, proteome, proteins):
    """
    :param transcriptome:
    :param proteome:
    :param proteins:
    :return: two dfs containing the value to plot correlation plots
    """
    transcriptome_corr = []
    proteome_corr = []

    for p in proteins.index:

        if p in transcriptome:
            donors = list(proteins.loc[p].index)
            vals = list(proteins.loc[p])
            for d in range(len(donors)):
                if transcriptome[p] == 0:
                    transcriptome_corr.append([donors[d], vals[d], 0])
                else:
                    transcriptome_corr.append([donors[d], vals[d], transcriptome[p]])

        if p in proteome:
            donors = list(proteins.loc[p].index)
            vals = list(proteins.loc[p])
            for d in range(len(donors)):
                if proteome[p] == 0:
                    proteome_corr.append([donors[d], vals[d], 0])
                else:
                    proteome_corr.append([donors[d], vals[d], proteome[p]])

    transcriptome_corr = list_to_df(transcriptome_corr,
                                    ['Sample', 'Peptidome Abundance', 'Transcriptome Abundance'], False)
    proteome_corr = list_to_df(proteome_corr, ['Sample', 'Peptidome Abundance', 'Proteome Abundance'], False)

    return transcriptome_corr, proteome_corr


def read_fa(ff):
    """
    :param ff: opened fasta file
    read_fa(file):Heng Li's Python implementation of his readfq function (tweaked to only bother with fasta)
    https://github.com/lh3/readfq/blob/master/readfq.py
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in ff:  # search for the start of the next record
                if l[0] in '>':  # fasta header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        # name, seqs, last = last[1:].partition(" ")[0], [], None # This version takes everything up to first space
        name, seqs, last = last[1:], [], None  # This version takes the whole line (post '>')
        for l in ff:  # read the sequence
            if l[0] in '>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:
            print "Input file does not appear to be a FASTA file - please check and try again"
            sys.exit()


def tidy_hla_name(short_str):
    """
    :param short_str: Shorter HLA gene name, as used by MHCflurry (e.g. HLA-A0101)
    :return: full proper gene gene with a linebreak for plotting, e.g. HLA-A*\n01:01
    """
    if len(short_str) == 9 and short_str.startswith('HLA-') and short_str[4] in ['A', 'B', 'C']:
        return short_str[:5] + '*\n' + short_str[5:7] + ':' + short_str[7:]
    else:
        raise IOError("Cannot convert HLA name, not in the format HLA-[ABC]xxxx")
