# -*- coding: utf-8 -*-

"""
get-JY-ome-data.py

Compile transcriptome and proteome data from published sources into a single csv.

Note that the TPM data comes from Kallisto mapping of publicly available RNAseq data.
The proteome data comes from the 2015 Bassani-Sternberg paper.

"""

from __future__ import division
import functions as fxn
import collections as coll
import gzip
import xlrd
import sys
import numpy as np
import urllib2

__version__ = '0.3.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def get_row_data(data, desired_row):
    """
    :param data: A worksheet of the data file, parsed by xlrd
    :param desired_row: The number of the row in question
    :return: a list of the fields of that row
    """
    return [data.cell_value(desired_row, x) for x in range(data.ncols)]


if __name__ == '__main__':

    fxn.check_scripts_dir()
    dat_path = '../Data/JY-omes/'

    tpm_file = 'jy-tpm.csv'
    map_file = 'Ensembl-mappings.tsv.gz'
    bs_file = 'BS-MCP2015-JY-data.xlsx'

    # First get the mappings of ensembl gene IDs to gene symbols
    print "Reading in Ensembl symbol mappings..."
    mapping = {}
    with gzip.open(dat_path + map_file, 'rU') as in_file:
        line_count = 0
        for line in in_file:
            bits = line.replace('\"', '').rstrip().split('\t')
            if line_count == 0:
                headers = bits
            else:
                ensg = bits[2]
                gs = bits[3]
                if ensg not in mapping.keys():
                    mapping[ensg] = gs
                else:
                    if mapping[ensg] == gs:
                        continue
                    else:
                        print 'Error'
                        sys.exit()
            line_count += 1

    # Let's read in the proteome data
    print "Reading in JY proteome intensity values..."
    xlsx_file = xlrd.open_workbook(dat_path + bs_file)
    sheet = xlsx_file.sheet_by_index(1)

    prot_dat = coll.defaultdict(list)
    for i, row in enumerate(range(sheet.nrows)):
        row_data = get_row_data(sheet, row)
        if i == 0:  # Get header information
            headers = row_data
        else:
            ensid = row_data[0]
            if ensid in mapping.keys():
                prot_dat[mapping[row_data[0]]].append(row_data[1])

    avg_prot_dat = {}
    for p in prot_dat:
        if len(prot_dat[p]) == 1:
            avg_prot_dat[p] = prot_dat[p][0]
        else:
            avg_prot_dat[p] = np.log2(np.mean([2**x for x in prot_dat[p]]))

    print "Reading in JY transcriptome values..."
    # Then finally the transcript data
    tpm_dat = {}
    with open(dat_path + tpm_file, 'rU') as in_file:
        line_count = 0
        for line in in_file:
            bits = line.rstrip().split(',')
            if line_count == 0:
                headers = bits
            else:
                val = float(bits[1])
                if val != 0:
                    tpm_dat[bits[0]] = np.log2(val)
            line_count += 1

    all_genes = list(set(avg_prot_dat.keys()) | set(tpm_dat.keys()))
    all_genes.sort()

    out_file_path = '../Data/JY-omes/JY-omes.csv'
    print "Writing JY-ome values out to", out_file_path + "..."
    with open(out_file_path, 'w') as out_file:
        out_file.write('Gene,Transcripts-log2TPM,Protein-log2int\n')
        for g in all_genes:
            if g in tpm_dat and g in avg_prot_dat:
                out_file.write(g + ',' + str(tpm_dat[g]) + ',' + str(avg_prot_dat[g]) + '\n')
            elif g in avg_prot_dat and g not in tpm_dat:
                out_file.write(g + ',,' + str(avg_prot_dat[g]) + '\n')
            elif g in tpm_dat and g not in avg_prot_dat:
                out_file.write(g + ',' + str(tpm_dat[g]) + ',\n')

    print "\t... done!"

    print "Downloading proteome files..."

    # os.chdir('../Data/')
    for proteome in [x for x in fxn.proteome_key if 'UP' in x]:
        print "\tGetting", fxn.proteome_key[proteome], "proteome file..."
        out_path = '../Data/' + fxn.get_date() + '_' + proteome + '_' + fxn.proteome_key[proteome] + '.fasta.gz'
        url = urllib2.urlopen('https://www.uniprot.org/uniprot/?query=proteome:'
                              + proteome + '&format=fasta&compress=yes')
        url_data = url.read()

        with open(out_path, 'wb') as out_file:
            out_file.write(url_data)
