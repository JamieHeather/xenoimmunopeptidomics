# -*- coding: utf-8 -*-

"""
generate-data-files.py

Read in the various original source peptidome spreadsheets and output roughly comparably-formatted tab-separated-values.

Aiming to do the minimum amount of data processing as possible.

Where available, express abundances in non-transformed values (as 'count-like' as possible).

In order to avoid the issue of ambiguous PTMs (e.g. XXX[ss]XXXX), or difference in representation, I'm currently just
  removing any non-amino acid characters and treating all sequences purely as amino acid sequences currently. We will
  just treat the phosphopeptidome-derived sequences as phosphorylated *somewhere*.

I'm also ignoring any sequences <8 or >14 amino acids long, just to weed out the wackier sequences.

"""

from __future__ import division
import os
import sys
import mygene
import xlrd
import collections as coll
import functions as fxn

__version__ = '0.7.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def get_outfile(out_dir, sharing, type, sample_id, abundance):
    """
    Quick function to give the appropriate location to save the processed file
    :param out_dir: base data directory
    :param sharing: whether the data is to be shared (confidential/public)
    :param type: peptide or phospho
    :param sample_id: name to save under
    :abundance: y/n, determines file type to save as (tsv or txt)
    :return: path to appropriate file
    """

    if abundance is True:
        out_ext = '.tsv'
    elif abundance is False:
        out_ext = '.txt'
    else:
        print "Abundance error!"
        sys.exit()

    if sharing == 'confidential':
        outfile = out_dir + type + '/' + sample_id + '-confidential' + out_ext
    else:
        outfile = out_dir + type + '/' + sample_id + out_ext

    return outfile


def remove_ptms(peptide):
    """
    Remove PTMs from data from John Castle's and Donald Hunt's labs.
    Contains mass modifiers in format "AAAA(+XX.XX)AAA" - just need to remove between brackets. Can be multiple present.
    The Hunt data may also contain more brackets to indicate ambiguous PTMs, and has lower case modified residues.
      All other brackets (not accompanying a '+') can be deleted, and all lower case characters can be capitalised.
    :param peptide: String from spreadsheet's peptide column
    :return: peptide sequence with PTM indicators removed
    """

    while "(+" in peptide:
        open = peptide.find('(+')
        close = peptide.find(')')
        peptide = peptide[:open] + peptide[close+1:]

    peptide = peptide.replace('(', '').replace(')', '').upper()
    return peptide


def fix_mice_seqs(peptide):
    """
    Fix sequences from the mice data, which have similar PTMs as the Hunt/Castle data (with square brackets),
      but also have extra padding residues bracketed with periods (e.g. 'X.XXX[+YY.YY]XXX.X')
    :param peptide: From spreadsheet
    :return: Just PTM-free, bound, amino acid sequence
    """
    while "[+" in peptide:
        open = peptide.find('[+')
        close = peptide.find(']')
        peptide = peptide[:open] + peptide[close+1:]
    if '.' in peptide:
        peptide = peptide[peptide.find('.')+1:]
        peptide = peptide[:peptide.find('.')]
    return peptide.replace('(', '').replace(')', '').upper()


def get_mouse_id(sheet_name):
    nam = str(sheet_name)
    # nam = nam.replace('M', 'Mouse').replace(' ', '-')
    return nam


if __name__ == '__main__':

    fxn.check_scripts_dir()
    data_dir = '../Data/'
    published_dir = 'Raw-Files-Published/'
    agenus_dir = 'Raw-Files-Agenus/'
    other_dir = 'Raw-Files-Other-Lines/'

    # Generate output folders to write into if not already there
    dirs_to_make = ['phospho', 'peptide', 'phospho-prot', 'peptide-prot']
    for out_dir in dirs_to_make:
        if not os.path.exists(data_dir + out_dir):
            os.mkdir(data_dir + out_dir)

    # Define upper/lower limits for allowed peptides
    lower = 8
    upper = 15

    # Loop through all the files in the appropriate directory
    # Using the hyphen delimited metadata in the filenames, determine how to handle the file
    print "Processing published data..."
    all_input_files = [x for x in os.listdir(data_dir + published_dir) if x.startswith('.') is False]

    # Define xlsx properties (i.e. appropriate worksheet and peptide column indexes)
    worksheets = {'Bourdetsky2014': 0, 'Hassan2013': 0, 'BassaniSternberg2015': 1,
                  'Hassan2013HCC': 0, 'Caron2015': 0}
    peptide_cols = {'Bourdetsky2014': 0, 'Hassan2013': 0, 'BassaniSternberg2015': 4,
                    'Hassan2013HCC': 0, 'Caron2015': 0}

    for f in all_input_files:

        # Check it's not a directory
        if os.path.isdir(data_dir + published_dir + f):
            continue

        print '\t' + f

        bits = f.split('-')
        if len(bits) != 7:
            print "File name error - please check file:", f
        else:
            ref, source, data_type, sharing, abundance, sample_id, ext = bits

        dat = []
        abundance = False

        outfile = get_outfile(data_dir, sharing, data_type, sample_id, abundance)

        # SysMHC csv data
        if ext == '.csv' and source == 'SysMHC':

            with open(data_dir + published_dir + f, 'rU') as in_file:
                line_count = 0
                for line in in_file:
                    bits = line.rstrip().split(',')
                    if line_count == 0:
                        headers = bits
                    else:
                        peptide = bits[3]
                        if lower <= len(peptide) <= upper:
                            dat.append(peptide)
                    line_count += 1

        # xlsx data
        elif ext == '.xlsx':

            # This is the most varied source of data, with differing kinds of info, in different places in xlsx files.
            input_xlsx = xlrd.open_workbook(data_dir + published_dir + f)
            sheet = input_xlsx.sheet_by_index(worksheets[ref])
            peptide_col = [str(sheet.cell_value(x, peptide_cols[ref])) for x in range(sheet.nrows)][1:]

            if ref == 'Caron2015':
                cell_types = [str(sheet.cell_value(x, 7)) for x in range(sheet.nrows)][1:]
                if 'JY' in f:
                    peptide_col = [peptide_col[x] for x in range(len(cell_types)) if 'JY' in cell_types[x]]
                elif 'Jurkat' in f:
                    peptide_col = [peptide_col[x] for x in range(len(cell_types)) if 'Jurkat' in cell_types[x]]

            # Length/PTM filter, take unique
            peptide_col = [remove_ptms(x) for x in peptide_col if lower <= len(x) <= upper]
            dat = list(set(peptide_col))

        # And output...
        with open(outfile, 'w') as out_file:
            for pep in dat:
                out_file.write(pep + '\n')

    # Then on to process the mouse data, which is in very different formats, and will be processed more manually
    print "Processing Agenus data..."
    mg = mygene.MyGeneInfo()
    agenus_raw_files = [x for x in os.listdir(data_dir + agenus_dir)]

    abundance_key = {'phospho': True, 'peptide': False}
    prot_info = {'phospho-culture': 3, 'phospho-mice': 3, 'peptide-culture': 11, 'peptide-mice': 4}
    swissprot_key = {'phospho': False, 'peptide': True}
    swissprot_index = {'culture': 1, 'mice': 3}

    phosphos = coll.defaultdict(list)
    for f in agenus_raw_files:

        if f.startswith('~') or not f.endswith('xlsx'):
            continue

        print '\t' + f

        bits = f.split('-')
        data_type = bits[2]
        growth_type = bits[3]
        conf_key = {'peptide': 'public', 'phospho': 'confidential'}

        ref, source, sharing, abundance = ['Ag-' + growth_type + '-' + data_type, 'Agenus',
                                           conf_key[data_type], abundance_key[data_type]]

        input_xlsx = xlrd.open_workbook(data_dir + agenus_dir + f)
        input_sheets = input_xlsx.sheet_names()

        for s in range(len(input_sheets)):
            # Loop through each sheet of the xslx and extract peptides
            sheet_name = str(input_sheets[s])
            sample_name = get_mouse_id(sheet_name)

            print '\t\t' + sample_name

            pep_outfile = get_outfile(data_dir, sharing, data_type, sample_name, abundance)
            prot_outfile = get_outfile(data_dir, sharing, data_type + '-prot', sample_name, abundance)
            prot_index = prot_info[data_type + '-' + growth_type]
            swissprot = swissprot_key[data_type]

            # Pull out the peptide and abundance columns from that worksheet, standardising the peptide sequence formats
            sheet = input_xlsx.sheet_by_index(s)
            peptide_col = [fix_mice_seqs(str(sheet.cell_value(x, 0))) for x in range(sheet.nrows)][1:]

            phosphos_seqs = [str(sheet.cell_value(x, 0)) for x in range(sheet.nrows)
                             if '79.' in str(sheet.cell_value(x, 0))][1:]

            for ps in phosphos_seqs:
                phosphos[ps].append(sample_name)

            # Need to treat phospho and peptides differently, as only phosphos have abundance
            if abundance:
                abundance_col = [str(sheet.cell_value(x, 1)) for x in range(sheet.nrows)][1:]
                abundance_col = [float(x) for x in abundance_col]

                dat_list = coll.defaultdict(list)
                for i in range(len(peptide_col)):
                    pep = fix_mice_seqs(peptide_col[i])
                    if lower <= len(pep) <= upper:
                        # Phospho data has a specific quirk - may have multiple abundance entries if co-eluting phosphos
                        if (data_type == 'peptide') or (abundance_col[i] not in dat_list[pep]):
                            dat_list[pep].append(abundance_col[i])

                # Sum the values of each entry for a given peptide
                pep_dat = coll.Counter()
                for pep in dat_list:
                    pep_dat[pep] = sum(dat_list[pep])

                # Finally write out to new, standardised files
                with open(pep_outfile, 'w') as out_file:
                    for pep in pep_dat.most_common():
                        out_file.write(pep[0] + '\t' + str(pep[1]) + '\n')

            else:
                with open(pep_outfile, 'w') as out_file:
                    for pep in peptide_col:
                        if lower <= len(pep) <= upper:
                            out_file.write(pep + '\n')

            # Output number of occurrences of peptides from specific proteins
            prot_col = [str(sheet.cell_value(x, prot_index)) for x in range(sheet.nrows)][1:]

            if swissprot:  # Add pipes to protein header to allow splitting cases with no protein matches
                prot_col = [(x + '|||').split('|')[swissprot_index[growth_type]] for x in prot_col]
                # Using mygene to convert Uniprot IDs to symbols to ensure standardised versioning

            # Then convert to protein sequences, allowing for server timeouts
            symbol_query = mg.querymany(prot_col, scopes='symbol,reporter,accession,entrezgene,uniprot',
                                        fields='symbol', species='human', returnall=True)
            symbols = coll.defaultdict(list)
            for query in symbol_query['out']:
                if 'symbol' in query:
                    if query['symbol'] not in symbols[query['query']]:
                        symbols[query['query']].append(query['symbol'])

            # Throw into a dict and output
            prot_dat = coll.Counter()
            for i in range(len(prot_col)):
                if prot_col[i]:
                    prot_list = symbols[prot_col[i]]
                    if abundance:
                        amount = abundance_col[i]
                    else:
                        amount = 1
                    for prot in prot_list:  # Allow for rare multi-mapping peptides
                        prot_dat[str(prot)] += amount

            with open(prot_outfile, 'w') as out_file:
                for prot in prot_dat.most_common():
                    out_file.write(prot[0] + '\t' + str(prot[1]) + '\n')

    # Finally process the other cell line phospho data...
    print "Processing additional peptide files..."
    pp_files = [x for x in os.listdir(data_dir + other_dir) if 'phospho' in x and x.endswith('.xlsx')]

    for fl in pp_files:
        nam = fl.split('-')[0]
        if nam == 'MDAMB436':
            nam = 'MDA-MB-436'

        input_xlsx = xlrd.open_workbook(data_dir + other_dir + fl)
        sheet_index = {'culture': 1, 'mice': 0}
        for wsheet in sheet_index:
            sheet = input_xlsx.sheet_by_index(sheet_index[wsheet])
            peptide_col = [remove_ptms(str(sheet.cell_value(x, 0))) for x in range(sheet.nrows)][1:]
            out_nam = wsheet[0].upper() + '-' + nam
            pep_outfile = get_outfile(data_dir, 'confidential', 'phospho', out_nam, False)
            with open(pep_outfile, 'w') as out_file:
                for pep in peptide_col:
                    out_file.write(pep + '\n')

    # ... and their non-phospho data
    print "Processing additional peptide files..."
    peptide_files = [x for x in os.listdir(data_dir + other_dir) if 'peptide-public' in x and x.endswith('.xlsx')]

    for fl in peptide_files:
        nam = fl.split('-')[0]
        if nam == 'MDAMB436':
            nam = 'MDA-MB-436'

        input_xlsx = xlrd.open_workbook(data_dir + other_dir + fl)
        sheet = input_xlsx.sheet_by_index(0)
        peptide_col = [remove_ptms(str(sheet.cell_value(x, 1))) for x in range(sheet.nrows)][1:]
        culture_col = [sheet.cell_value(x, 5) for x in range(sheet.nrows)][1:]
        mouse_col = [sheet.cell_value(x, 6) for x in range(sheet.nrows)][1:]
        c_outfile = get_outfile(data_dir, 'public', 'peptide', 'C-'+ nam, False)
        m_outfile = get_outfile(data_dir, 'public', 'peptide', 'M-' + nam, False)

        with open(c_outfile, 'w') as c_out_file, open(m_outfile, 'w') as m_out_file:
            for i in range(len(peptide_col)):
                if lower <= len(peptide_col[i]) <= upper:
                    if culture_col[i] in ['Yes', 'MS1']:
                        c_out_file.write(fix_mice_seqs(peptide_col[i]) + '\n')
                    if mouse_col[i] in ['Yes', 'MS1']:
                        m_out_file.write(fix_mice_seqs(peptide_col[i]) + '\n')
