# -*- coding: utf-8 -*-

"""
run-analysis.py

Run each of the relevant scripts in turn, in the right order, to produce the data for analysis.

Note that due to intellectual property interests that phosphopeptide sequences have been omitted from the repo.
The generate-data-files.py script has thus been run in advance, and is kept in here for reference.

"""

import os

if not os.path.exists('../Plots/'):
    os.mkdir('../Plots/')

scripts_order = ['generate-data-files.py',
                 'get-omics-data.py',
                 'plotting-weights.py',
                 'peptide-analysis.py',
                 'published-comparison.py',
                 # 'netmhc-check.py',  # Requires netMHC installed/pre-run locally first
                 'mouse-proteome-check.py',
                 # 'phospho-analysis.py'  # Requires phosphopeptide data not included in repository
                 ]

for script in scripts_order:
    print script
    execfile(script)
