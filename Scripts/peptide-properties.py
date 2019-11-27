# -*- coding: utf-8 -*-

"""
peptide-properties.py

Analyse the electrochemical properties of the peptides not predicted to bind to any of the relevant MHC alleles
See https://modlamp.org/modlamp.html#modlamp.descriptors.GlobalDescriptor.charge_density

"""

from __future__ import division
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import functions as fxn
import os
import collections as coll
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def mwu(list1, list2):
    """
    :param list1: One list of numbers
    :param list2: A second list of numbers
    :return: P value for a Mann-Whitney U test (unpaired, non-parametric)
    """
    return stats.mannwhitneyu(list(list1), list(list2))[1]


def get_peptide_values(list_peptides, descriptor_name):
    """
    :param list_peptides: List of amino acid peptides
    :param descriptor_name: MODLamp-prescribed descriptor name
    :return: corresponding values for that descriptor for each of the peptides in the input list
    """
    properties = PeptideDescriptor(list_peptides, descriptor_name)
    properties.calculate_moment()
    return [x[0] for x in properties.descriptor]


data_dir = '../Data/NonPredictedBinders/'

plot_dir = fxn.plot_dir + fxn.get_date() + '-NPB-peptide-properties/'
if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)

pep_files = [x for x in os.listdir(data_dir) if x.endswith('.txt')]

peptides = coll.defaultdict(list)

for f in pep_files:
    nam = f.split('-')[0]
    with open(data_dir + f, 'rU') as in_file:
        for line in in_file:
            peptides[nam].append(line.rstrip())


eisenbergs = coll.defaultdict(list)
eisenbergs_long = []
charges = coll.defaultdict(list)
charges_long = []
charge_densities = coll.defaultdict(list)
charge_densities_long = []
polarities = coll.defaultdict(list)
polarities_long = []
gravy = coll.defaultdict(list)
gravy_long = []

for gp in peptides:
    #
    eisenbergs[gp] = get_peptide_values(peptides[gp], 'eisenberg')
    for val in eisenbergs[gp]:
        eisenbergs_long.append([gp, val])
    #
    properties = GlobalDescriptor(peptides[gp])
    properties.calculate_charge(ph=7.4, amide=True)
    charges[gp] = [x[0] for x in properties.descriptor]
    for val in charges[gp]:
        charges_long.append([gp, val])
    #
    properties = GlobalDescriptor(peptides[gp])
    properties.charge_density(ph=7.4, amide=True)
    charge_densities[gp] = [x[0] for x in properties.descriptor]
    for val in charge_densities[gp]:
        charge_densities_long.append([gp, val])
    #
    polarities[gp] = get_peptide_values(peptides[gp], 'polarity')
    for val in polarities[gp]:
        polarities_long.append([gp, val])
    #
    gravy[gp] = get_peptide_values(peptides[gp], 'gravy')
    for val in gravy[gp]:
        gravy_long.append([gp, val])


eisenbergs_long = fxn.list_to_df(eisenbergs_long, ['Group', 'Hydrophobicity'], False)
charges_long = fxn.list_to_df(charges_long, ['Group', 'Charge'], False)
charge_densities_long = fxn.list_to_df(charge_densities_long, ['Group', 'Charge density'], False)
polarities_long = fxn.list_to_df(polarities_long, ['Group', 'Polarity'], False)
gravy_long = fxn.list_to_df(gravy_long, ['Group', 'Hydrophobicity'], False)

ordr = ['Cultured', 'Mouse', 'Both']

sns.violinplot(data=eisenbergs_long, x='Group', y='Hydrophobicity', order=ordr)
plt.xlabel("")
plt.ylabel("Hydrophobicity (Eisenberg scale)")
plt.savefig(plot_dir + 'eisenberg-hydrophobicity.png', dpi=300, bbox_inches='tight')
plt.close()

sns.violinplot(data=gravy_long, x='Group', y='Hydrophobicity', order=ordr)
plt.xlabel("")
plt.ylabel("Hydrophobicity (GRAVY scale)")
plt.savefig(plot_dir + 'gravy-hydrophobicity.png', dpi=300, bbox_inches='tight')
plt.close()

sns.violinplot(data=charges_long, x='Group', y='Charge', order=ordr)
plt.xlabel("")
plt.savefig(plot_dir + 'charge.png', dpi=300, bbox_inches='tight')
plt.close()

sns.violinplot(data=charge_densities_long, x='Group', y='Charge density', order=ordr)
plt.xlabel("")
plt.savefig(plot_dir + 'charge-densities.png', dpi=300, bbox_inches='tight')
plt.close()

sns.violinplot(data=polarities_long, x='Group', y='Polarity', order=ordr)
plt.xlabel("")
plt.savefig(plot_dir + 'polarity.png', dpi=300, bbox_inches='tight')
plt.close()


# Stats
# >>> mwu(eisenbergs['Cultured'], eisenbergs['Mouse'])
# 0.2880000876829417
# >>> mwu(eisenbergs['Cultured'], eisenbergs['Both'])
# 0.49520230961139955
# >>> mwu(eisenbergs['Mouse'], eisenbergs['Both'])
# 0.30555391628097484

# >>> mwu(charges['Cultured'], charges['Mouse'])
# 5.28838236429353e-07
# >>> mwu(charges['Cultured'], charges['Both'])
# 0.00427644151209229
# >>> mwu(charges['Mouse'], charges['Both'])
# 0.11709912122209087

# >>> mwu(charge_densities['Cultured'], charge_densities['Mouse'])
# 6.568512587149028e-06
# >>> mwu(charge_densities['Cultured'], charge_densities['Both'])
# 0.0281842424179716
# >>> mwu(charge_densities['Mouse'], charge_densities['Both'])
# 0.05657544597005422

# >>> mwu(gravy['Cultured'], gravy['Mouse'])
# 0.12405181679526123
# >>> mwu(gravy['Cultured'], gravy['Both'])
# 0.4580397469761439
# >>> mwu(gravy['Mouse'], gravy['Both'])
# 0.1737808782744939

# >>> mwu(polarities['Cultured'], polarities['Mouse'])
# 0.026223491556441856
# >>> mwu(polarities['Cultured'], polarities['Both'])
# 0.21350762590659111
# >>> mwu(polarities['Mouse'], polarities['Both'])
# 0.24293670370628512
