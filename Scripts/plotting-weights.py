# -*- coding: utf-8 -*-

"""
plotting-weights.py

Extra supplementary plotting, regarding the yields of the different growth types

"""

import functions as fxn
import seaborn as sns
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy import stats


__version__ = '0.2.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def get_regression_str(x, y):
    """
    :param x: list of X values
    :param y: list of y values
    :return: a string detailing the R squared and P values of the linear regression of X and Y
    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_2 = r_value * r_value
    return "R2: " + "{0:.3f}".format(r_2) + ' P: ' + format_p_vals(p_value)


def format_p_vals(number):
    """
    :param number: a float between 0 and 1 (a p value)
    :return: string of formatted p value (to 3 sig figs if > 0.001, otherwise in scientific notation)
    """
    if number > 0.001:
        return "{0:.3f}".format(number)
    else:
        return "{:.2E}".format(number)

# if __name__ == '__main__':

# Initialise directory
fxn.check_scripts_dir()
sns.set(font="Arial", font_scale=1.5)
plot_dir = fxn.plot_dir + fxn.get_date() + '-weights/'
if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)

# Read in data
weights = []
with open('../Data/sample-weights.csv', 'rU') as in_file:
    line_count = 0
    for line in in_file:
        bits = line.rstrip().split(',')
        if line_count == 0:
            headers = bits
        else:
            weights.append([int(x) for x in bits[:2]] + [bits[2]])
        line_count += 1

headers[0] = 'Day'
headers[2] = 'Growth type'
weights = fxn.list_to_df(weights, headers, False)
weights = weights.sort_values(by='Growth type')

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
sns.scatterplot(data=weights, x='Day', y='mg', hue='Growth type')

plt.ylabel("Cumulative sample weight (mg)")
plt.xlabel("Days of growth")
plt.savefig(plot_dir + 'weight-over-time.png', dpi=300, bbox_inches='tight')
plt.close()

# Read in number of peptides to check correlation with input xenograft tumor weight
with open('../Data/mouse-weights-v-peptides.csv', 'rU') as in_file:
    wvp = []
    line_count = 0
    for line in in_file:
        bits = line.rstrip().split(',')
        if line_count == 0:
            headers = bits
        else:
            bits[1:] = [int(x) for x in bits[1:]]
            wvp.append(bits)
        line_count += 1

wvp = fxn.list_to_df(wvp, headers, True)

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
g = sns.lmplot(data=wvp, x='Tumor weight (mg)', y='Number peptides',
               line_kws={'color': 'grey'},  scatter_kws={'color': fxn.orange})

for ax in g.axes.flat:
    ax.text(0.7, 1.02, get_regression_str(wvp['Tumor weight (mg)'], wvp['Number peptides']),
            color=fxn.orange, transform=ax.transAxes, fontsize=10)


plt.ylabel("Number of peptides")
plt.savefig(plot_dir + 'weight-over-peptides.png', dpi=300, bbox_inches='tight')
plt.close()
