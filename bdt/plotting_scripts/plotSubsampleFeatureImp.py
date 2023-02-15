import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()
year = args.year

import os, sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas, root_numpy
import argparse
import pickle
import sklearn
from   sklearn.externals import joblib
from   sklearn.ensemble  import GradientBoostingClassifier
from   sklearn.model_selection import train_test_split
from   pdb import set_trace
from   array       import array

from   xgboost import XGBClassifier, plot_importance
import pdb

def get_xgb_imp(xgb, feat_names):
    from numpy import array
    imp_vals = xgb.booster().get_fscore()
    imp_dict = {feat_names[i]:float(imp_vals.get('f'+str(i),0.)) for i in range(len(feat_names))}
    total = array(imp_dict.values()).sum()
    return {k:v/total for k,v in imp_dict.items()}

plt.clf()
merge_dict = {}
for isample in range(11):

    tags = ['_noIP2D'] 
    bdttag = tags[0]

    clfs = []
    dicts = {}
    for ii, itag in enumerate(tags):
        clfs.append(joblib.load(open('../results/classifier_%s_%s_%s.pkl'%(itag,year,isample)))  )

        feat_imp = clfs[ii].get_booster().get_fscore().items()
        vals = np.array(feat_imp)
        feat_dict = dict(feat_imp)
        vals = [v for k,v in feat_dict.items()]
        max_val = float(max(vals))
        
        for k,v in feat_dict.items():
            feat_dict[k] = feat_dict[k]/max_val
        
        for k,v in feat_dict.items():
            try:
                merge_dict[k].append(v)
            except:
                merge_dict[k] = [v]

## now compute average over samples
for k,v in merge_dict.items():
    merge_dict[k] = sum(merge_dict[k]) / len(merge_dict[k])

## sort
from collections import OrderedDict
merge_dict_s = OrderedDict(sorted(merge_dict.items(), key=lambda x: x[1]))

feats = list(merge_dict_s.keys())
y_pos = np.arange(len(feats))
score = list(merge_dict_s.values())

plt.rcdefaults()
fig, ax = plt.subplots()

ax.barh(y_pos, score, align='center', color='MediumSeaGreen', lw=0)
ax.set_yticks(y_pos)
ax.set_yticklabels(feats)
# ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Average norm. feature importance')
ax.set_title(year)

for i, v in enumerate(score):
    ax.text(v-0.1, i, s='%.2f'%v, ha='center', va='center', color='white',fontweight='bold')

ax.set_xlim(0,1.1)
ax.set_ylim(-1., len(y_pos))
plt.show()

plt.tight_layout()
plt.savefig('../results/feature_imp_allsamples%s_%s_noIP2D.pdf'%(bdttag, year) )
