# import argparse
# parser = argparse.ArgumentParser(description="")
# parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
# args = parser.parse_args()
# year = args.year
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

import os, sys
# sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import pandas, root_numpy
from copy import deepcopy
import argparse
import pickle
import sklearn
from   sklearn.externals import joblib
from   sklearn.ensemble  import GradientBoostingClassifier
from   sklearn.metrics   import roc_curve
from   sklearn.model_selection import train_test_split
from   pdb import set_trace
from   array       import array

from ROOT      import TFile, TTree, TH1F, gROOT, TChain

mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

gROOT.SetBatch(True)

##########################################################################################
#####   SIGNAL AND BACKGROUND SELECTION
##########################################################################################
sig_mass_str = '(tagged_mass > {M} - 2.5*{S}) & (tagged_mass < {M} + 2.5*{S})'.format( M=mc_mass,S=mc_sigma)    

truth_match_str = '( (truthMatchMum==1) & (truthMatchMup ==1) & (truthMatchTrkm==1) & (truthMatchTrkp==1) )'

ct_str = '( ( (tagB0==1) & (genSignal==1)) | ( (tagB0==0) & (genSignal==2) ) )'

bkg_mass_str = '( ( (tagged_mass > {M}-7.*{S} ) & (tagged_mass < {M}-3.*{S} ) ) | ( (tagged_mass > {M}+3.*{S} ) & (tagged_mass < {M}+7*{S} ) ) )' .format( M=mc_mass,S=mc_sigma) 

## consider mass preselection which is different in 2016 and 2018 -> it is modified below
mass_range_str = '(mumuMass < 2.702)'

years = ['2016', '2017', '2018' ]

for year in years:
  if year == '2016': 
    mass_range_str = '(mumuMass < 2.702 || mumuMass > 4.)'


  sig_selection_cutbased = sig_mass_str + ' && ' + \
                           truth_match_str + ' && ' +\
                           ct_str + ' && ' +\
                           'trig == 0 && ' +\
                           mass_range_str
  
  bkg_selection_cutbased = bkg_mass_str + ' && ' + \
                           mass_range_str
  

  sig_list = [
            '../sub_samples/sample_%s_MC_LMNR_0_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_1_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_2_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_3_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_4_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_5_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_6_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_7_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_8_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_9_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_MC_LMNR_10_add_vars_passPreselection_addBDT_noIP2D_%s.root'%(year,year),
  ]
  
  bkg_list = [
            '../sub_samples/sample_%s_data_LMNR_0_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_1_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_2_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_3_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_4_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_5_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_6_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_7_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_8_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_9_add_vars_addBDT_noIP2D_%s.root'%(year,year),
            '../sub_samples/sample_%s_data_LMNR_10_add_vars_addBDT_noIP2D_%s.root'%(year,year),
  ]

  branches = [
      'pass_preselection',
      'bdt_prob',
      'tagged_mass',
      'mumuMass',
  ]


  branches = list(set(branches))

  sig = pandas.DataFrame(
      root_numpy.root2array(
          sig_list, 
          'ntuple',
          branches  = branches + ['weight', 'truthMatchMum', 'truthMatchMup', 'truthMatchTrkm', 'truthMatchTrkp', 'genSignal', 'tagB0', 'trig'],
          selection = sig_selection_cutbased,
      )
  )

  bkg = pandas.DataFrame(
      root_numpy.root2array(
          bkg_list, 
          'ntuple',
          branches  = branches,
          selection = bkg_selection_cutbased,
      )
  )

  sig['target'] = np.ones (sig.shape[0]).astype(np.int)
  bkg['target'] = np.zeros(bkg.shape[0]).astype(np.int)
  sig.insert(2, 'label', 'signal', True)
  bkg.insert(2, 'label', 'background', True)
  data_all = pandas.concat([sig, bkg])
  data = data_all[data_all.pass_preselection==1]

#     plt.plot(data.bdt_prob, 'data.target=0')
#     data.plot('bdt_prob', ['target'])#, color='limegreen')
#     ax = data.hist(column='bdt_prob' bins=120, alpha=0.5, by='target')

  fig, ax = plt.subplots()
  for name, group in data.groupby('target'):
      plt.hist(group['bdt_prob'], bins=200, alpha=0.5, label=set(group.label), histtype='stepfilled', normed=True, log=True)
      plt.yscale('log', nonposy='clip')
  
  plt.legend(loc='best')
  # ax.set_ylim(0.03, 100)
  ax.set_title(year)
  plt.xlabel('BDT output')
  plt.ylabel('Arbitrary units')
  plt.savefig('../results/bdt_output_%s.pdf' %year )
