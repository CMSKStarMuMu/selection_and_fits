import argparse
parser = argparse.ArgumentParser(description="")
# parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--tag"      , help = "", default = 'punzi_removeTkMu_fixBkg_2017')
args = parser.parse_args()
# year = args.year

import os, sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pandas, root_numpy
from copy import deepcopy
import argparse
import pickle
import sklearn
from   sklearn.externals import joblib
from   sklearn.ensemble  import GradientBoostingClassifier
from   sklearn.metrics   import roc_curve, roc_auc_score
from   sklearn.model_selection import train_test_split
from   array import array

from   xgboost import XGBClassifier, plot_importance

from palettable.matplotlib import Viridis_15
from ROOT      import TFile, TTree, TH1F, gROOT, TChain
import pdb

mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

gROOT.SetBatch(True)

## weight input background mass
min_m  = 4.8
max_m  = 5.8
n_bins = 500

colors = [
    'darkslateblue',
    'plum',
    'lightcoral',
    'lightsalmon',
    'navajowhite',
    'gold',
    'khaki',
    'yellow',
    'aliceblue',          
    'antiquewhite',       
    'aqua',               
    'aquamarine',         
    'azure',              
    'beige',              
    'bisque',             
    'black',              
    'blanchedalmond',     
    'mediumslateblue',
    'mediumpurple',
    'mediumorchid',
]

# years = [ '2018' ]
years = ['2016', '2017', '2018' ]
nsamples = 11


tags = [
    '_noIP2D',
#     '_noIP2D_noTrkSign'
] 


tag_labels = [
    'no trkMinIP',
#     'no trkMinIP, no Trk Sign',
] 
kinds = ['test', 'train']


metrics = ['auc', 'fpr', 'tpr', 'thresholds']
samples = ['train', 'test']

results_tags = dict(zip(tags, [ 
                              dict(zip(samples, [
                                                 dict( zip(metrics, [ [] for ii in range(len(metrics))])) 
                              for jj in range(len(samples)) ] ) ) 
               for kk in range(len(tags))] ))



sig_mass_str = '(tagged_mass > {M} - 2.5*{S}) & (tagged_mass < {M} + 2.5*{S})'.format( M=mc_mass,S=mc_sigma)    

truth_match_str = '( (truthMatchMum==1) & (truthMatchMup ==1) & (truthMatchTrkm==1) & (truthMatchTrkp==1) )'

ct_str = '( ( (tagB0==1) & (genSignal==1)) | ( (tagB0==0) & (genSignal==2) ) )'

bkg_mass_str = '( ( (tagged_mass > {M}-7.*{S} ) & (tagged_mass < {M}-3.*{S} ) ) | ( (tagged_mass > {M}+3.*{S} ) & (tagged_mass < {M}+7*{S} ) ) )' .format( M=mc_mass,S=mc_sigma) 

## consider mass preselection which is different in 2016 and 2018 -> it is modified below
mass_range_str = '(mumuMass < 2.702)'


plt.clf()

import itertools
xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
plt.plot(xy, xy, color='grey', linestyle='--')
plt.xlim([10**-3, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')


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
            '../sub_samples/sample_%s_MC_LMNR_0_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_1_add_vars_passPreselection.root'%year,
            '../sub_samples/sample_%s_MC_LMNR_2_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_3_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_4_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_5_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_6_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_7_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_8_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_9_add_vars_passPreselection.root'%year, 
            '../sub_samples/sample_%s_MC_LMNR_10_add_vars_passPreselection.root'%year,
  ]
  
  bkg_list = [
            '../sub_samples/sample_%s_data_LMNR_0_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_1_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_2_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_3_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_4_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_5_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_6_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_7_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_8_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_9_add_vars.root'%year,
            '../sub_samples/sample_%s_data_LMNR_10_add_vars.root'%year,
  ]


  for isample in range(nsamples):

    branches = [
      'bCosAlphaBS',
      'bLBS/bLBSE',
      'trkmDCASign',
      'trkpDCASign',
      'bVtxCL',
      'bDCABS/bDCABSE',
      'kstarmass',
      'sum_isopt_04',
      'tagged_mass',
      'pass_preselection',
    ]

    branches = list(set(branches))
    
    sig = pandas.DataFrame(
        root_numpy.root2array(
            sig_list[:isample] + sig_list[(isample + 1):], 
            'ntuple',
            branches  = branches + ['weight'],
            selection = sig_selection_cutbased,
        )
    )
    
    bkg = pandas.DataFrame(
        root_numpy.root2array(
            bkg_list[:isample] + bkg_list[(isample + 1):], 
            'ntuple',
            branches  = branches,
            selection = bkg_selection_cutbased,
        )
    )

    ##########################################################################################
    #####   DEFINE THE TARGETS
    ##########################################################################################
    
    sig['target'] = np.ones (sig.shape[0]).astype(np.int)
    bkg['target'] = np.zeros(bkg.shape[0]).astype(np.int)
    
    sig['normfactor'] = sig.weight#*sig.m_weight #/len(sig) *10000
    bkg['normfactor'] = 1. ##bkg.m_weight# /len(bkg) #* 10000
    
    data_all = pandas.concat([sig, bkg])
    
    data = data_all[data_all.pass_preselection==1]
    train, test = train_test_split(data, test_size=0.30, random_state = 17)

    ##########################################################################################
    # #####   LOAD and PREDICT ON THE TEST UNBIAS SAMPLE
    # ##########################################################################################
    clfs = []
    pred = []
    pred_train = []
    
    for ii, itag in enumerate(tags):

        if '_noIP2D_noTrkSign' in itag:
          features = [
            'bCosAlphaBS',
            'bLBS/bLBSE',
            'bVtxCL',
            'bDCABS/bDCABSE',
            'kstarmass',
            'sum_isopt_04',
            ]
        elif '_noIP2D' in itag:
          features = [
            'bCosAlphaBS',
            'bLBS/bLBSE',
            'trkmDCASign',
            'trkpDCASign',
            'bVtxCL',
            'bDCABS/bDCABSE',
            'kstarmass',
            'sum_isopt_04',
            ]

        clfs.append(joblib.load(open('results/classifier_%s_%s_%s.pkl'%(itag,year,isample)))  )
        best_iter = clfs[-1].best_iteration
        pred.append(       clfs[ii].predict_proba(test[features], ntree_limit=best_iter)[:, 1])
        pred_train.append( clfs[ii].predict_proba(train[features], ntree_limit=best_iter)[:, 1])


#     pdb.set_trace()
    ##########################################################################################
    #####   ROC CURVE
    ##########################################################################################
    
    # ## draw ROC 
    fprs = []
    tprs = []
    thresholds = []
    fprs_train = []
    tprs_train = []
    thresholds_train = []
    
    for ii, itag in enumerate(tags):
        print 'saving train/test results for tag=%s for sample %s'%(itag,isample) + '  auc = %s'%roc_auc_score(test.target, pred[ii])
        fpr, tpr , threshold = roc_curve(test.target, (pred[ii] * test.pass_preselection) ) ### !! changed from test unbias
        fprs.append(fpr)
        tprs.append(tpr)
        thresholds.append(threshold)

        results_tags[itag]['test']['auc'].append(roc_auc_score(test.target, pred[ii]))
        results_tags[itag]['test']['fpr'].append(fpr)
        results_tags[itag]['test']['tpr'].append(tpr)
        results_tags[itag]['test']['thresholds'].append(threshold)


        fpr_tr, tpr_tr, threshold_tr = roc_curve(train.target, pred_train[ii])
        fprs_train.append(fpr_tr)
        tprs_train.append(tpr_tr)
        thresholds_train.append(threshold_tr)

        results_tags[itag]['train']['auc'].append(roc_auc_score(train.target, pred_train[ii]))
        results_tags[itag]['train']['fpr'].append(fpr_tr)
        results_tags[itag]['train']['tpr'].append(tpr_tr)
        results_tags[itag]['train']['thresholds'].append(threshold_tr)

#         results_tags[itag]['train']['auc'].append(roc_auc_score(train.target, pred_train[ii]))
#         results_tags[itag]['train']['fpr'].append(fpr_tr)
#         results_tags[itag]['train']['tpr'].append(tpr_tr)
#         results_tags[itag]['train']['thresholds'].append(threshold_tr)

#         plt.plot(fpr_tr, tpr_tr, color=colors[isample], linestyle='dotted', label = 'train sample %s'%isample )
#         plt.plot(fpr_tr, tpr_tr, color=colors[isample+ii], linestyle='dotted', label = 'train sample, %s'%itag.replace('_', '') )


# pdb.set_trace()

plt.clf()

c_fill      = 'rgba(52, 152, 219, 0.2)'
c_line      = 'rgba(52, 152, 219, 0.5)'
c_line_main = 'rgba(41, 128, 185, 1.0)'
c_grid      = 'rgba(189, 195, 199, 0.5)'
c_annot     = 'rgba(149, 165, 166, 0.5)'
c_highlight = 'rgba(192, 57, 43, 1.0)'
fpr_mean    = np.linspace(0, 1, 1000)


for ik,kind in enumerate(kinds):
    print 'ikind: ', kind
    for ij, itag in enumerate(tags):
      interp_tprs = []
      results = results_tags[itag]
      
    
      print 'itag: ', itag
      for i in range(nsamples*len(years)):
          print 'isample: ', i
          fpr           = results[kind]['fpr'][i]
          tpr           = results[kind]['tpr'][i]
          interp_tpr    = np.interp(fpr_mean, fpr, tpr)
          interp_tpr[0] = 0.0
          interp_tprs.append(interp_tpr)
    
#         plt.plot(
#           fpr,
#           tpr,
#           color=colors[ij],
#           label=r"ROC subsample %s (AUC = %0.4f $\pm$ %0.4f) %s" % (i,results[kind]['auc'][i]),
#           lw=2,
#           alpha=0.2,
#         )

      tpr_mean     = np.mean(interp_tprs, axis=0)
      tpr_mean[-1] = 1.0
      auc          = np.mean(results[kind]['auc'])
      std_auc = np.std(results[kind]['auc'])
    
#       pdb.set_trace()
      print 'plot'
      plt.plot(
        fpr_mean,
        tpr_mean,
        color=colors[ik*3+ij],
        label=r"AveROC (AUC = %0.4f $\pm$ %0.4f) %s" % (auc, std_auc, kind),
#         label=r"AveROC (AUC = %0.4f $\pm$ %0.4f) %s %s" % (auc, std_auc, tag_labels[ij], kind),
        lw=2,
        alpha=0.8,
      )
    
      tpr_std      = np.std(interp_tprs, axis=0)
    
      # std_tpr = np.std(interp_tprs, axis=0)
      tprs_upper = np.minimum(tpr_mean + tpr_std, 1)
      tprs_lower = np.maximum(tpr_mean - tpr_std, 0)
      print 'plt fill between'
      plt.fill_between(
        fpr_mean,
        tprs_lower,
        tprs_upper,
        color=colors[ik*3+ij],  # was ij
        alpha=0.15,
#         label=r"$\pm$ 1 std. dev.",
      )


plt.legend(loc="lower right")
plt.grid()
plt.tight_layout()

plt.xscale('log')
# plt.show()
xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,3,4,5,6,7,8,9])]+[1]
plt.plot(xy, xy, color='grey', linestyle='--')
plt.xlim([10**-3, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('background efficiency')
plt.ylabel('signal efficiency')

if len(years) == 1:
  plt.title(years[0])
else:
  plt.title('full dataset')
  
plt.tight_layout()
plt.savefig('results/comparisons/compare_roc_%s_%s_Average_noIP2D_forAN.pdf'%(kind, (len(years) == 1)*years[0] + (len(years) > 1)*"allYears" ) )
    
    
# plt.xlim([0.03, 0.1])
# plt.ylim([0.6, 1.0])
# plt.savefig('results/comparisons/compare_roc_%s_allYears_Average_IP_zoom.pdf'%(kind)  )
