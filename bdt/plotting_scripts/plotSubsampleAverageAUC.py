import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--tag"      , help = "", default = 'punzi_removeTkMu_fixBkg_2017')
args = parser.parse_args()
year = args.year

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

years = [ year ]
nsamples = 11

tags = [
#     '_punzi_removeTkMu_fixBkg',
    ## isolation
#     '_useMaxIso',
    ## kstar mass
#     '_deltakstMass',
#     '_removekstMass',
#     ## IP pars
#     '_useMinMaxIP',
    '_noIP2D',
#     '_noIP2D_noTrkSign'
] 


tag_labels = [
#     'nominal',
    ## isolation
#     'max iso',
    ## kstar mass
#     'use |m-0.892|',
#     'remove kst mass',
# #     ## IP pars	
#     'min/max IP',
    'no trkMinIP',
#     'no trkMinIP, newHP 6then3',
#     'no trkMinIP, no Trk Sign',
] 
kinds = ['test', 'train']


metrics = ['xaxis', 'auc']#, 'tpr', 'thresholds']
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


## consider mass preselection which is different in 2016 and 2018
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
    auc = []
    auc_train = []
    auc_xaxis = []
    
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

        auc.append(       clfs[ii].evals_result()['validation_1']['auc'])#[:, 1])
        auc_train.append( clfs[ii].evals_result()['validation_0']['auc'] )#[:, 1])
        auc_xaxis.append( np.arange(len(auc_train[ii])))  ## to be cross checked


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
    delta_train_test = []
    
    for ii, itag in enumerate(tags):
        
        print 'saving train/test results for tag=%s for sample %s'%(itag,isample) #+ '  auc = %s'%roc_auc_score(test.target, pred[ii])

        train_arr = np.array(auc_train[ii])
        test_arr = np.array(auc[ii])
        x_arr = np.array(auc_xaxis[ii])
        x_arr = x_arr / float(len(x_arr))
        
        delta_train_test.append( list (np.subtract(train_arr, test_arr)) )

        results_tags[itag]['test']['auc'].append( test_arr)
        results_tags[itag]['test']['xaxis'].append(x_arr)

        results_tags[itag]['train']['auc'].append( train_arr)
        results_tags[itag]['train']['xaxis'].append(x_arr)


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
          fpr           = results[kind]['xaxis'][i]
          tpr           = results[kind]['auc'][i]
          interp_tpr    = np.interp(fpr_mean, fpr, tpr)
          interp_tprs.append(interp_tpr)

      tpr_mean     = np.mean(interp_tprs, axis=0)

      print 'plot'
      plt.plot(
        fpr_mean,
        tpr_mean,
        color=colors[ik*3+ij],
        label="Ave AUC %s"%(kind),#(AUC = %0.4f $\pm$ %0.4f) %s %s" % (auc, std_auc, tag_labels[ij], kind),
        lw=2,
        alpha=0.8,
      )

plt.legend(loc="lower right")
plt.grid()
plt.tight_layout()

plt.xlim([10**-3, 1.0])
# plt.ylim([0.0, 0.006])
plt.xlabel('% of n trees')
plt.ylabel('AUC')

if len(years) == 1:
  plt.title(years[0])
else:
  plt.title('full dataset')
  
plt.tight_layout()
plt.savefig('results/comparisons/compare_AUC_%s_%s_newNtuples_Average.pdf'%(kind, (len(years) == 1)*years[0] + (len(years) > 1)*"allYears" ) )
    