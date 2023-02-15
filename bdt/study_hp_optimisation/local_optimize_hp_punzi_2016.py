import os, sys
# sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import pandas
import pickle
import sklearn
import joblib
from sklearn.ensemble  import GradientBoostingClassifier
from sklearn.metrics   import roc_curve
from sklearn.model_selection import train_test_split
from math  import sqrt 
import pdb
import uproot
from   xgboost import XGBClassifier, plot_importance
import numpy as np

from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt.domain_reduction import SequentialDomainReductionTransformer


mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

a_punzi = 2
b_punzi = 5

from os import path


def xgboost_fom(
                max_depth,
                learning_rate,
                n_estimators,
                min_child_weight,
                subsample,
                colsample_bytree,
              ):
              
    ## add scale for imbalanced datasets
    len_target0 = float(len(train[train.target == 0]))
    len_target1 = float(len(train[train.target == 1]))
    scale_pos_weight = float(len_target0/ len_target1  )      

    clf = XGBClassifier(
        max_depth        = int(round(max_depth)), 
        learning_rate    = learning_rate,
        n_estimators     = int(round(n_estimators)),
        min_child_weight = min_child_weight, 
        subsample        = subsample,
        colsample_bytree = colsample_bytree,
        seed             = 1986,
        verbosity        = 0,
        use_label_encoder = False,
        scale_pos_weight = scale_pos_weight
        )
    
    clf.fit(
        train[features], 
        train.target,
        eval_set              = [(train[features], train.target), (test[features], test.target)],
        early_stopping_rounds = 20,
        eval_metric           = 'auc', #'auc',
        verbose               = False,
        sample_weight         = train['normfactor'],
    )

    best_iter = clf.best_iteration
    true_bkg = test[(test.target==0) & (test.pass_preselection==1)]
    true_sig = test[(test.target==1) & (test.pass_preselection==1)]

    true_bkg['bdt_score']  = clf.predict_proba(true_bkg.loc[:, features_tuple], ntree_limit=best_iter)[:, 1] ## should be probability to be signal (bdt_score) for each data entry
    true_sig['bdt_score']  = clf.predict_proba(true_sig.loc[:, features_tuple], ntree_limit=best_iter)[:, 1] ## should be probability to be signal (bdt_score) for each data entry


    punzis = {}
    print ('score \t S (on test only) \t B (on test s.) \t punzi')
    for wp in np.arange(0.7,1,0.01):
        n_sig = float( len(true_sig[true_sig.bdt_score > wp])) 
        n_bkg = float( len(true_bkg[true_bkg.bdt_score > wp])) 
        den = pow(a_punzi,2)/8 + 9*pow(b_punzi,2)/13 + a_punzi*sqrt(n_bkg) + 0.5*b_punzi*sqrt(pow(b_punzi,2) + 4*a_punzi*sqrt(n_bkg) + 4*n_bkg)
        if den > 0:
            punzis[wp] =  n_sig/ den
        else:
            punzis[wp] = 0.0
        
        print ('%.3f\t'%wp, n_sig, '(%.3f)'%(float(n_sig/len(true_sig))), '\t', n_bkg, '(%.3f)'%(float(n_bkg/len(true_bkg))), '\t', punzis[wp])
    return max(punzis.values())


##########################################################################################
#####   SIGNAL AND BACKGROUND SELECTION
##########################################################################################
sig_mass_str = '(tagged_mass > {M} - 2.5*{S}) & (tagged_mass < {M} + 2.5*{S})'.format( M=mc_mass,S=mc_sigma)    

truth_match_str = '( (truthMatchMum==1) & (truthMatchMup ==1) & (truthMatchTrkm==1) & (truthMatchTrkp==1) )'

ct_str = '( ( (tagB0==1) & (genSignal==1)) | ( (tagB0==0) & (genSignal==2) ) )'

bkg_mass_str = '( ( (tagged_mass > {M}-7.*{S} ) & (tagged_mass < {M}-3.*{S} ) ) | ( (tagged_mass > {M}+3.*{S} ) & (tagged_mass < {M}+7*{S} ) ) )' .format( M=mc_mass,S=mc_sigma) 


sig_selection_cutbased = sig_mass_str + ' & ' + \
                         truth_match_str + ' & ' +\
                         ct_str + ' & ' +\
                         '(trig == 0) & ((mumuMass < 2.702) | (mumuMass > 4.))'

bkg_selection_cutbased = bkg_mass_str + ' & ' + \
                         '((mumuMass < 2.702) | (mumuMass > 4.)) '



sig_list = [
          'sub_samples/sample_2016_MC_LMNR_0_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_2_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_3_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_4_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_5_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_6_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_7_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_8_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_9_add_vars_passPreselection.root:ntuple',  
          'sub_samples/sample_2016_MC_LMNR_10_add_vars_passPreselection.root:ntuple', 
]

bkg_list = [
          'sub_samples/sample_2016_data_LMNR_0_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_2_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_3_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_4_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_5_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_6_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_7_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_8_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_9_add_vars.root:ntuple',
          'sub_samples/sample_2016_data_LMNR_10_add_vars.root:ntuple',
]


tag = '_optimize_punzi_'

##########################################################################################
#####   FEATURES AND BRANCHES
##########################################################################################
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

branches = features + [
    'tagged_mass',
    'pass_preselection',
]

branches = list(set(branches))

sig_uproot = uproot.concatenate(sig_list, 
                   expressions=branches + ['weight'],
                   cut = sig_selection_cutbased, 
                   library="np"
                   )

sig = pandas.DataFrame(sig_uproot)

bkg_uproot = uproot.concatenate(bkg_list, 
                   expressions=branches,
                   cut = bkg_selection_cutbased, 
                   library="np"
                   )
bkg = pandas.DataFrame(bkg_uproot)


##########################################################################################
#####   DEFINE THE TARGETS
##########################################################################################

sig['target'] = np.ones (sig.shape[0]).astype(int)
bkg['target'] = np.zeros(bkg.shape[0]).astype(int)

sig['normfactor'] = sig.weight
bkg['normfactor'] = 1. 


##########################################################################################
#####   SPLIT TRAIN AND TEST SAMPLES
##########################################################################################
data_all = pandas.concat([sig, bkg])

train, test = train_test_split(data_all, test_size=0.3, random_state = 17)

features_tuple = tuple(features)

par_bounds = {
  'max_depth': (3, 7),
  'learning_rate': (0.005, 0.1),
  'n_estimators': (100, 1000),
  'min_child_weight': (1, 7),
  'subsample': (0.5, 0.9),
  'colsample_bytree' :(0.55, 1.),
  }


bounds_transformer = SequentialDomainReductionTransformer(
     minimum_window=[
                      0.05, ## colsample_by_tree
                      0.01, ## learning rate
                      0.5,  ## max depth 
                      0.5,  ## min child w
                      0.5,  ## n estimators
                      0.05, ## subsample 
                     ])

optimizer = BayesianOptimization(xgboost_fom,
                                 par_bounds,
                                 verbose=2,
                                 random_state=12,
                                 bounds_transformer=bounds_transformer

                                )
## save optimization steps
logger = JSONLogger(path="logs_punzi_2016_fixBkg_noIP2D.json")

optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

my_init_points = 5
optimizer.maximize(
#     n_iter = 10,
    n_iter = 50, #50
    init_points = my_init_points,
    xi=1e-1,
    acq="ei",
#     alpha = 1e-1,
)

print('-'*53)

print('Final Results')
print(optimizer.max)

print ('\n ---- all iterations ------ \n')
for i, res in enumerate(optimizer.res):
    print("Iteration {}: \n\t{}".format(i, res))
    
    
for i, ipar in enumerate(optimizer.space.keys):  

  par_min_bound = [b[i][0] for b in bounds_transformer.bounds]
  par_max_bound = [b[i][1] for b in bounds_transformer.bounds]
  x = [x[i] for x in optimizer.space.params]
  bounds_transformers_iteration = list(range(my_init_points, len(x)))    

  plt.plot(bounds_transformers_iteration, par_min_bound[1:], label='%s lower bound'%ipar)
  plt.plot(bounds_transformers_iteration, par_max_bound[1:], label='%s upper bound'%ipar)
  plt.plot(x, label=ipar)
  plt.legend()
  plt.savefig('test_%s_bounds_6pars_punzi_2016_fixBkg_noIP2D.pdf'%ipar)
  plt.clf()

