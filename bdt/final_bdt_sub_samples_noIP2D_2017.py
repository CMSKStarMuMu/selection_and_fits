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

from   xgboost import XGBClassifier, plot_importance

from ROOT      import TFile, TTree, TH1F, gROOT, TChain


mc_sigma = 0.0400 ## was 49
mc_mass  = 5.27783 

gROOT.SetBatch(True)

## weight input background mass
min_m  = 4.8
max_m  = 5.8
n_bins = 500


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
                         '(trig == 0) & ((mumuMass < 2.702) )'

bkg_selection_cutbased = bkg_mass_str + ' & ' + \
                         '((mumuMass < 2.702) ) '



sig_list = [
          'sub_samples/sample_2017_MC_LMNR_0_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_1_add_vars_passPreselection.root',
          'sub_samples/sample_2017_MC_LMNR_2_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_3_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_4_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_5_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_6_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_7_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_8_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_9_add_vars_passPreselection.root', 
          'sub_samples/sample_2017_MC_LMNR_10_add_vars_passPreselection.root',
]

bkg_list = [
          'sub_samples/sample_2017_data_LMNR_0_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_1_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_2_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_3_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_4_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_5_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_6_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_7_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_8_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_9_add_vars.root',
          'sub_samples/sample_2017_data_LMNR_10_add_vars.root',
]

import pdb
for isample in range(0,11):

#     tag = '_useSimpleIP_test_%s'%isample
#     tag = '_useMinMaxIP_2017_%s'%isample
#     tag = '_useMaxIso_2017_%s'%isample
    tag = '_noIP2D_noNan_2017_%s'%isample
#     tag = '_removekstMass_2017_%s'%isample
#     tag = '_deltakstMass_2017_%s'%isample

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
    
    sig = pandas.DataFrame(
        root_numpy.root2array(
#             sig_list[isample], 
            sig_list[:isample] + sig_list[(isample + 1):], 
            'ntuple',
            branches  = branches + ['weight'],
            selection = sig_selection_cutbased,
        )
    )
    
    bkg = pandas.DataFrame(
        root_numpy.root2array(
#             bkg_list[isample], 
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
    
    ## add b mass calculation for plotting only
    bkg.hist(column='tagged_mass',bins=400)
#     plt.savefig('mass_bkg_onesample.pdf')
    sig.hist(column='tagged_mass',bins=400)   
#     plt.savefig('mass_sig_onesample.pdf')
    
    sig['normfactor'] = sig.weight#*sig.m_weight #/len(sig) *10000
    bkg['normfactor'] = 1. ##bkg.m_weight# /len(bkg) #* 10000
    
    ##########################################################################################
    #####   SPLIT TRAIN AND TEST SAMPLES
    ##########################################################################################
    data_all = pandas.concat([sig, bkg])
    data_all = data_all.dropna(subset=['bDCABS/bDCABSE','bLBS/bLBSE'])
    data = data_all[data_all.pass_preselection==1]
#     pdb.set_trace()
#     exit()
    train, test = train_test_split(data, test_size=0.30, random_state = 17)

    ## add scale for imbalanced datasets
    len_target0 = float(len(train[train.target == 0]))
    len_target1 = float(len(train[train.target == 1]))
    scale_pos_weight = float(len_target0/ len_target1  )      

    # Final Results
    # {'target': 1229.8436579065483, 'params': {'colsample_bytree': 0.9652754948134503, 'learning_rate': 0.09776091659728355, 
    # 'max_depth': 3.5133750725542874, 'min_child_weight': 4.210391142964952, 'n_estimators': 715.5798160205072, 'subsample': 0.6387998087738809}}
    ##limining n estimators: 'n_estimators': 510.3865440136714

    clf = XGBClassifier(
        colsample_bytree = 0.965,
        learning_rate    = 0.098,
        min_child_weight = 4.210,
        n_estimators     = 700,
        subsample        = 0.638,
        max_depth        = int(round(3.5133)),
        seed             = 1986,
        silent           = False,
        scale_pos_weight = scale_pos_weight
        )


    clf.fit(
        train[features], 
        train.target,
        eval_set              = [(train[features], train.target), (test[features], test.target)],
        early_stopping_rounds = 20,
        eval_metric           = 'auc',
        verbose               = True,
        sample_weight         = train['normfactor'],
    )
    
    best_iter = clf.best_iteration
    joblib.dump(clf, 'results/classifier_%s.pkl' %(tag), compress=True)
    
    ##########################################################################################
    # #####   PREDICT ON THE TEST UNBIAS SAMPLE
    # ##########################################################################################
    pred = clf.predict_proba(test[features],  ntree_limit=best_iter)[:, 1] #### !!! changed from test_unbias
#     limit_pred = clf.predict_proba(test[features], ntree_limit=best_iter)[:, 1] 
    
    train_sig = clf.predict_proba( train[features][train.target>0.5], ntree_limit=best_iter )[:,1]
    train_bkg = clf.predict_proba( train[features][train.target<0.5], ntree_limit=best_iter )[:,1]
    
    test_sig = clf.predict_proba(test[features][test.target>0.5], ntree_limit=best_iter)[:,1]
    test_bkg = clf.predict_proba(test[features][test.target<0.5], ntree_limit=best_iter)[:,1]
    
    ##########################################################################################
    #####   ROC CURVE
    ##########################################################################################
    plt.clf()
    
    import itertools
    xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.xlim([10**-3, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    
    ## draw baseline point
    #   plt.plot([fpr], [tpr], label='baseline', markerfacecolor='red', marker='o', markersize=10)
    # 
    # ## draw ROC         
    fpr, tpr, threshold = roc_curve(test.target, (pred * test.pass_preselection) ) ### !! changed from test unbias
    plt.plot(fpr, tpr, color='r', label='ROC test')
    
    pred2 = clf.predict_proba(train[features], ntree_limit=best_iter)[:, 1]
    fpr2, tpr2, sara2 = roc_curve(train.target, pred2)
    plt.plot(fpr2, tpr2, color='b', label='ROC train')

    plt.xscale('log')
    plt.grid()
    
    plt.legend(loc='best')
    plt.grid()
    plt.title('ROC')
    plt.tight_layout()
    plt.savefig('results/roc_train_test_%s.pdf' %(tag))
    plt.clf()

    roc_file = open('results/roc_%s.pck' %(tag), 'w+')
    pickle.dump((tpr, fpr), roc_file)
    roc_file.close()
    
    ##########################################################################################
    #####   OVERTRAINING TEST
    ##########################################################################################
    
    low  = 0
    high = 1
    low_high = (low,high)
    bins = 50
    
    #################################################
    plt.hist(
        train_sig,
        color='r', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='S (train)'
    )
    
    #################################################
    plt.hist(
        train_bkg,
        color='b', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='B (train)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_sig,
        bins=bins, 
        range=low_high, 
        normed=True,
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_sig) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='r', 
        label='S (test)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_bkg,
        bins=bins, 
        range=low_high, 
        normed=True
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_bkg) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='b', 
        label='B (test)'
    )
    
    #################################################
    plt.xlabel('BDT output')
    plt.ylabel('Arbitrary units')
    plt.legend(loc='best')
    ks_sig = ks_2samp(train_sig, test_sig)
    ks_bkg = ks_2samp(train_bkg, test_bkg)
    plt.suptitle('KS p-value: sig = %.3f%s - bkg = %.2f%s' %(ks_sig.pvalue * 100., '%', ks_bkg.pvalue * 100., '%'))
    
    # plt.tight_layout()
    plt.savefig('results/overtrain_%s.pdf' %(tag))
    
    ##########################################################################################
    #####   FEATURE IMPORTANCE
    ##########################################################################################
    plot_importance(clf)
    plt.tight_layout()
    plt.savefig('results/feat_importance_%s.pdf' %(tag))
    
    plt.close()        
    
    ##########################################################################################
    #####   OVERTRAINING SCORE
    ##########################################################################################
    plt.clf()
    
    auc_train = clf.evals_result()['validation_0']['auc']
    auc_test  = clf.evals_result()['validation_1']['auc']
    
    n_estimators = np.arange(len(auc_train))
    
    plt.plot(n_estimators, auc_train, color='r', label='auc train')
    plt.plot(n_estimators, auc_test , color='b', label='auc test' )
    
    plt.xlabel('# tree')
    plt.ylabel('Area Under ROC')
    
    plt.xscale('log')
    plt.grid()
    
    # plt.xlim([1, 1000])
    # plt.ylim([0.985, 1.0])
    
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('results/auc_score_%s.pdf' %(tag))
    
    
    
    ######################################
    # Correlation Matrix Plot
    ######################################
    sig  ['bdt'] = clf.predict_proba(sig[features], ntree_limit=best_iter)[:, 1] 
    bkg  ['bdt'] = clf.predict_proba(bkg[features], ntree_limit=best_iter)[:, 1] 

    sig_corr = sig[features+ ['bdt', 'tagged_mass']]
    bkg_corr = bkg[features+ ['bdt', 'tagged_mass']]
    
    sig_correlations = sig_corr.corr()
    bkg_correlations = bkg_corr.corr()
    
    fig = plt.figure()
    ticks = np.arange(0,len(features)+2,1)
    
    ax  = fig.add_subplot(121)
    cax = ax.matshow(sig_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 9))


    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(features+ ['bdt', 'tagged_mass'],fontsize=5)
    ax.set_yticklabels(features+ ['bdt', 'tagged_mass'],fontsize=5)
    plt.xticks(rotation=90)
    
    # plot correlation matrix
    ax2  = fig.add_subplot(122)
    cax2 = ax2.matshow(bkg_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 11))
    ax2.set_xticks(ticks)
    ax2.set_yticks(ticks)
    ax2.set_xticklabels(features+ ['bdt', 'tagged_mass'],fontsize=5)
    ax2.set_yticklabels(features+ ['bdt', 'tagged_mass'],fontsize=0)
    plt.xticks(rotation=90)
    
    cbar = fig.colorbar(cax, orientation="horizontal", pad=0.045)
    cbar.ax.tick_params(labelsize=5)
    # cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
    # cb = plt.colorbar(cax, cax = cbaxes)  
    # fig.colorbar(cax)
    # fig.colorbar(cax2)
    
    plt.show()
    plt.savefig('results/bkg_correlation_%s.pdf' %(tag))
    
    plt.close()
    