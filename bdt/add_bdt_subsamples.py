import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()
year = args.year

import os, sys
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

# sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')
import pandas, root_numpy

import sklearn
from sklearn.externals import joblib

import ROOT, pdb
import root_pandas

samples = [
           'data_LMNR',
           'data_Charmonium',
           'MC_LMNR', 
           'MC_JPSI', 
           'MC_PSI',
#            'MC_otherb' 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_BJPSIK',
#            'MC_HBJPSIX',    # 2018
#              'MC_B0ZC', ##2018
#            'data_sameSign', 
          ]


tags = [
#             '_punzi_removeTkMu_fixBkg_%s'%year,
#             '_noIP2D_noTrkSign_%s'%year,
            '_noIP2D_%s'%year,
#             '_noIP2D_noNan_%s'%year,
       ] 

tag = tags[0] ### remember to update pass_preselection

for str_file in samples:
    for i in range(11):  
#         pdb.set_trace()
        ifile = 'sub_samples/sample_%s_%s_%s_add_vars.root'%(args.year, str_file, str(i))  
        if 'MC' in str_file:
          ifile = 'sub_samples/sample_%s_%s_%s_add_vars_passPreselection.root'%(args.year, str_file, str(i))  

        classifier_name = 'results/classifier_%s_%s.pkl' %(tag,i)
        print 'adding bdt score from classifier: %s'%(classifier_name)
        classifier = joblib.load(classifier_name)
        
        feat_names = [
         'bCosAlphaBS',
         'bLBS/bLBSE',
         'trkmDCASign',
         'trkpDCASign',
         'bVtxCL',
         'bDCABS/bDCABSE',
         'kstarmass',
         'sum_isopt_04',
        ]
    
        additional = [
          'tagged_mass',
          'pass_preselection',
        ]

        print 'loading support dataset %s' %ifile
        dataset_support = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
                branches=feat_names + additional,
            )
        )
        print '\t...done'
        
        print 'loading dataset...'
        dataset = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
            )
        )
        print '\t...done'
        

        ## add branch for BDT 
        print 'computing probabilities...'
        bdt_prob_array = classifier.predict_proba(dataset_support[feat_names])[:,1]
        print '\t...done'
        
        print 'adding new column to the dataset...'
        dataset['bdt_prob'] = bdt_prob_array
        print '\t...done'	
        
        # https://github.com/scikit-hep/root_pandas
        dataset.to_root( ifile.replace('.root', '_addBDT%s.root'%tag), key='ntuple', store_index=False)
