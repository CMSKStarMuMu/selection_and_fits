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

import ROOT
import pandas, root_numpy, root_pandas

samples = [
           'MC_JPSI', 
           'MC_LMNR', 
           'MC_PSI',
          ]


for str_file in samples:
    for i in range(11):  
      
        ifile = 'sub_samples/sample_%s_%s_%s_add_vars.root'%(args.year, str_file, str(i))  
          
        print 'loading dataset %s'%ifile
        dataset = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
            )
        )
        print '\t...done'
        
        dataset = dataset[(dataset.pass_preselection > 0) & (dataset.tagged_mass > 4.5)] 
        
        # https://github.com/scikit-hep/root_pandas
        dataset.to_root( ifile.replace('.root', '_passPreselection.root'), key='ntuple', store_index=False)
