import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("inputfile" , help = "", default = '../final_ntuples/2016data_noIP2D_addxcutvariable.root')
args = parser.parse_args()
year = args.year

import os, sys
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

import pandas, root_numpy
import pdb
import sklearn

import ROOT
import root_pandas


# ifile = '../final_ntuples/%sdata_noIP2D.root'%(args.year)  
ifile = args.inputfile        
print 'loading dataset %s'%ifile
dataset = pandas.DataFrame(
    root_numpy.root2array(
        ifile,
        'ntuple',
    )
)
print '\t...done'

dataset = dataset[(dataset.tagged_mass > 5) & (dataset.tagged_mass < 5.6) &  \
         (dataset.mumuMass*dataset.mumuMass > 8.68) & (dataset.mumuMass*dataset.mumuMass < 10.09) & \
         (dataset.pass_preselection == 1) & \
         (dataset.passB0Psi_jpsi == 1) & \
         (dataset.xcut == 0) #& \
#          (dataset.l1_00_OS == 1 ) & \
#          ((dataset.runN > 278808)  )
#          & (dataset.kstTrk1Pt > 1.2) 
         ] 

## for Psi bin
# dataset = dataset[(dataset.tagged_mass > 5) & (dataset.tagged_mass < 5.6) &  \
#          (dataset.mumuMass*dataset.mumuMass > 12.86) & (dataset.mumuMass*dataset.mumuMass < 14.18) & \
#          (dataset.pass_preselection == 1) & \
#          (dataset.passB0Psi_psip == 1)  ] 



# https://github.com/scikit-hep/root_pandas
dataset.to_root( ifile.replace('.root', '_passSPlotCuts.root'.replace('final_ntuples/', 'final_ntuples/for_splot/')), key='ntuple', store_index=False)
# dataset.to_root( ifile.replace('.root', '_passSPlotCuts_PsiP.root'), key='ntuple', store_index=False)
