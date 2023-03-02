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

##preselection
# ifile = '../bdt/ntuples_pre_selection/%sData_passPreselection_passSPlotCuts.root'%year
# w_data = '../../../compare_data_mc/preBDT/out_distribution_JPSI_2018A_2018D_L1all_preBDT_Nov21_pol2.root'
# if year =='2017':
#   w_data = '../../../compare_data_mc/preBDT/out_distribution_JPSI_2017B_2017F_L1all_preBDT_Nov21_pol2.root'
# if year =='2016':
#   w_data = '../../../compare_data_mc/preBDT/out_distribution_JPSI_2016B_2016H_L1all_preBDT_Nov21_pol2.root'

## post selection
ifile = args.inputfile.replace('.root', '_passSPlotCuts.root')

w_data = 'out_distribution_JPSI_2018A_2018D_L1all_postBDT_Nov21.root'
if year =='2017':
  w_data = 'out_distribution_JPSI_2017B_2017F_L1all_postBDT_Nov21.root'
if year =='2016':
  w_data = 'out_distribution_JPSI_2016B_2016H_L1all_postBDT_Nov21.root'

## for PSI bin
# w_data = 'out_distribution_PSI_2018A_2018D_L1all_postBDT_Nov21.root'
# if year =='2017':
#   w_data = 'out_distribution_PSI_2017B_2017F_L1all_postBDT_Nov21.root'
# if year =='2016':
#   w_data = 'out_distribution_PSI_2016B_2016H_L1all_postBDT_Nov21.root'

print 'loading dataset %s'%ifile
dataset = pandas.DataFrame(
    root_numpy.root2array(
        ifile,
        'ntuple',
    )
)
print '\t...done'

print 'loading support dataset %s'%w_data

dataset_support = pandas.DataFrame(
    root_numpy.root2array(
        w_data,
        'fulldata',
    )
)
print '\t...done'
print '\t...have n events: ', len(dataset)

dataset = dataset.merge(dataset_support,how='left',on=['tagged_mass', 'eventD'], suffixes=('', '_remove'))
dataset.drop([icol for icol in dataset.columns if 'remove' in icol], axis=1, inplace=True)
print '\t...have n events after merge: ', len(dataset)

# https://github.com/scikit-hep/root_pandas
dataset.to_root( ifile.replace('.root', '_mergeSweights.root'), key='ntuple', store_index=False)
