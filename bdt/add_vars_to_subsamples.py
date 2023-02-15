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
import ROOT
import pandas, root_numpy, root_pandas
import pdb


samples = [
           'data_LMNR',
           'data_Charmonium',
           'MC_JPSI', 
           'MC_LMNR', 
           'MC_PSI',
          ]


def addTrk1DCA (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkpDCABS']
   else: return row['kstTrkmDCABS']
def addTrk1DCAE (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']:      return row['kstTrkpDCABSE']
   else: return row['kstTrkmDCABSE']
def addTrk1Eta (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkpEta']
   else: return row['kstTrkmEta']
def addTrk1Phi (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkpPhi']
   else: return row['kstTrkmPhi']
def addTrk1Pt (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkpPt']
   else: return row['kstTrkmPt']

def addTrk2DCA (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkmDCABS']
   else: return row['kstTrkpDCABS']
def addTrk2DCAE (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']:      return row['kstTrkmDCABSE']
   else: return row['kstTrkpDCABSE']
def addTrk2Eta (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkmEta']
   else: return row['kstTrkpEta']
def addTrk2Phi (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkmPhi']
   else: return row['kstTrkpPhi']
def addTrk2Pt (row):
   if row['kstTrkpPt'] > row['kstTrkmPt']: return row['kstTrkmPt']
   else: return row['kstTrkpPt']


def addMu1Eta (row):
   if row['mupPt'] > row['mumPt']: return row['mupEta']
   else: return row['mumEta']
def addMu1Pt (row):
   if row['mupPt'] > row['mumPt']: return row['mupPt']
   else: return row['mumPt']
def addMu1Phi (row):
   if row['mupPt'] > row['mumPt']: return row['mupPhi']
   else: return row['mumPhi']
def addMu2Eta (row):
   if row['mupPt'] > row['mumPt']: return row['mumEta']
   else: return row['mupEta']
def addMu2Phi (row):
   if row['mupPt'] > row['mumPt']: return row['mumPhi']
   else: return row['mupPhi']
def addMu2Pt (row):
   if row['mupPt'] > row['mumPt']: return row['mumPt']
   else: return row['mupPt']


def addTaggedMass (row):
   if row['tagB0'] == 1: return row['bMass']
   else: return row['bBarMass']
   

for str_file in samples:
    for i in range(11):  
        
        ifile = 'sub_samples/sample_%s_%s_%s_newphi.root'%(args.year, str_file, str(i))  
        
        print 'loading dataset %s ...'%ifile
        dataset = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
            )
        )
        print '\t...done'
        
        dataset['tagged_mass']     = dataset.apply (lambda row: addTaggedMass(row), axis=1)
        dataset = dataset[(dataset.tagged_mass > 4.5)] 

        ## get distributions for leading/trailing trks
        dataset['kstTrk1DCABS']  = dataset.apply (lambda row: addTrk1DCA(row),  axis=1)
        dataset['kstTrk1DCABSE'] = dataset.apply (lambda row: addTrk1DCAE(row), axis=1)
        dataset['kstTrk1Eta']    = dataset.apply (lambda row: addTrk1Eta(row), axis=1)
        dataset['kstTrk1Phi']    = dataset.apply (lambda row: addTrk1Phi(row), axis=1)
        dataset['kstTrk1Pt']     = dataset.apply (lambda row: addTrk1Pt(row), axis=1)

        dataset['kstTrk2DCABS']  = dataset.apply (lambda row: addTrk2DCA(row),  axis=1)
        dataset['kstTrk2DCABSE'] = dataset.apply (lambda row: addTrk2DCAE(row), axis=1)
        dataset['kstTrk2Eta']    = dataset.apply (lambda row: addTrk2Eta(row), axis=1)
        dataset['kstTrk2Phi']    = dataset.apply (lambda row: addTrk2Phi(row), axis=1)
        dataset['kstTrk2Pt']     = dataset.apply (lambda row: addTrk2Pt(row), axis=1)

        dataset['mu1Eta']    = dataset.apply (lambda row: addMu1Eta(row), axis=1)
        dataset['mu1Phi']    = dataset.apply (lambda row: addMu1Phi(row), axis=1)
        dataset['mu1Pt']     = dataset.apply (lambda row: addMu1Pt(row),  axis=1)
        dataset['mu2Eta']    = dataset.apply (lambda row: addMu2Eta(row), axis=1)
        dataset['mu2Phi']    = dataset.apply (lambda row: addMu2Phi(row), axis=1)
        dataset['mu2Pt']     = dataset.apply (lambda row: addMu2Pt(row),  axis=1)
        

        dataset['trkpDCASign']  = abs(dataset.kstTrkpDCABS/dataset.kstTrkpDCABSE)
        dataset['trkmDCASign']  = abs(dataset.kstTrkmDCABS/dataset.kstTrkmDCABSE)

        dataset['min_trkMinIP2D'] = dataset[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].min(axis=1) 
        dataset['max_trkMinIP2D'] = dataset[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].max(axis=1) 

        ## define isolation: # tracks with pt in a cone
        dataset['isopt_mum_04']  = dataset.mumIsoPt_dr04    /dataset.mumPt
        dataset['isopt_mup_04']  = dataset.mupIsoPt_dr04    /dataset.mupPt
        dataset['isopt_trkm_04'] = dataset.kstTrkmIsoPt_dr04/dataset.kstTrkmPt
        dataset['isopt_trkp_04'] = dataset.kstTrkpIsoPt_dr04/dataset.kstTrkpPt
        dataset['sum_isopt_04']  = dataset.isopt_mum_04 + dataset.isopt_mup_04 + dataset.isopt_trkm_04 + dataset.isopt_trkp_04

        ## define tagged kstar mass
        dataset['kstarmass']       = dataset.tagB0*dataset.kstMass +(1- dataset.tagB0)*dataset.kstBarMass


        if args.year == '2016':
            dataset['pass_preselection'] =  ( dataset.mumNTrkLayers >= 6)  & ( dataset.mupNTrkLayers >= 6 ) & \
                                            ( dataset.mumNPixLayers >= 1)  & ( dataset.mupNPixLayers >= 1 ) & \
                                            ( dataset.mumHighPurity == 1 ) & ( dataset.mupHighPurity == 1 ) & \
                                            ( dataset.mumTMOneStationTight == 1 ) & ( dataset.mupTMOneStationTight == 1 ) & \
                                            ( dataset.kstTrkmHighPurity == 1 )    & ( dataset.kstTrkpHighPurity == 1    ) & \
                                            ( dataset.kkMass > 1.035 ) & \
                                            ( dataset.kstTrkmTrackerMuon == 0 ) & ( dataset.kstTrkpTrackerMuon == 0 ) & \
                                            (~((dataset.kstTrkmGlobalMuon == 1) & ( dataset.kstTrkmNTrkLayers > 5 ) & ( dataset.kstTrkmNPixHits > 0))) & \
                                            (~((dataset.kstTrkpGlobalMuon == 1) & ( dataset.kstTrkpNTrkLayers > 5 ) & ( dataset.kstTrkpNPixHits > 0))) 
        
        ## for 2017 and 2018, adding significance of the track wrt BS (required at trigger level)
        else:
            dataset['pass_preselection'] = ( dataset.mumTMOneStationTight == 1 ) & ( dataset.mupTMOneStationTight == 1 ) & \
                                            ( dataset.kkMass > 1.035 ) & \
                                            ( dataset.kstTrkmTrackerMuon == 0 ) & ( dataset.kstTrkpTrackerMuon == 0 ) & \
                                            (~((dataset.kstTrkmGlobalMuon == 1) & ( dataset.kstTrkmNTrkLayers > 5 ) & ( dataset.kstTrkmNPixHits > 0))) & \
                                            (~((dataset.kstTrkpGlobalMuon == 1) & ( dataset.kstTrkpNTrkLayers > 5 ) & ( dataset.kstTrkpNPixHits > 0))) & \
                                            (( (dataset.charge_trig_matched ==  1) & (dataset.kstTrkpPt > 1.2) & (dataset.trkpDCASign > 2) ) | \
                                             ( (dataset.charge_trig_matched == -1) & (dataset.kstTrkmPt > 1.2) & (dataset.trkmDCASign > 2) ) )
                
        
        dataset['pass_preselection'] = dataset['pass_preselection'].astype(np.int32)
        dataset['eventD'] = dataset['eventN'].astype(np.double)
        
        if 'data' in str_file:
            dataset = dataset[(dataset.pass_preselection > 0)] 

        if 'MC' in str_file:
          dataset['weight']   = dataset['weight'].astype(np.float32)
          if year =='2016': 
            dataset['weightBF'] = dataset['weightBF'].astype(np.float32)
            dataset['weightGH'] = dataset['weightGH'].astype(np.float32)
          
        # https://github.com/scikit-hep/root_pandas
        dataset.to_root( ifile.replace('.root', 'add_vars.root').replace('newphi', ''), key='ntuple', store_index=False)
