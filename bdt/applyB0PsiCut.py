import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--tag"     , help = "", default = '_noIP2D')
args = parser.parse_args()
year = args.year

import os, sys
import ROOT, math
import numpy as np
import pandas, root_numpy
from ROOT import TLorentzVector
from copy import deepcopy as dc

from PhysicsTools.HeppyCore.utils.deltar import deltaR
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath('../utils/utils.py') ) ) )
from utils.utils import *

tag = args.tag + '_' + year

samples = [
           'data',
           'MC_LMNR', 
           'MC_JPSI', 
           'MC_PSI', 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_BJPSIK',
#            'MC_B0ZC', ##2018
#            'MC_HBJPSIX'
          ]

muonmass_ = 0.1056583745
kaonmass_ = 0.493677
pionmass_ = 0.139570
JPsiMass_ = 3.096916

tkp_lv = TLorentzVector()
tkm_lv = TLorentzVector()
mum_lv = TLorentzVector()
mup_lv = TLorentzVector()

pion_lv = TLorentzVector()
kaon_lv = TLorentzVector()

@np.vectorize
def addRejectPsi(
            mumuMass, mumuMassE, deltaB0M, deltaJpsiM, deltaPsiPM
          ):

    passSel = ( (abs(mumuMass - JPsiMass_) > 3*mumuMassE) & \
                (abs(mumuMass - PsiPMass_) > 3*mumuMassE) &  \
           (( (mumuMass < JPsiMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.19)) ) | \
            ( (mumuMass > PsiPMass_) & ~( (abs(deltaB0M - deltaPsiPM) < 0.08)) ) | \
            ( (mumuMass > JPsiMass_) & (mumuMass < PsiPMass_) & \
              ~( (abs(deltaB0M - deltaJpsiM) < 0.09) | (abs(deltaB0M - deltaPsiPM) < 0.07)) ))) 

    return passSel


def addRejectPsi2016(
            mumuMass, mumuMassE, deltaB0M, deltaJpsiM, deltaPsiPM
          ):

    passSel = ( (abs(mumuMass - JPsiMass_) > 3*mumuMassE) & \
                (abs(mumuMass - PsiPMass_) > 3*mumuMassE) &  \
           (( (mumuMass < JPsiMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.20) | (abs(deltaB0M - deltaPsiPM) < 0.0 )) ) | \
            ( (mumuMass > PsiPMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.0)  | (abs(deltaB0M - deltaPsiPM) < 0.11)) ) | \
            ( (mumuMass > JPsiMass_) & (mumuMass < PsiPMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.10) | (abs(deltaB0M - deltaPsiPM) < 0.08)) ))) 

    return passSel

# was 190
@np.vectorize
def addKeepJpsi(
            mumuMass, mumuMassE
          ):

    passSel = (abs(mumuMass - JPsiMass_) < 3*mumuMassE)
    return passSel

@np.vectorize
def addKeepPsip(
            mumuMass, mumuMassE
          ):

    passSel = (abs(mumuMass - PsiPMass_) < 3*mumuMassE)
    return passSel


@np.vectorize
def addDR(
            mumEta,  mumPhi,  
            tkpEta,  tkpPhi
          ):
        
    if mumEta == -99:
        return -99    

    return deltaR(mumEta, mumPhi, tkpEta, tkpPhi )

@np.vectorize
def addPsi2sMass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    tkm_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (mum_lv + mup_lv + tkp_lv + tkm_lv).M()
    return opt1

@np.vectorize
def addpipiMass(
            tkmPt,  tkmEta,  tkmPhi,  
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if tkmPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)

    opt1 = (mum_lv + tkp_lv ).M()
    return opt1


@np.vectorize
def addpiKMass(
            mumPt,  mumEta,  mumPhi,  
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, pionmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)

    opt1 = (mum_lv + tkp_lv ).M()
    return opt1


@np.vectorize
def addmmpi2(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi,
            tagB0
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    if tagB0 == 1:
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
    else:
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)

    opt1 = (mum_lv + mup_lv + pion_lv ).M()
    opt2 = (mum_lv + mup_lv + kaon_lv ).M()
    return opt1, opt2


@np.vectorize
def addmmpi2Paolo(
            mumuPt,  mumuEta,  mumuPhi, mumuMass, 
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi,
            tagB0
          ):
        
    mum_lv.SetPtEtaPhiM(mumuPt, mumuEta, mumuPhi, mumuMass)
    if tagB0 == 1:
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
    else:
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)

    opt1 = (mum_lv  + pion_lv ).M()
    opt2 = (mum_lv  + kaon_lv ).M()
    return opt1, opt2


@np.vectorize
def addkstarmass(
            tkmPt, tkmEta, tkmPhi,  
            tkpPt, tkpEta, tkpPhi,
            tagB0
          ):
        
    
    if tagB0 == 1:
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    else:
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (kaon_lv + pion_lv ).M()
    return opt1

@np.vectorize
def addbwtmass(
            mumuPt,  mumuEta,  mumuPhi, mumuMass, 
            tkmPt, tkmEta, tkmPhi,  
            tkpPt, tkpEta, tkpPhi,
            tagB0
          ):
        
    
    mum_lv.SetPtEtaPhiM(mumuPt, mumuEta, mumuPhi, mumuMass)
    if tagB0 == 1:
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    else:
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (mum_lv + kaon_lv + pion_lv ).M()
    return opt1




def addmmpiKaon (row):
   if row['tagB0'] == 1 :
      return row['bBarMass']
   else:      
      return row['bMass']
def addkst2 (row):
   if row['tagB0'] == 1 :
      return row['kstBarMass']
   else:      
      return row['kstMass']

def addpi1Pt (row):
   if row['tagB0'] == 1 :
      return row['kstTrkpPt']
   else:      
      return row['kstTrkmPt']
def addpi2Pt (row):
   if row['tagB0'] == 1 :
      return row['kstTrkmPt']
   else:      
      return row['kstTrkpPt']



@np.vectorize
def addmmkkmass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
    kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)

    opt1 = (mum_lv + mup_lv + pion_lv + kaon_lv).M()
    return opt1


@np.vectorize
def addmmpipimass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
    kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)

    opt1 = (mum_lv + mup_lv + pion_lv + kaon_lv).M()
    return opt1


@np.vectorize
def addxcut(wt_mass,  wt_kstarmass,  kaonPt, pionPt, mmpiMass, mmkMass):
        
    bool1 =  ( (5.2791 - wt_mass) - 0.3 ) / (-0.1-0.3)< (((wt_kstarmass-0.896)--0.4) / (0.6--0.4))
    bool2 = kaonPt > pionPt
    bool3 = (wt_kstarmass-0.896)>0
    bool4 = (mmpiMass > 3.2) & (mmpiMass < 3.6)
    bool5 = (mmkMass >  4.7) & (mmkMass  < 4.9)
    bool6 = ((mmkMass - 3.8) / (4.8 - 3.8)) > ((mmpiMass-3)/(3.6-3))
    
    xcut = bool1 & bool2 & bool3 & bool4 & bool5 & bool6
    return xcut



nSigma_psiRej = 3.

for str_file in samples:

    input_files = []
    print (str_file)

    for i in range(11):
        ofile = '../final_ntuples/%s%s%s_part%s_addB0Psi.root'%(year, str_file, tag, i)

        input_files = []
        input_files.append('../final_ntuples/%s%s%s_passBDT_part%s.root'%(args.year, str_file, tag, i ))        

        print ('loading dataset %s'%input_files[0])
        dataset_all = pandas.DataFrame(
            root_numpy.root2array(
                input_files,
                'ntuple',
            )
        )
        print ('\t...done. n events: ', len(dataset_all))
    
    
        dataset_all['deltaB0M']   = dataset_all.tagged_mass - B0Mass_  
        dataset_all['deltaJpsiM'] = dataset_all.mumuMass - JPsiMass_   
        dataset_all['deltaPsiPM'] = dataset_all.mumuMass - PsiPMass_   
    
        dataset_all['deltaM_jpsi'] = dataset_all.deltaB0M - dataset_all.deltaJpsiM
        dataset_all['deltaM_psi']  = dataset_all.deltaB0M - dataset_all.deltaPsiPM
    
        dataset_all['passB0Psi_jpsi'] = addKeepJpsi (dataset_all.mumuMass, dataset_all.mumuMassE)
        dataset_all['passB0Psi_psip'] = addKeepPsip (dataset_all.mumuMass, dataset_all.mumuMassE)


    
        if year == '2016':
            dataset_all['passB0Psi_lmnr'] = addRejectPsi2016(dataset_all.mumuMass, dataset_all.mumuMassE, dataset_all.deltaB0M, dataset_all.deltaJpsiM, dataset_all.deltaPsiPM )
        else:
            dataset_all['passB0Psi_lmnr'] = addRejectPsi(dataset_all.mumuMass, dataset_all.mumuMassE, dataset_all.deltaB0M, dataset_all.deltaJpsiM, dataset_all.deltaPsiPM )

        dataset_all['passB0Psi_jpsi'] = dataset_all['passB0Psi_jpsi'].astype(np.int32)
        dataset_all['passB0Psi_psip'] = dataset_all['passB0Psi_psip'].astype(np.int32)
        dataset_all['passB0Psi_lmnr'] = dataset_all['passB0Psi_lmnr'].astype(np.int32)

        ### add more variables
        dataset_all['wt_mass']= dataset_all.apply (lambda row: addmmpiKaon(row), axis=1) 

        dataset_all['wt_kstarmass']= dataset_all.apply (lambda row: addkst2(row), axis=1) 
        
        dataset_all['kaonPt']   = dataset_all.apply (lambda row: addpi1Pt(row), axis=1)
        dataset_all['pionPt']   = dataset_all.apply (lambda row: addpi2Pt(row), axis=1)
        

        dataset_all['mmpiMass'], dataset_all['mmkMass'] = addmmpi2(
                                        dataset_all.mumPt,  dataset_all.mumEta,  dataset_all.mumPhi,  
                                        dataset_all.mupPt,  dataset_all.mupEta,  dataset_all.mupPhi,  
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                        dataset_all.tagB0
                                      )
                                      
        dataset_all['mmkkMass']= addmmkkmass(
                                        dataset_all.mumPt,  dataset_all.mumEta,  dataset_all.mumPhi,  
                                        dataset_all.mupPt,  dataset_all.mupEta,  dataset_all.mupPhi,  
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                      )
        dataset_all['mmpipiMass']= addmmpipimass(
                                        dataset_all.mumPt,  dataset_all.mumEta,  dataset_all.mumPhi,  
                                        dataset_all.mupPt,  dataset_all.mupEta,  dataset_all.mupPhi,  
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                      )
        dataset_all['pipiMass']= addpipiMass(
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                      )
        dataset_all['xcut'] = addxcut(dataset_all.wt_mass, dataset_all.wt_kstarmass, dataset_all.kaonPt, dataset_all.pionPt, dataset_all.mmpiMass, dataset_all.mmkMass )
        dataset_all['xcut'] = dataset_all['xcut'].astype(np.int32)        
        dataset_all['muptrkm_pik'] = addpiKMass(
                                        dataset_all.mupPt,  dataset_all.mupEta,  dataset_all.mupPhi,  
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                      )
        dataset_all['muptrkm_kp'] = addpiKMass(
                                        dataset_all.kstTrkmPt,  dataset_all.kstTrkmEta,  dataset_all.kstTrkmPhi,
                                        dataset_all.mupPt,  dataset_all.mupEta,  dataset_all.mupPhi,  
                                      )
        dataset_all['mumtrkp_pik'] = addpiKMass(
                                        dataset_all.mumPt,  dataset_all.mumEta,  dataset_all.mumPhi,  
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                      )
        dataset_all['mumtrkp_kpi'] = addpiKMass(
                                        dataset_all.kstTrkpPt,  dataset_all.kstTrkpEta,  dataset_all.kstTrkpPhi,
                                        dataset_all.mumPt,  dataset_all.mumEta,  dataset_all.mumPhi,  
                                      )
         

    
        import root_pandas
        dataset_all.to_root(ofile, key='ntuple')#, store_index=False)
    