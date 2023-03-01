import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("dimusel"    , help = "Define if keep or remove dimuon resonances", choices=['rejectPsi', 'keepJpsi', 'keepPsiP'])
parser.add_argument("year"       , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--mcw"      , help = "", action='store_true')
parser.add_argument("--scaleErr" , help = "", action='store_true')
args = parser.parse_args()

'''
code to fit the B0 mass distribution:
- unbinned fit
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
'''

import os, sys, inspect
from os import path
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import ROOT
from ROOT import gSystem
ROOT.gROOT.SetBatch(True)

gSystem.Load('libRooFit')
gSystem.Load('../utils/func_roofit/libRooDoubleCBFast')
gSystem.Load('../utils/func_roofit/libRooGaussDoubleSidedExp')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial, RooExtendPdf
import sys, math, pdb
from uncertainties import ufloat
import random
import numpy as np

ROOT.RooMsgService.instance().setGlobalKillBelow(4)
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(50000)


def _getFittedVar(varName, w=None):
    if w is not None:
        return ufloat (w.var(varName).getVal() , w.var(varName).getError())
    else :
        return ufloat (varName.getVal()        , varName.getError())

def _goodFit(r):
    return (r.status()==0 and r.covQual() == 3)

def _accFit(r):
    return (r.status()==4 and r.covQual() == 3)

def _writeFitStatus(r):
    str_status = "GOOD" if r.status()==0 else "NOT CONV"
    txt = ROOT.TLatex(.16,.7, "fit status: " + str_status + ", covQ = %s" %r.covQual() )
    txt . SetNDC() ;
    txt . SetTextSize(0.033) ;
    txt . SetTextFont(42)
    return txt

def _writeChi2(chi2):
    txt = ROOT.TLatex(.16,.8, "fit #chi^{2}/ndf: %.1f "%chi2 )
    txt . SetNDC() ;
    txt . SetTextSize(0.033) ;
    txt . SetTextFont(42)
    return txt
    
def _constrainVar(var, nsigma):
    
    constr = _getFittedVar(var.GetName(), w)
    gauss_constr = RooGaussian(  "c_%s" %var.GetName() , 
                                 "c_%s" %var.GetName() , 
                                var         ,  
                                ROOT.RooFit.RooConst( constr.n ), 
                                ROOT.RooFit.RooConst( nsigma*constr.s )
                                ) 
    print 'constraining var',   var.GetName(), ': ',     constr.n , ' with uncertainty ' , nsigma*constr.s                          
    return gauss_constr                        


sys.path.append("../utils")
from utils.utils import *
from utils.fit_functions import *

nbins = 80
nSigma_psiRej = 3.
cut_base      = 'passB0Psi_lmnr == 1 '
## entries per bin in data
n_bin = n_data[args.year]
n_bin = 7.1112E5
## uncertainty on fRT
fm_sigma = fM_sigmas[args.year]

q2binning = [
                1,
                2, 
                4.3,
                6,
                8.68,
                10.09,
                12.86,
                14.18,
                16,
#                 19,
]

from collections import OrderedDict
fitStats = OrderedDict()
covStats = OrderedDict()
chi2s    = OrderedDict()

  
def fitMC(fulldata, correctTag, ibin):

    print 'now fitting: ', ibin, ' for ', correctTag*'correctTag ', (1-correctTag)*'misTag'  
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    if ibin == 4:
        cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s) && xcut==0'%(q2binning[ibin], q2binning[ibin+1])
    fulldata_v2 = fulldata.reduce(RooArgSet(thevarsMC), cut)
    data = fulldata_v2.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)


    pol_c1      = RooRealVar ("p1"           , "coeff x^0 term" ,  -0.5,   -10, 10);
    pol_c2      = RooRealVar ("p2"           , "coeff x^0 term" ,  0.5,   -10, 10);
    bkg_pol     = RooChebychev("bkg_pol"     , "bkg_pol" ,  tagged_mass, RooArgList(pol_c1));
    slope       = RooRealVar    ("slope"   , "slope"           ,   -4,   -10, 10);
    bkg_exp     = RooExponential("bkg_exp" , "exponential"     ,  slope,   tagged_mass  );
    fsig        = RooRealVar("fsig"         , "fsig"   ,   0.9, 0, 1)

    nsig        = RooRealVar("Yield"         , "nsig"   ,   500000,     0,    1E7)
    nbkg        = RooRealVar("nbkg"          , "nbkg"   ,     1000,     0,    1E6 )
    
    doextended = False
    fitrange   = "mcrange"
    tag        = 'RT' if correctTag else 'WT'
    if correctTag:
        fitrange   = "datarange"
        mean        = RooRealVar ("mean_{RT}^{%s}"%ibin,        "massRT"         , B0Mass_,     5,    6, "GeV")
        if ibin < 4:  #4  
            ### double CB, 1 sigma
            sigmaCB     = RooRealVar ("#sigma_{RT1}^{%s}"%ibin, "sigmaCB"        ,  0.03  ,     0,   1, "GeV"  )#1
            alpha1      = RooRealVar ("#alpha_{RT1}^{%s}"%ibin,  "#alpha_{1}"     ,  0.5   ,    0,  10  )#10
            alpha2      = RooRealVar ("#alpha_{RT2}^{%s}"%ibin,  "#alpha_{2}"     ,  2.5   ,    0,  10  )#10
            n1          = RooRealVar ("n_{RT1}^{%s}"%ibin,       "n_1"            ,  1     ,    0.01,  20  )
            n2          = RooRealVar ("n_{RT2}^{%s}"%ibin,       "n_2"            ,  1     ,    0.01,  100  )
            doublecb = ROOT.RooDoubleCBFast("doublecb_RT%s"%ibin, "doublecb", tagged_mass, mean, sigmaCB, alpha1, n1, alpha2, n2)	
            signalFunction = doublecb

        elif ibin == 7:  #4  
            ### CB + gaus
            crystalBall( mean  , initial_sigma1,  1.5,  1, tagged_mass, w, 'RT1', ibin, [0, 10])
            singleG(mean, initial_sigma2, tagged_mass, w, 'RT2' , ibin)
            gausCB( w.pdf("cbshape_RT1_%s"%ibin) , w.pdf("gaus_RT2_%s"%ibin), 0.5, tagged_mass, w, 'RT', ibin )
            signalFunction = w.pdf("gauscb_RT_%s"%ibin)   

        else:
            ### double CB, 2 sigmas
            crystalBall( mean  , initial_sigma1,  1.5,  1, tagged_mass, w, 'RT1', ibin, [0, 10])
            crystalBall( mean  , initial_sigma2,  -2 ,  1, tagged_mass, w, 'RT2', ibin, [-10, 0])
            doubleCB     ( w.pdf("cbshape_RT1_%s"%(ibin)), w.pdf("cbshape_RT2_%s"%ibin), 0.8  , tagged_mass, w, "RT%s"%ibin)
            signalFunction = w.pdf("doublecb_RT%s"%ibin)   
          
        fitFunction = signalFunction  ## sara 04.11 to uncomment

    else:
        ### double CB, 1 sigma
        mean        = RooRealVar ("mean_{WT}^{%s}"%ibin,     "massWT"         , B0Mass_,    5,   6, "GeV")
        sigmaCB     = RooRealVar ("#sigma_{WT1}^{%s}"%ibin,  "sigmaCB"        ,  0.03  ,    0,   1, "GeV"  )
        alpha1      = RooRealVar ("#alpha_{WT1}^{%s}"%ibin,  "#alpha_{1}"     ,  0.5   ,    0,  10  )
        alpha2      = RooRealVar ("#alpha_{WT2}^{%s}"%ibin,  "#alpha_{2}"     ,  2.5   ,    0,  10  )
        n1          = RooRealVar ("n_{WT1}^{%s}"%ibin,       "n_1"            ,  1     ,    0.01,  100  )
        n2          = RooRealVar ("n_{WT2}^{%s}"%ibin,       "n_2"            ,  1     ,    0.01,  100  )
        doublecb = ROOT.RooDoubleCBFast("doublecb_%s"%ibin, "doublecb", tagged_mass, mean, sigmaCB, alpha1, n1, alpha2, n2)	
        signalFunction = doublecb
        fitFunction    = doublecb  

#     fitFunction = RooAddPdf("fitfunction%s"%ibin, "fitfunction%s"%ibin , RooArgList(signalFunction, bkg_exp), RooArgList(fsig))

    getattr(w,"import")(signalFunction)
    
    r = fitFunction.fitTo(data, RooFit.Extended(doextended), RooFit.Save(), RooFit.Range(fitrange))
    print 'fit status: ', r.status(), r.covQual() 
    total_integral = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('fitrange')).getVal()
    center = mean.getVal()
    initial_sigma = 0.032
    for dsigma in np.linspace(0, 0.02, 200):
        tmp_sigma = initial_sigma + dsigma
        tagged_mass.setRange('tmp_range', center-tmp_sigma, center+tmp_sigma)
        tmp_int = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('tmp_range')).getVal()
        if tmp_int/total_integral > 0.68:
          print tag, ' reached sigma:', tmp_sigma, ' for which we have ', tmp_int/total_integral 
          break
    
    ## draw everything 
    params = signalFunction.getParameters(RooArgSet(tagged_mass)) 
    w.saveSnapshot("reference_fit_%s_%s"%(tag, ibin),params,ROOT.kTRUE) 
    frame = tagged_mass.frame(RooFit.Range(fitrange))#, RooFit.Title('correctly'*correctTag + 'wrongly'*(1-correctTag) + ' tagged events'))
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    
    drawPdfComponents(fitFunction, frame, ROOT.kGreen if correctTag else ROOT.kViolet, RooFit.NormRange(fitrange), RooFit.Range(fitrange), isData=False)
    fitFunction.plotOn(frame, RooFit.NormRange(fitrange), RooFit.Range(fitrange) )
    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))
    
    frame.Draw()
    niceFrame(frame, 'correctly '*correctTag + 'mis'*(1-correctTag) + 'tagged events')
    frame. addObject(_writeFitStatus(r))
    
    if ibin >= 4 and correctTag:
      s2 = w.var("#sigma_{RT2}^{%s}"%(ibin))
      tmp_sigma_err = tmp_sigma * s2.getError() / s2.getVal()
      txt_width = ROOT.TLatex(.62,.3, "#sigma_{eff} = %.4f #pm %.4f GeV"%(tmp_sigma, tmp_sigma_err) )
      txt_width . SetNDC() ;
      txt_width . SetTextSize(0.026) ;
      txt_width . SetTextFont(42)
      frame. addObject(txt_width)
    
    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams    = fitFunction.getParameters(observables)
    nparam      = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    pdfstring = ''

    if correctTag:
#         pdfstring = "doublegaus_RT%s_Norm[tagged_mass]_Comp[doublegaus_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)
        pdfstring = "doublecb_RT%s_Norm[tagged_mass]_Comp[doublecb_RT%s]_Range[datarange]_NormRange[datarange]"%(ibin,ibin)
#         pdfstring = "gauscb_RT%s_Norm[tagged_mass]_Comp[gauscb_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)
#         pdfstring = "expGaussExp_RT%s_Norm[tagged_mass]_Comp[expGaussExp_RT%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)

        if ibin == 7:   
            pdfstring = "gauscb_RT_%s_Norm[tagged_mass]_Comp[gauscb_RT_%s]_Range[datarange]_NormRange[datarange]"%(ibin,ibin)
        if doextended:
            dict_s_rt[ibin]   = _getFittedVar(nsig)
        else:
            dict_s_rt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
        nRT = RooRealVar ("nRT_%s"%ibin, "yield of RT signal",0,1.E9)
        nRT.setVal(  dict_s_rt[ibin].n)
        nRT.setError(dict_s_rt[ibin].s)
        print 'setting nRT to ', dict_s_rt[ibin].n
        getattr(w,"import")(nRT)

    else:
        pdfstring = "doublecb_%s_Norm[tagged_mass]_Comp[doublecb_%s]_Range[mcrange]_NormRange[mcrange]"%(ibin,ibin)
        
        dict_s_wt[ibin]    = ufloat(data.sumEntries(), math.sqrt(data.sumEntries()))
        nWT = RooRealVar ("nWT_%s"%ibin, "yield of WT signal",0,1.E7)
        nWT.setVal(  dict_s_wt[ibin].n)
        nWT.setError(dict_s_wt[ibin].s)
        print 'setting nWT to ', dict_s_wt[ibin].n
        getattr(w,"import")(nWT)
    
    ## eval and save goodness of fit indicators
    chi2s['%s%s'%(tag,ibin)] = frame.chiSquare(pdfstring, "h_fullmc",  nparam)
    frame. addObject(_writeChi2( chi2s['%s%s'%(tag,ibin)] ))
    fitStats['%s%s'%(tag,ibin)] = r.status()
    covStats['%s%s'%(tag,ibin)] = r.covQual()

    c1 = ROOT.TCanvas() 
    upperPad = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)

    upperPad.Draw()
    lowerPad.Draw()

    upperPad.cd()
    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], -1)
    frame.Draw()
    
    ## add plot of pulls
    lowerPad.cd()
    hpull  = frame.pullHist("h_fullmc", pdfstring)
    frame2 = tagged_mass.frame(RooFit.Range(fitrange), RooFit.Title(''))
    frame2.addPlotable(hpull,"P") 
    niceFrameLowerPad(frame2, 'pull')
    frame2.Draw()
    line = ROOT.TLine(4.9,0,5.7,0)
    if correctTag:
        line = ROOT.TLine(5.0,0,5.6,0)
    line.SetLineColor(ROOT.kGreen+3)
    line.Draw()
    
    ## save to pdf and root files
    for ilog in [True,False]:
        upperPad.SetLogy(ilog)
        c1.SaveAs('MC_fit_results_mass/noIP2D/save_fit_mc_{IBIN}_{YEAR}_noIP2D_{TAG}{LOG}{CHAN}{MCW}{SCALE}.pdf'.format(IBIN = ibin, 
                                                                                                                        YEAR = args.year, 
                                                                                                                        TAG = tag, 
                                                                                                                        LOG = '_logScale'*ilog, 
                                                                                                                        CHAN = '_Jpsi'*(args.dimusel=='keepJpsi'),
                                                                                                                        MCW = '_MCw_xgbv8'*(args.mcw==True), 
                                                                                                                        SCALE = '_scaleErr'*(args.scaleErr==True)))
    out_f.cd()
    r.Write('results_%s_%s'%(tag, ibin))
    
    return dict_s_rt[ibin].n if correctTag else dict_s_wt[ibin].n
   
   




# tData = ROOT.TChain('ntuple')
tMC = ROOT.TChain('ntuple')
tFriend = ROOT.TChain('wTree')

string_nonan = '' + ('_noNan' * (args.year == '2017'))

if args.dimusel == 'rejectPsi':
    tMC.Add('/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/%sMC_LMNR_noIP2D%s_addxcutvariable.root'%(args.year,string_nonan))
    if args.mcw == True:
        tFriend.Add('/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-LMNR-XGBv8/%s.root'%args.year)
elif args.dimusel == 'keepJpsi':
    tMC.Add('/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/%sMC_JPSI_noIP2D%s_addxcutvariable.root'%(args.year,string_nonan))
    if args.mcw == True:
        tFriend.Add('/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-Jpsi-XGBv8/%s.root'%args.year)
elif args.dimusel == 'keepPsiP':
    tMC.Add('/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/%sMC_PSI_noIP2D%s_addxcutvariable.root'%(args.year,string_nonan))
    if args.mcw == True:
        tFriend.Add('/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-Psi-XGBv8/%s.root'%args.year)

if args.year == 'test':
    tMC = ROOT.TChain('ntuple')
    tMC.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016MC_LMNR_100k.root')

print 'mc file name:',  tMC.GetFile().GetName()

tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.7, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);
passB0Psi_lmnr  = RooRealVar("passB0Psi_lmnr" , "passB0Psi_lmnr", -200, 2);
passB0Psi_jpsi  = RooRealVar("passB0Psi_jpsi" , "passB0Psi_jpsi", -200, 2);
passB0Psi_psip  = RooRealVar("passB0Psi_psip" , "passB0Psi_psip", -200, 2);
weight          = RooRealVar("weight" , "weight", 0,10);
sf_to_data      = RooRealVar("sf_to_data" , "sf_to_data", 0,999999);
xcut            = RooRealVar("xcut" , "xcut", -1, 5);
MCw             = RooRealVar("MCw" , "MCw", 0,999999);

tagged_mass.setRange("datarange", 5.0,5.6) ;
# tagged_mass.setRange("datarange"  , 4.9,5.7) ;
tagged_mass.setRange("mcrange"  , 4.9,5.7) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)
thevars.add(passB0Psi_lmnr)
thevars.add(passB0Psi_jpsi)
thevars.add(passB0Psi_psip)
thevars.add(xcut)

genSignal       = RooRealVar("genSignal"      , "genSignal"      , 0, 10);
thevarsMC   = thevars; 
thevarsMC.add(genSignal)
thevarsMC.add(weight)
fullmc = RooDataSet()
if args.mcw == False:
  fullmc      = RooDataSet('fullmc', 'fullmc', tMC,  RooArgSet(thevarsMC), "", "weight")
else:
  thevarsMC.add(MCw)
  fullmc   = RooDataSet('fullmc', 'fullmc', tMC,  RooArgSet(thevarsMC), "", "weight*MCw")


# c10 = ROOT.TCanvas() 
# MCw.setRange("mcwrange"  , -2,4) ;
# wframe = MCw.frame(RooFit.Range("mcwrange"))
# fullmc.plotOn(wframe)
# wframe.Draw()
# c10.SaveAs('weights.pdf')

print 'RooDataSet created'
nMCEntries = fullmc.numEntries()
print '->  ', nMCEntries, ' Num entries'
print '->  ', fullmc.sumEntries(), ' Sum entries'
randVar  = RooRealVar("rand","rand",0,1) 
p0       = RooPolynomial("px","px",randVar) ;
rDataset = p0.generate(RooArgSet(randVar),int(nMCEntries))
fullmc.merge(rDataset) 
thevarsMC.add(randVar)

# nDesired = n_bin/nMCEntries
# nev_cut = 'rand < %f'%nDesired
# part_mc = fullmc.reduce(thevarsMC, nev_cut)


### define correct and wrong tag samples
rt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2))')
wt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '((tagB0==0 && genSignal==1) || (tagB0==1 && genSignal==2))')

dict_s_rt  = {}
dict_s_wt  = {}

out_fname = "MC_fit_results_mass/noIP2D/results_fits_%s_fM%s_noIP2D%s%s.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP'), '_MCw_xgbv8'*(args.mcw==True), '_scaleErr'*(args.scaleErr==True))
out_f = TFile (out_fname, "RECREATE") 
w = ROOT.RooWorkspace("w")
initial_n_1 =  3.
initial_n_2 =  1.
initial_a_1 =  1.
initial_a_2 = -1.
initial_sigma1 = 0.028
initial_sigma2 = 0.048
initial_sigmaCB = 0.048
nRT_fromMC = 10
nWT_fromMC = 10


for ibin in range(len(q2binning)-1):

    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
    if args.dimusel == 'keepJpsi' and \
       (q2binning[ibin] < 8.68 or q2binning[ibin] > 10): 
           continue
    if args.dimusel == 'keepPsiP' and \
       (q2binning[ibin] < 12 or q2binning[ibin] > 14): 
           continue
    if   ibin ==4:  cut_base = 'passB0Psi_jpsi== 1 '
    elif ibin ==6:  cut_base = 'passB0Psi_psip== 1 '

    nRT_fromMC = fitMC(rt_mc, True, ibin)
    nWT_fromMC = fitMC(wt_mc, False, ibin)
    print ' --------------------------------------------------------------------------------------------------- '


print '--------------------------------------------------------------------------------------------------- '
print 'bin\t\t fit status \t cov. matrix \t\t chi2'
# for i,k in enumerate(fitStats.keys()):    
#     if i%3==0:  print '------------------------------------------------------'
#     print k , '\t\t', fitStats[k], '\t\t', covStats[k], '\t\t', chi2s[k]
# 
for i,k in enumerate(dict_s_rt.keys()):    
    print 'bin ' , i , '\t\t n rt events: ', dict_s_rt[k], '\t\t n wt events: ', dict_s_wt[k], '\t\t n tot events: ', dict_s_wt[k]+dict_s_rt[k]

out_f.Close()
w.writeToFile(out_f.GetName(), False)
