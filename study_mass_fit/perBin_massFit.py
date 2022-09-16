import argparse
parser = argparse.ArgumentParser(description="")
# parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi", choices=['rejectPsi', 'keepJpsi', 'keepPsiP'])
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()


'''
code to fit the B0 mass distribution:
- unbinned fit
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
'''

import os, sys, inspect
import numpy as np
from os import path
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

from copy import deepcopy as dc
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import ROOT
from ROOT import gSystem
ROOT.gROOT.SetBatch(True)

gSystem.Load('libRooFit')
gSystem.Load('../utils/func_roofit/libRooDoubleCBFast')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial, RooExtendPdf
import sys, math, pdb
from uncertainties import ufloat

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

nbins = 60
nSigma_psiRej = 3.
cut_base      = 'passB0Psi_lmnr == 1 '

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



   
def fitData(fulldata, ibin, w):

    cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s) && (tagged_mass > 5. && tagged_mass < 5.6) '%(q2binning[ibin], q2binning[ibin+1])
    if ibin == 4:
        cut  = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s) && (tagged_mass > 5. && tagged_mass < 5.6) && xcut==0'%(q2binning[ibin], q2binning[ibin+1])
    data = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)
    data.Print("V")
    tagged_mass.setMin(5.0) ;
    tagged_mass.setMax(5.6) ;
    
    nrt_mc = _getFittedVar("nRT_%s"%ibin, w)
    nwt_mc = _getFittedVar("nWT_%s"%ibin, w)
    fraction = nwt_mc / (nrt_mc + nwt_mc)
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , '+/-', fraction.s 
    
    ### creating WT component
    w.loadSnapshot("reference_fit_WT_%s"%ibin)
    mean_wt      = w.var("mean_{WT}^{%s}"%ibin)
    sigma_wt     = w.var("#sigma_{WT1}^{%s}"%ibin)
    alpha_wt1    = w.var("#alpha_{WT1}^{%s}"%ibin)
    alpha_wt2    = w.var("#alpha_{WT2}^{%s}"%ibin)
    n_wt1        = w.var("n_{WT1}^{%s}"%ibin)
    n_wt2        = w.var("n_{WT2}^{%s}"%ibin)
    theWTgauss   = w.pdf("doublecb_%s"%ibin)   

    n_wt1.setMin(0.01)  
    n_wt2.setMin(0.01)  

    c_mean_wt     = _constrainVar(mean_wt,  1)
    c_sigma_wt    = _constrainVar(sigma_wt,  1)
    c_alpha_wt1   = _constrainVar(alpha_wt1, 1)
    c_alpha_wt2   = _constrainVar(alpha_wt2, 1)
    c_n_wt1       = _constrainVar(n_wt1, 1)
    c_n_wt2       = _constrainVar(n_wt2, 1)

    c_pdfs_wt = RooArgSet(c_mean_wt, c_sigma_wt, c_alpha_wt1, c_alpha_wt2, c_n_wt1, c_n_wt2)
    c_vars_wt = RooArgSet(  mean_wt,   sigma_wt,   alpha_wt1,   alpha_wt2,   n_wt1,   n_wt2)

    constr_list_wt = RooArgList(c_pdfs_wt)
    constr_list_wt.add(theWTgauss)
    c_theWTgauss = RooProdPdf ("c_theWTgauss%s"%ibin, "c_theWTgauss", constr_list_wt)     

    c_vars = RooArgSet()
    c_vars.add(c_vars_wt);


    ### creating RT component
    w.loadSnapshot("reference_fit_RT_%s"%ibin)
    mean_rt      = w.var("mean_{RT}^{%s}"%ibin)

    sigma_rt1    = w.var("#sigma_{RT1}^{%s}"%ibin)
    c_sigma_rt1  = _constrainVar(sigma_rt1, 1)

    alpha_rt1    = w.var("#alpha_{RT1}^{%s}"%ibin)
    c_alpha_rt1  = _constrainVar(alpha_rt1, 1)

    n_rt1        = w.var("n_{RT1}^{%s}"%ibin)
    c_n_rt1      = _constrainVar(n_rt1, 1)
    n_rt1.setMin(0.01)  

    sigma_rt2 = RooRealVar()
    alpha_rt2 = RooRealVar()
    n_rt2     = RooRealVar()
    f1rt     = RooRealVar()

    ### creating variable for the difference between the two peaks
    # deltaPeaks = RooFormulaVar("deltaPeaks%s"%ibin, "@0 - @1", RooArgList(mean_rt, mean_wt))  
    mRT_ufloat = ufloat(mean_rt.getVal(), mean_rt.getError())
    mWT_ufloat = ufloat(mean_wt.getVal(), mean_wt.getError())
    deltaPeakValue = mRT_ufloat - mWT_ufloat 
    deltaPeakVar = RooRealVar("deltaPeakVar%s"%ibin, "deltaPeakVar%s"%ibin, deltaPeakValue.n, 0, 0.2)  

    mRT_data = RooFormulaVar("mRT_data%s"%ibin, "@0 + @1", RooArgList(mean_wt , deltaPeakVar))  
    ### creating constraints for the difference between the two peaks (not used)
    c_deltaPeaks = RooGaussian("c_deltaPeaks%s"%ibin , "c_deltaPeaks", deltaPeakVar, ROOT.RooFit.RooConst( deltaPeakValue.n ), 
                                ROOT.RooFit.RooConst( deltaPeakValue.s )  ## value to be checked
                                ) 


    if ibin in range(0,4):
        alpha_rt2    = w.var("#alpha_{RT2}^{%s}"%ibin)
        c_alpha_rt2  = _constrainVar(alpha_rt2, 1)
        n_rt2        = w.var("n_{RT2}^{%s}"%ibin)
        c_n_rt2      = _constrainVar(n_rt2, 1)
        n_rt2.setMin(0.01)  
        c_vars_rt = RooArgSet(  sigma_rt1,   alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2)
        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2)
        theRTgauss = ROOT.RooDoubleCBFast("doublecb_RT%s"%ibin, "doublecb", tagged_mass, mRT_data, sigma_rt1, alpha_rt1, n_rt1, alpha_rt2, n_rt2)	
        
    elif ibin in range(4,7):
        sigma_rt2    = w.var("#sigma_{RT2}^{%s}"%ibin)
        c_sigma_rt2  = _constrainVar(sigma_rt2, 1)
        alpha_rt2    = w.var("#alpha_{RT2}^{%s}"%ibin)
        c_alpha_rt2  = _constrainVar(alpha_rt2, 1)
        n_rt2        = w.var("n_{RT2}^{%s}"%ibin)
        c_n_rt2      = _constrainVar(n_rt2, 1)
        n_rt2.setMin(0.01)  
        f1rt         = w.var("f^{RT%s}"%ibin)
        c_f1rt       = _constrainVar(f1rt, 1)
        c_vars_rt = RooArgSet(  sigma_rt1,   sigma_rt2,   alpha_rt1,   alpha_rt2,   n_rt1,   n_rt2,   f1rt)
        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2, c_f1rt)

        cbshape1      = RooCBShape ("cbshape_RT1_%s"%(ibin), "cbshape_RT1_%s"%(ibin),  tagged_mass, mRT_data, sigma_rt1, alpha_rt1, n_rt1)
        cbshape2      = RooCBShape ("cbshape_RT2_%s"%(ibin), "cbshape_RT2_%s"%(ibin),  tagged_mass, mRT_data, sigma_rt2, alpha_rt2, n_rt2)
        theRTgauss =  RooAddPdf ("doublecb_RT%s"%ibin, "doublecb"  ,  RooArgList(cbshape1,cbshape2), RooArgList(f1rt))
  

    elif ibin == 7:
        sigma_rt2    = w.var("#sigma_{RT2}^{%s}"%ibin)
        c_sigma_rt2  = _constrainVar(sigma_rt2, 1)
        f1rt         = w.var("f^{RT%s}"%ibin)
        c_f1rt       = _constrainVar(f1rt, 1)
        c_pdfs_rt = RooArgSet(c_sigma_rt1, c_sigma_rt2, c_alpha_rt1, c_n_rt1, c_f1rt)
        c_vars_rt = RooArgSet(  sigma_rt1,   sigma_rt2,   alpha_rt1,   n_rt1,   f1rt)

        cbshape1      = RooCBShape ("cbshape_RT1_%s"%(ibin), "cbshape_RT1_%s"%(ibin),  tagged_mass, mRT_data, sigma_rt1, alpha_rt1, n_rt1)
#         singleG(mRT_data, sigma_rt2, tagged_mass, out_w, 'RT2' , ibin)
        singlegaus   = RooGaussian("gaus_RT2_%s"%bin, "singlegaus", tagged_mass,  mRT_data, sigma_rt2)
        theRTgauss =  RooAddPdf ("doublecb_RT%s"%ibin, "doublecb"  ,  RooArgList(cbshape1,singlegaus), RooArgList(f1rt))
#         gausCB( out_w.pdf("cbshape_RT1_%s"%ibin) , out_w.pdf("gaus_RT2_%s"%ibin), f1rt, tagged_mass, out_w, 'RT', ibin )
#         theRTgauss = out_w.pdf("gauscb_RT_%s"%ibin)   
 
  
#     the_string = 'double' if ibin !=7 else 'gaus'
#     theRTgauss  = w.pdf("%scb_RT%s%s"%(the_string,'_'*(ibin==7), ibin))   

#     c_vars_rt.add(mRT_data) 
#     c_pdfs_rt.add(c_deltaPeaks) 
    constr_list_rt = RooArgList(c_pdfs_rt)
    constr_list_rt.add(theRTgauss)
    c_theRTgauss = RooProdPdf ("c_theRTgauss%s"%ibin, "c_theRTgauss", constr_list_rt)    
    c_vars.add(c_vars_rt);

    fm              = RooRealVar ("f_{M}^{%s}"%ibin , "fm"             , fraction.n , 0, 1)
    signalFunction  = RooAddPdf  ("sumgaus%s"%ibin  , "rt+wt"           , RooArgList(c_theWTgauss,c_theRTgauss), RooArgList(fm))

    
    ### creating constraints on mistag fraction
    c_fm = RooGaussian("c_fm%s"%ibin, "c_fm" , fm,  ROOT.RooFit.RooConst(fraction.n) , ROOT.RooFit.RooConst(fm_sigma[ibin]) )
    c_pdfs = RooArgSet(c_fm)
    c_vars.add(fm)

    c_pdf_list = RooArgList(signalFunction)
    c_pdf_list.add(c_pdfs)
    c_signalFunction = RooProdPdf ("c_signalFunction%s"%ibin, "c_signalFunction", c_pdf_list)     


    ### now create background parametrization
    slope         = RooRealVar    ("slope"      , "slope"           ,    -6,   -10, 10);
#     slope         = RooRealVar    ("slope"      , "slope"           ,    0.5,   -10, 10);
    bkg_exp       = RooExponential("bkg_exp"    , "exponential"     ,  slope,   tagged_mass  );
    pol_c1        = RooRealVar    ("p1"         , "coeff x^0 term"  ,    0.5,   -10, 10);
    pol_c2        = RooRealVar    ("p2"         , "coeff x^1 term"  ,    0.5,   -10, 10);
    bkg_pol       = RooChebychev  ("bkg_pol"    , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1));
   
    nsig          = RooRealVar("Yield"         , "signal yield"     ,     1000,     0,   1000000);
    nbkg          = RooRealVar("nbkg"          , "bkg fraction"     ,     1000,     0,   5500000);
    if ibin == 4 :
      if args.year =='2018':
        nsig          = RooRealVar("Yield"         , "signal frac"    ,     1500000,     0,   3E6);
        nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,      500000,     0,   2E6);
      if args.year =='2017' or args.year =='2016':
        nsig          = RooRealVar("Yield"         , "signal frac"    ,      719000,     500000,   1E6);
        nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,      164000,     120000,   3E5);
#         nsig          = RooRealVar("Yield"         , "signal frac"    ,      700000,     500000,   1E6);
#         nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,      170000,     100000,   5E5);
    elif ibin == 6:
        nsig          = RooRealVar("Yield"         , "signal frac"    ,     100000,     0,   1E6);
        nbkg          = RooRealVar("nbkg"          , "bkg fraction"   ,      60000,     0,   1E6);


    fitFunction  = RooAddPdf  ("fitfunction%s"%ibin   , "fitfunction%s"%ibin , RooArgList(c_signalFunction,bkg_exp), RooArgList(nsig,nbkg))
#     fitFunction.fitTo(data, 
#                           RooFit.Extended(True), 
# #                           RooFit.Save(), 
#                           RooFit.Range("full"), 
#                           RooFit.Verbose(False),
#                           ROOT.RooFit.Constrain(c_vars),
#                           ROOT.RooFit.NumCPU(64)
#                          )
#     print 'finished first tmp fit ------------------------------------------------------------------------------'                     
    r = fitFunction.fitTo(data, 
                          RooFit.Extended(True), 
                          RooFit.Save(), 
                          RooFit.Range("full"), 
                          RooFit.Verbose(False),
                          ROOT.RooFit.Constrain(c_vars),
                          ROOT.RooFit.NumCPU(64)
                         )
    print 'finished FINAL tmp fit ------------------------------------------------------------------------------'                     
    print 'from fit: ', nsig.getVal() + nbkg.getVal()
    print 'data:', data.numEntries()

    r.Print()
    print("c_vars:  ")
    c_vars.Print()

#     r.correlationMatrix().Print()
    fitStats['data%s'%(ibin)] = r.status()
    covStats['data%s'%(ibin)] = r.covQual()
    frame = tagged_mass.frame( RooFit.Range("full") )
    data.plotOn(frame, RooFit.Binning(nbins), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, RooFit.NormRange("full"), RooFit.Range("full"))

    ## evaluate sort of chi2 and save number of RT/WT events
    observables = RooArgSet(tagged_mass)
    flparams = fitFunction.getParameters(observables)
    nparam = int(flparams.selectByAttrib("Constant",ROOT.kFALSE).getSize())
    pdfstring = "fitfunction%s_Norm[tagged_mass]_Range[full]_NormRange[full]"%ibin
    chi2s['data%s'%ibin] = frame.chiSquare(pdfstring, "h_fulldata",  nparam)
    frame. addObject(_writeChi2( chi2s['data%s'%ibin] ))

    drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("full"), RooFit.Range("full"), isData = True)
#     fitFunction.paramOn(frame, RooFit.Layout(0.62,0.86,0.89))

#     parList = RooArgSet (nsig, mean_rt, mean_wt, sigma_rt1, alpha_rt1, n_rt1, n_rt2, alpha_rt2, mean_wt)
    parList = RooArgSet (nsig, nbkg, mean_rt, mean_wt, sigma_rt1, mean_wt)
    parList.add(sigma_wt)

    if ibin >= 4:
        parList.add(sigma_rt2)
    parList.add(fm)
#     fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.56,0.86,0.89))

    frame.Draw()
    niceFrame(frame, '')
    frame. addObject(_writeFitStatus(r))

#     pdb.set_trace()
    ## get the signal component
#     center = mean_rt.getVal()
#     initial_sigma = 0.02
# 
# #     set_tm = ROOT.RooArgSet(tagged_mass)
# #     normSet = ROOT.RooFit.NormSet(set_tm)
# #     inttot = theWTgauss.createIntegral(set_tm, normSet, ROOT.RooFit.Range('tmp_range1'))
# #     int2   = theWTgauss.createIntegral(set_tm, normSet, ROOT.RooFit.Range('tmp_range2'))
# 
#     tot_int = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.Range('full')).getVal()
#     for dsigma in np.linspace(0,0.1,100):
#         tmp_sigma = initial_sigma + dsigma
#         tagged_mass.setRange('tmp_range', center-tmp_sigma, center+tmp_sigma)
#         tmp_int = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.Range('tmp_range')).getVal()
#         tmp_int2 = signalFunction.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('tmp_range')).getVal()
# #         pdb.set_trace()
# #         print tmp_sigma, tmp_int, tmp_int2, tot_int
#         if tmp_int/tot_int > 0.68:
#           print ' overall signal sigma:', tmp_sigma, ' for which we have ', tmp_int/tot_int 
#           break
# 
# #     pdb.set_trace()
#     tot_wt = theWTgauss.createIntegral(RooArgSet(tagged_mass), RooFit.Range('full')).getVal()
#     for dsigma in np.linspace(0,0.1,100):
#         tmp_sigma = initial_sigma + dsigma
#         tagged_mass.setRange('tmp_range', center-tmp_sigma, center+tmp_sigma)
#         tmp_int = theWTgauss.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('tmp_range')).getVal()
#         if tmp_int/tot_wt > 0.68:
#           print ' theWTgauss signal sigma:', tmp_sigma, ' for which we have ', tmp_int/tot_wt 
#           break
# 
#     tot_rt = theRTgauss.createIntegral(RooArgSet(tagged_mass), RooFit.Range('full')).getVal()
#     for dsigma in np.linspace(0,0.1,100):
#         tmp_sigma = initial_sigma + dsigma
#         tagged_mass.setRange('tmp_range', center-tmp_sigma, center+tmp_sigma)
#         tmp_int = theRTgauss.createIntegral(RooArgSet(tagged_mass), RooFit.NormSet(RooArgSet(tagged_mass)), RooFit.Range('tmp_range')).getVal()
#         if tmp_int/tot_rt > 0.68:
#           print ' theRTgauss signal sigma:', tmp_sigma, ' for which we have ', tmp_int/tot_rt 
#           break
# 
#     set_tm = ROOT.RooArgSet(tagged_mass)
#     normSet = ROOT.RooFit.NormSet(set_tm)
#     inttot = theWTgauss.createIntegral(set_tm, normSet, ROOT.RooFit.Range('tmp_range1'))
#     int2   = theWTgauss.createIntegral(set_tm, normSet, ROOT.RooFit.Range('tmp_range2'))
    

    c1 = ROOT.TCanvas() 
    upperPad = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)

    upperPad.Draw(); lowerPad.Draw()
    upperPad.cd()
    frame.Draw()
    if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ], 0)
    leg = ROOT.TLegend( 0.7,0.65, 0.9,0.9, '','brNDC')
    extrastring_0 = 'fitfunction%s_Norm[tagged_mass]'%ibin
    extrastring_1 = '_Range[full]_NormRange[full]'
    leg.AddEntry(frame.findObject(extrastring_0+extrastring_1),   'Fit', 'l')
    leg.AddEntry(frame.findObject(extrastring_0+"_Comp[c_signalFunction%s]"%ibin+extrastring_1),  'total signal', 'l')
    leg.AddEntry(frame.findObject(extrastring_0+"_Comp[c_theRTgauss%s]"%ibin+extrastring_1),  'RT signal', 'l')
    leg.AddEntry(frame.findObject(extrastring_0+"_Comp[c_theWTgauss%s]"%ibin+extrastring_1),  'WT signal', 'l')
    leg.AddEntry(frame.findObject(extrastring_0+"_Comp[bkg_exp]"+extrastring_1),       'combinatorial bkg', 'l')
    frame.Draw()
    leg.Draw()

    ## add plot of pulls
    lowerPad.cd()
    hpull = frame.pullHist("h_fulldata", pdfstring)
    frame2 =  tagged_mass.frame(RooFit.Range("full"), RooFit.Title(''))
    frame2.addPlotable(hpull,"P") 
    niceFrameLowerPad(frame2, 'pull')
    frame2.Draw()
    line = ROOT.TLine(5.0,0,5.6,0)
    line.SetLineColor(ROOT.kGreen+3)
    line.Draw()

    for ilog in [False]:
        upperPad.SetLogy(ilog)
        c1.SaveAs('fit_results_mass/newbdt_puw/save_fit_data_%s_%s%s_fixBkg_deltaPeakConstr.pdf'%(ibin, args.year, '_logScale'*ilog))


    out_f.cd()
    r.Write('results_data_%s'%(ibin))

    params = fitFunction.getParameters(RooArgSet(tagged_mass)) 
    out_w.saveSnapshot("reference_fit_data_%s"%(ibin),params,ROOT.kTRUE) 
    getattr(out_w, 'import')(fitFunction)
#     getattr(out_w, 'import')(signalFunction)
#     getattr(out_w, 'import')(bkg_pol)






tData = ROOT.TChain('ntuple')

if args.year == 'test':
    tData.Add('/gwteray/users/fiorendi/final_ntuples_p5prime_allyears/2016Data_100k.root')
    fname_mcresults = 'fit_results_mass_checkOnMC/results_fits_2016_newSigmaFRT_Jpsi.root'
else:    
    tData.Add('/gwdata/y/users/fiorendi/final_ntuples_p5prime_allyears/%sdata_newphi_punzi_removeTkMu_fixBkg_B0Psicut_addxcutvariable.root'%args.year)

    fname_mcresults = '/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%s_fM.root'%args.year
    if args.dimusel == 'keepJpsi':
        fname_mcresults = '/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%s_fM_Jpsi.root'%args.year
    elif args.dimusel == 'keepPsiP':
        fname_mcresults = '/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%s_fM_Psi.root'%args.year


fo = ROOT.TFile()
try:
  fo = ROOT.TFile(fname_mcresults,'open')
except:
  print ('file %s or not found'%(fo))
w = fo.Get('w')

tagged_mass  = w.var("tagged_mass")

# tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.7, "GeV")
# tagged_mass     = RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 5., 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"       , "tagB0"    , 0, 2);
passB0Psi_lmnr  = RooRealVar("passB0Psi_lmnr" , "passB0Psi_lmnr", -200, 2);
passB0Psi_jpsi  = RooRealVar("passB0Psi_jpsi" , "passB0Psi_jpsi", -200, 2);
passB0Psi_psip  = RooRealVar("passB0Psi_psip" , "passB0Psi_psip", -200, 2);
xcut            = RooRealVar("xcut" , "xcut", -1, 5);

tagged_mass.setRange("full",   5.0,5.6) ;
thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)
thevars.add(passB0Psi_lmnr)
thevars.add(passB0Psi_jpsi)
thevars.add(passB0Psi_psip)
thevars.add(xcut)


fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))
out_f = TFile ("results_data_fits_%s%s_fM_fixBkg_deltaPeakConstr.root"%(args.year, '_Jpsi'*(args.dimusel=='keepJpsi') + '_Psi'*(args.dimusel=='keepPsiP')),"RECREATE") 
out_w = ROOT.RooWorkspace("data_w")

for ibin in range(len(q2binning)-1):

    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
    if args.dimusel == 'keepJpsi' and \
       (q2binning[ibin] != 8.68): 
           continue
    if args.dimusel == 'keepPsiP' and \
       (q2binning[ibin] != 12.86): 
           continue

    if   ibin ==4:  cut_base = 'passB0Psi_jpsi== 1 '
    elif ibin ==6:  cut_base = 'passB0Psi_psip== 1 '

    fitData(fulldata, ibin, w)
    print ' --------------------------------------------------------------------------------------------------- '


print '--------------------------------------------------------------------------------------------------- '
print 'bin\t\t fit status \t cov. matrix \t\t chi2'
for i,k in enumerate(fitStats.keys()):    
    if i%3==0:  print '------------------------------------------------------'
    print k , '\t\t', fitStats[k], '\t\t', covStats[k], '\t\t', chi2s[k]

out_f.Close()
out_w.writeToFile(out_f.GetName(), False)
