import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("inputfile" , nargs='*', help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: Jpsi, LMNR")
parser.add_argument("-d", "--doubleg"  , dest = "doubleg" , help = "Define if use 1 or 2 gaussian for the signal model"   , default = '2')
parser.add_argument("-e", "--era"      , dest = "era"     , help = "Define which data taking period to use (2017B,2017F)" , default = '2018A,2018D')
parser.add_argument("-r", "--range"    , dest = "range"   , help = "Define B mass range"                                  , default = '5.0,5.6')
parser.add_argument("-l", "--l1seed"   , dest = "l1seed"  , help = "L1 seed: l1_11_4, l1_12_5, l1_10_0, l1_00, l1_00_OS"  , default = 'all')
parser.add_argument("--save_workspace"      , help = "", action='store_true')

args = parser.parse_args()

import pdb
import ROOT
from ROOT import gSystem

gSystem.Load('libRooFit')
from ROOT import RooFit, RooStats, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar, RooWorkspace, RooAbsData
from ROOT import RooGaussian, RooExponential, RooChebychev, TCanvas
import sys
import math

if args.save_workspace:
  print ('will exit after saving the workspace!!!!!')  

print ('FIXME: Need to manually set the desired era')  
from miniB0KstarMuMu.miniKstarMuMu.eras_allYears import *
era = args.era
if '2018' in args.inputfile[0]:
  era = '2018A,2018D'
elif '2017' in args.inputfile[0]:
#   era = 'DCDC,DCDC'
  era = '2017B,2017F'
elif '2016' in args.inputfile[0]:
  era = '2016B,2016H'
#   era = '2016G,2016H'
#   era = '2016B,2016H'
# pdb.set_trace()
print ('using full dataset per year: %s '%era)  


firstRun = eras[era.split(',')[0]][0] - 1
lastRun  = eras[era.split(',')[1]][1] + 1

lowRange  = float( args.range.split(',')[0] )
highRange = float( args.range.split(',')[1] )


B0Mass_   = 5.27958
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

B0Mass     = RooRealVar("B0Mass"    , "B0Mass"  , B0Mass_  )
JPsiMass   = RooRealVar("JPsiMass"  , "JPsiMass", JPsiMass_)
PsiPMass   = RooRealVar("PsiPMass"  , "PsiPMass", PsiPMass_)
KStMass    = RooRealVar("KStMass"   , "KStMass" , KStMass_ )

nSigma_psiRej = 3.


RooAbsData.setDefaultStorageType(RooAbsData.Tree)

print 'reading ntuple'
tree = ROOT.TChain('ntuple')
for filename in args.inputfile:
# for filename in args.inputfile.split(' '):
    print(filename)
    tree.AddFile(filename)

print 'loaded'

tagged_mass = RooRealVar("tagged_mass"    , "tagged_mass", 2, 10, "GeV")
mumuMass    = RooRealVar("mumuMass" , "mumuMass" , 0, 6)
mumuMassE   = RooRealVar("mumuMassE", "mumuMassE", 0, 10000)
tagB0       = RooRealVar("tagB0"    , "tagB0"    , 0, 2)
runN        = RooRealVar("runN"     , "runN"     , firstRun, lastRun)  
eventD      = RooRealVar("eventD"   , "eventD"     , 0, 99999999999)  

l1_11_4   = RooRealVar("l1_11_4"  , "l1_11_4"     , -2, 2000)  
l1_12_5   = RooRealVar("l1_12_5"  , "l1_12_5"     , -2, 2000)  
l1_10_0   = RooRealVar("l1_10_0"  , "l1_10_0"     , -2, 2000)  
l1_00     = RooRealVar("l1_00"    , "l1_00"       , -2, 2000)  
l1_00_OS  = RooRealVar("l1_00_OS" , "l1_00_OS"    , -2, 2000)  
# 
# ## variables for pre and bdt-selection
pass_preselection    = RooRealVar("pass_preselection"   , "pass_preselection",  1,  2);
passB0Psi_jpsi       = RooRealVar("passB0Psi_jpsi"   , "passB0Psi_jpsi",  0,  3);
passB0Psi_psip       = RooRealVar("passB0Psi_psip"   , "passB0Psi_psip",  0,  3);
xcut                 = RooRealVar("xcut"   , "xcut",  -1,  2);
bdt_prob             = RooRealVar("bdt_prob"            , "bdt_prob"         , -2,  2);

bPt                 = RooRealVar('bPt','bPt',                              0,1000)
bEta                = RooRealVar('bEta','bEta',                            -3,3)
kstTrk1Pt           = RooRealVar('kstTrk1Pt','kstTrk1Pt',                  0, 1000)


print 'preparing rooArgSet'

thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)
thevars.add(runN)
thevars.add(eventD)
thevars.add(passB0Psi_jpsi)
thevars.add(passB0Psi_psip )
thevars.add(kstTrk1Pt )
# thevars.add(l1_10_0 )
# thevars.add(l1_00   )
thevars.add(l1_00_OS)
thevars.add(l1_11_4)
thevars.add(l1_12_5)
# 
thevars.add(pass_preselection)
thevars.add(bdt_prob         )

thevars.add(bPt  )
thevars.add(bEta )
thevars.add(xcut )


fulldata   = RooDataSet('fulldata', 'fulldata', tree,  RooArgSet(thevars))
tagged_mass.setRange(lowRange,highRange);


# frame = tagged_mass.frame()
# fulldata.plotOn(frame, RooFit.Binning(100), RooFit.MarkerSize(.5))
# c1 = ROOT.TCanvas()
# frame.Draw()
# c1.SaveAs('fitJpsiKstar_forSplot_%s_%s_%s_preBDT_mcweights.pdf'%(era.split(',')[0],era.split(',')[1],args.l1seed))




if args.dimusel == 'keepJpsi':
  cut = '(mumuMass*mumuMass>8.68 && mumuMass*mumuMass < 10.09) && \
          pass_preselection==1 && \
          (tagged_mass > 5.0 && tagged_mass < 5.6) && \
          xcut == 0 && \
          passB0Psi_jpsi == 1 ' #&& \
#           kstTrk1Pt > 1.2'
elif args.dimusel == 'keepPsiP':
  cut = '(mumuMass*mumuMass>12.86 && mumuMass*mumuMass < 14.18) && \
          pass_preselection==1 && \
          (tagged_mass > 5.0 && tagged_mass < 5.6) && \
          passB0Psi_psip == 1'
# elif args.dimusel == 'rejectPsi':
#   cut = 'passB0Psi_lmnr == 1'
# elif args.dimusel == 'keepPsi':
#   cut = '(passB0Psi_jpsi==1 || passB0Psi_psi==1)'
# elif args.dimusel == 'nocut':
#   cut = 'mumuMass > 0'
else:
  print '\nYou should define which dimuon mass to consider. Please choose between following options: \nkeepPsiP, keepJpsi, rejectPsi, keepPsi'
  sys.exit(0)

  
if 'all' not in args.l1seed:
   cut = cut + '&& %s == 1'%args.l1seed 

print cut    
data = fulldata.reduce(thevars, cut)
print data.sumEntries()
## mass model 
mean        = RooRealVar ("mass"          , "mean"          ,  B0Mass_,   5.24,    5.3, "GeV")
sigma       = RooRealVar ("sigma"         , "sigma"         ,  0.024,     0,   0.04, "GeV")
signalGauss = RooGaussian("signalGauss"   , "signal gauss"  ,  tagged_mass,  mean,sigma)

sigma2       = RooRealVar ("sigma2"      , "sigma2"         ,  0.046,     0.03,   0.09, "GeV")
signalGauss2 = RooGaussian("signalGauss2" , "signal gauss2"  ,  tagged_mass,  mean,sigma2)
f1           = RooRealVar ("f1"           , "f1"             ,  0.4  ,     0.,   1.)
gaus         = RooAddPdf  ("gaus"         , "gaus1+gaus2"    , RooArgList(signalGauss,signalGauss2), RooArgList(f1))

## make bkg model
pol_c1      = RooRealVar ("p1"           , "coeff x^0 term",    -1,   -10, 10);
pol_c2      = RooRealVar ("p2"           , "coeff x^1 term",    0.2,   -10, 10);
pol_c3      = RooRealVar ("p3"           , "coeff x^2 term",    0.,   -10, 10);
bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  tagged_mass, RooArgList(pol_c1,pol_c2));
if args.dimusel == 'keepPsiP':
  bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  tagged_mass, RooArgList(pol_c1,pol_c2, pol_c3));
slope         = RooRealVar    ("slope"   , "slope"           ,   -4,   -10, 10);
bkg_exp       = RooExponential("bkg_exp" , "exponential"     ,  slope,   tagged_mass  );

## combined model
nsig        = RooRealVar("nsig"           , "signal frac"   ,    1700000,     0,    5000000)
nbkg        = RooRealVar("nbkg"           , "bkg fraction"  ,    1700000,     0,    10000000)
fitFunction = RooAddPdf ("fitFunction"    , "fit function"  ,  RooArgList(gaus, bkg_pol), RooArgList(nsig, nbkg))
# fitFunction = RooAddPdf ("fitfunction"    , "fit function"  ,  RooArgList(signalGauss, bkg_exp), RooArgList(nsig, nbkg))
if args.doubleg=='2':
    fitFunction = RooAddPdf ("fitFunction"    , "fit function"  ,  RooArgList(gaus, bkg_pol), RooArgList(nsig, nbkg))


print 'Calculate sWeights'
## fit the model to the data.
r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(lowRange,highRange))

frame = tagged_mass.frame()
data.plotOn(frame, RooFit.Binning(100), RooFit.MarkerSize(.5))
fitFunction.plotOn(frame, );
fitFunction.plotOn(frame, RooFit.Components("bkg_exp"),      RooFit.LineStyle(ROOT.kDashed));
fitFunction.plotOn(frame, RooFit.Components("signalGauss"),  RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1));
if args.doubleg=='2':
    fitFunction.plotOn(frame, RooFit.Components("gaus"),         RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));
    fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+3));
fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))
frame.SetTitle('')
frame.GetYaxis().SetTitleOffset(1.35)
frame.getAttText().SetTextSize(0.022) 
frame.getAttText().SetTextFont(42) 
frame.getAttLine().SetLineColor(0) 
c1 = ROOT.TCanvas()
frame.Draw()
c1.SaveAs('fit%sKstar_forSplot_%s_%s_%s_postBDT.pdf'%(('Jpsi'*(args.dimusel=='keepJpsi') + 'Psi'*(args.dimusel=='keepPsiP')), 
                                                             era.split(',')[0],era.split(',')[1],
                                                             args.l1seed))

if args.save_workspace:
    out_fname = "out_workspace_%s_%s_%s_L1%s_postBDT.root"%(('JPSI'*(args.dimusel=='keepJpsi') + 'PSI'*(args.dimusel=='keepPsiP'), 
                                                                  era.split(',')[0],era.split(',')[1],args.l1seed))
    
    ws = ROOT.RooWorkspace("w")
    params = fitFunction.getParameters(RooArgSet(tagged_mass))
    ws.saveSnapshot("reference_fit_data",params,ROOT.kTRUE)
    getattr(ws,"import")(fitFunction)
    getattr(ws,"import")(data)
    ws.writeToFile(out_fname, False)
    exit(0)


## Now we use the SPlot class to add SWeights to our data set based on our model and our yield variables
sData = RooStats.SPlot("sData","An SPlot", data, fitFunction, RooArgList(nsig, nbkg) )

## Check that our weights have the desired properties
print 'Check SWeights:'
print 'Yield of B0 is ' , nsig.getVal(), '.  From sWeights it is ', sData.GetYieldFromSWeight('nsig')
print 'Yield of bkg is ', nbkg.getVal(), '.  From sWeights it is ', sData.GetYieldFromSWeight('nbkg')

outfile = ROOT.TFile("out_distribution_%s_%s_%s_L1%s_postBDT.root"%(('JPSI'*(args.dimusel=='keepJpsi') + 'PSI'*(args.dimusel=='keepPsiP')), 
                                                                     era.split(',')[0], 
                                                                     era.split(',')[1], 
                                                                     args.l1seed), 
                                                                     "recreate")
outfile . cd();
thetree = data.GetClonedTree() ### was .tree()
print 'cloned'
thetree . Write();
print 'written'
# outfile . Close();
# 'prind closed'


