import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("inputfile" , nargs='*', help = "Path to the input ROOT file")
parser.add_argument("-e", "--era"     , dest = "era"    ,  help = "Define BH, GH, BF"                 , default='BH')
args = parser.parse_args()

import pdb
import ROOT
from ROOT import gSystem
import math

gSystem.Load('libRooFit')
from ROOT import RooFit, RooStats, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar, RooWorkspace, RooAbsData
from ROOT import RooGaussian, RooExponential, RooChebychev, TCanvas
import os, sys, inspect
from os import path
ROOT.gROOT.SetBatch(True)

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

sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')
sys.path.append("../utils")

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from utils.utils import *
from utils.fit_functions import *


from miniB0KstarMuMu.miniKstarMuMu.eras_allYears import *
era = args.era
if '2018' in args.inputfile[0]:
  era = '2018A,2018D'
elif '2017' in args.inputfile[0]:
#   era = 'DCDC,DCDC'
  era = '2017B,2017F'
elif '2016' in args.inputfile[0]:
  era = '2016B,2016H'
#   era = '2016B,2016F'
#   era = '2016G,2016H'
# pdb.set_trace()
print ('using full dataset per year: %s '%era)  

if era == '2016B,2016F':
  year = 'ALLAPV'
  year_string = ', 2016B-F'
elif era == '2016G,2016H':
  year = 'GOOD'
  year_string = ', 2016GH' 
elif era == '2016B,2016H':
  year = '2016'
  year_string = ', 2016' 
elif '2017' in era:
  year = '2017'
  year_string = ', 2017'
elif '2018' in era:
  year = '2018'
  year_string = ', 2018'
else:
  print 'era not defined!'  

firstRun = eras[era.split(',')[0]][0] - 1
lastRun  = eras[era.split(',')[1]][1] + 1

RooAbsData.setDefaultStorageType(RooAbsData.Tree)

print 'filename:', args.inputfile[0]
file = ROOT.TFile(args.inputfile[0], 'read')
print 'reading workspace'
ws = file.Get('w')
ws.Print()

data = ws.data('fulldata')
tagged_mass = ws.var('tagged_mass')
fitFunction = ws.pdf('fitFunction')

frame = tagged_mass.frame()
data.plotOn(frame, RooFit.Binning(100), RooFit.MarkerSize(.5))
fitFunction.plotOn(frame, );
fitFunction.plotOn(frame, RooFit.Components("bkg_pol"),      RooFit.LineStyle(ROOT.kDashed));
fitFunction.plotOn(frame, RooFit.Components("signalGauss"),  RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1));
fitFunction.plotOn(frame, RooFit.Components("gaus"),         RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));
fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+3));
frame.GetXaxis().SetTitle('#mu^{+}#mu^{-}K^{+}#pi^{-} mass [GeV]')
frame.SetTitle('sara')
frame.Draw()
if 'JPSI' in args.inputfile[0]:
  writeCMS(frame, year, [q2binning[4], q2binning[5]], 0, yearString = year_string)
elif 'PSI' in args.inputfile[0]:
  writeCMS(frame, year, [q2binning[6], q2binning[7]], 0, yearString = year_string)
niceFrame(frame, '')

leg = ROOT.TLegend( 0.67, 0.75, 0.88,0.88, '','brNDC')
extrastring_0 = 'fitfunction%s_Norm[tagged_mass]'#%ibin
extrastring_1 = '_Range[full]_NormRange[full]'
leg.AddEntry(frame.findObject('fitFunction_Norm[tagged_mass]_Comp[gaus]_Range[fit_nll_fitFunction_fulldata]_NormRange[fit_nll_fitFunction_fulldata]'),     'signal component', 'l')
leg.AddEntry(frame.findObject('fitFunction_Norm[tagged_mass]_Comp[bkg_pol]_Range[fit_nll_fitFunction_fulldata]_NormRange[fit_nll_fitFunction_fulldata]'),  'background component', 'l')

c1 = ROOT.TCanvas()
frame.Draw()
leg.Draw()
frame.getAttText().SetTextSize(0.027) 
frame.getAttText().SetTextFont(42) 

c1.SaveAs(args.inputfile[0].replace('out_workspace_', 'splot_mass_').replace('.root', '.pdf'))
