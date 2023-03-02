import argparse

parser = argparse.ArgumentParser(description="")
# parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("-c", "--channel" , dest = "channel",  help = "Define if JPSI or LMNR or PsiPrime", default='Jpsi')
parser.add_argument("-e", "--era"     , dest = "era"    ,  help = "Define BH, GH, BF"                 , default='BH')
parser.add_argument("-l", "--l1"      , dest = "l1"     , help = "L1 seed: l1_11_4, l1_12_5, l1_10_0, l1_00, l1_00_OS", default = 'all')
parser.add_argument("-y", "--year"    , dest = "year"   , help = "year", default = '2018')

args = parser.parse_args()

import ROOT
import sys, pdb
import math
from array import array
from collections import OrderedDict


B0Mass_   = 5.27963
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

b0_width   = 0.03601 
kaonMass_  = 0.493677

mc_sigma = 0.040
mc_mass  = 5.27783 

nSigma_psiRej = 3.

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)


class PlotContainer(object):
    def __init__(self, var, sel, mcsel, xtitle, norm, pltopt, nbins=0, xmin=0, xmax=0, mybins=0, label=None, logx=False, logy=False, fill=False):
        self.var         = var 
        self.sel         = sel 
        self.mcsel       = mcsel 
        self.xtitle      = xtitle 
        self.norm        = norm 
        self.nbins       = nbins 
        self.xmin        = xmin 
        self.xmax        = xmax
        self.pltopt      = pltopt
        self.logy        = logy
        self.logx        = logx
        self.fill        = fill
        self.mybins      = mybins 
        if label:
            self.label = label
        else:
            self.label = self.var

def setHistoProperties(h, color, xtitle, ytitle, xoff, yoff, markerSt, markerSize ):
    h.SetLineColor(color )
#     h.SetFillColor(color )
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.GetXaxis().SetTitleOffset(xoff)
    h.GetYaxis().SetTitleOffset(yoff)
    h.SetMarkerStyle(markerSt)
    h.SetMarkerSize(markerSize)
    h.SetMarkerColor(color)
    h.GetXaxis().SetLabelColor(ROOT.kWhite)
    h.GetXaxis().SetLabelSize(0.)
    h.SetNdivisions(0)

def set_baseline_histo(rp_th1):
    rp_th1 .SetTitle("");
    rp_th1 .SetDirectory(0);
    rp_th1 .SetStats(0);
    rp_th1 .GetXaxis().SetTitle(xtitle);
    rp_th1 .GetYaxis().SetTitle('data/MC');
    
    rp_th1 .GetXaxis().SetTitleSize(0.08);
    rp_th1 .GetYaxis().SetTitleSize(0.07);
    rp_th1 .GetYaxis().SetTitleOffset(0.7);
    rp_th1 .GetYaxis().SetNdivisions(4+800)
    rp_th1 .GetXaxis().SetTickLength(0.05)
    rp_th1 .GetXaxis().SetLabelSize( 0.07)
    rp_th1 .GetYaxis().SetLabelSize( 0.06)

    rp_th1.GetYaxis().SetRangeUser(0.05,2.1)
    rp_th1.SetMaximum(2.1)
    rp_th1.SetMinimum(0.05)



t0  = ROOT.TChain('ntuple')
t7  = ROOT.TChain('ntuple') ##mc

if args.channel=='Jpsi':

    t0 .Add('../final_ntuples/for_splot/%sdata_noIP2D%s_addxcutvariable_passSPlotCuts_mergeSweights.root'%(args.year, ''+'_noNan'*(args.year=='2017')))
    t7 .Add('/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/MC_with_xgbv8_weights/%sMC_JPSI_noIP2D%s_addxcutvariable_withMCw.root'%(args.year, ''+'_noNan'*(args.year=='2017')))

    psiselection = '((tagged_mass > 5 & tagged_mass < 5.6) && (mumuMass*mumuMass > 8.68 && mumuMass*mumuMass < 10.09) && pass_preselection == 1 && passB0Psi_jpsi == 1 && xcut == 0)'
    trig_index = 1

print 'data', t0.GetName()  
print 'MC', t7.GetName()  

selMC   = '((truthMatchMum == 1 && truthMatchMup == 1 && truthMatchTrkm == 1 && truthMatchTrkp == 1) &&\
            trig == {TR}) && \
            '.format(M=mc_mass,S=mc_sigma, TR=trig_index) + psiselection


selData = 'tagged_mass > 0'

pt_mu1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(10)]
pt_mu2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(1)] + [i*4+22 for i in range(3)]
pt_tk1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] + [i*2.5+15 for i in range(3)]
pt_tk2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] 

L_bins     = [i*0.01 for i in range(30)] + [i*0.05+0.3 for i in range(14)] + [i*0.1+1. for i in range(6)]
LE_bins    = [0] + [0.002+i*0.0004 for i in range(5)] + [i*0.0002+0.004 for i in range(20)] + [i*0.0004+0.008 for i in range(5)] + [i*0.001+0.01 for i in range(4)] + [i*0.002+0.014 for i in range(4)]
LS_bins    = [i*2 for i in range(50)]    + [i*4+100 for i in range(20)]   + [i*10+180. for i in range(3)]
DCA_bins   = [i*0.05 for i in range(50)] + [i*0.1+2.5 for i in range(10)] + [i*0.5+3.5 for i in range(4)]

cos_bins   = [0.994, 0.997, 0.998, 0.999, 0.9992, 0.9994, 0.9996, 0.9997, 0.9998, 0.9999, 0.99995, 1.]
pt_b_bins  = [i*10  for i in range(1)] + [i+10 for i in range(10)] + [i*2+20 for i in range(10)] + [i*5+40 for i in range(5)] + [i*10+70 for i in range(2)]


scales = [1]

toplot = [

    PlotContainer(var = 'bPt'                   , sel = selData, mcsel = selMC, xtitle = 'p_{T}(B^{0}) [GeV]'           , norm = True, nbins = 150 , xmin =  0            , xmax =  80  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bEta'                  , sel = selData, mcsel = selMC, xtitle = '#eta(B^{0})'                  , norm = True, nbins = 100 , xmin =  -2.5         , xmax =  2.5 , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bVtxCL'                , sel = selData, mcsel = selMC, xtitle = 'B^{0} vtx CL'                 , norm = True, nbins = 50  , xmin =  0            , xmax =  1   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstVtxCL'              , sel = selData, mcsel = selMC, xtitle = 'K* vtx CL'                    , norm = True, nbins = 50  , xmin =  0            , xmax =  1   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mu1Pt'                 , sel = selData, mcsel = selMC, xtitle = 'p_{T}(#mu_{1}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu1_bins[0], xmax =  pt_mu1_bins[-1] , mybins = pt_mu1_bins, fill = True, pltopt = 'HIST', label = 'leadingMuPt' ),
    PlotContainer(var = 'mu2Pt'                 , sel = selData, mcsel = selMC, xtitle = 'p_{T}(#mu_{2}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu2_bins[0], xmax =  pt_mu2_bins[-1] , mybins = pt_mu2_bins, fill = True, pltopt = 'HIST', label = 'trailingMuPt'),
    PlotContainer(var = 'kstTrk1Pt'             , sel = selData, mcsel = selMC, xtitle = 'leading trk p_{T} [GeV]'      , norm = True, nbins = 80  , xmin = pt_tk1_bins[0], xmax =  pt_tk1_bins[-1] , mybins = pt_tk1_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrk2Pt'             , sel = selData, mcsel = selMC, xtitle = 'trailing trk p_{T} [GeV]'     , norm = True, nbins = 80  , xmin = pt_tk2_bins[0], xmax =  pt_tk2_bins[-1] , mybins = pt_tk2_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mu1Eta'                , sel = selData, mcsel = selMC, xtitle = '#eta(#mu_{1})'                , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingMuEta' ),
    PlotContainer(var = 'mu2Eta'                , sel = selData, mcsel = selMC, xtitle = '#eta(#mu_{2})'                , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingMuEta'),
    PlotContainer(var = 'kstTrk1Eta'            , sel = selData, mcsel = selMC, xtitle = '#eta(tk_{1})'                 , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingTkEta' ),
    PlotContainer(var = 'kstTrk2Eta'            , sel = selData, mcsel = selMC, xtitle = '#eta(tk_{2})'                 , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingTkEta'),
    PlotContainer(var = 'bCosAlphaBS'           , sel = selData, mcsel = selMC, xtitle = 'cos#theta_{BS}'               , norm = True, nbins = 100 , xmin = cos_bins[0]  , xmax = cos_bins[-1], mybins = cos_bins, fill = True, pltopt = 'HIST', logy = True, label = 'bCosThetaBS'),
    PlotContainer(var = 'bLBS'                  , sel = selData, mcsel = selMC, xtitle = 'L_{BS} [cm]'                  , norm = True, nbins = 150 , xmin = L_bins[0]    , xmax = L_bins[-1]  , mybins = L_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bLBSE'                 , sel = selData, mcsel = selMC, xtitle = 'L_{BS} uncertainty [cm]'      , norm = True, nbins = 100 , xmin = LE_bins[0]   , xmax = LE_bins[-1] , mybins = LE_bins,  fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bLBS/bLBSE'            , sel = selData, mcsel = selMC, xtitle = 'L_{BS}/#sigma_{L}'            , norm = True, nbins = 100 , xmin = LS_bins[0]   , xmax = LS_bins[-1] , mybins = LS_bins,  fill = True, pltopt = 'HIST', label = 'LSigmaBS'),
    PlotContainer(var = 'bDCABS/bDCABSE'        , sel = selData, mcsel = selMC, xtitle = 'b DCA from BS significance'   , norm = True, nbins = 100 , xmin = DCA_bins[0]  , xmax = DCA_bins[-1], mybins = DCA_bins, fill = True, pltopt = 'HIST', label = 'bDCASign'),
    PlotContainer(var = 'kstarmass'             , sel = selData, mcsel = selMC, xtitle = 'kstar mass [GeV]'             , norm = True, nbins = 60  , xmin =  0.7         , xmax =  1.1 , fill = True, pltopt = 'HIST', label = 'kstMass'),
    PlotContainer(var = 'sum_isopt_04'          , sel = selData, mcsel = selMC, xtitle = 'B^{0} isolation'              , norm = True, nbins = 50  , xmin =  0           , xmax =  5   , fill = True, pltopt = 'HIST'),
#     PlotContainer(var = 'bdt_prob'              , sel = selData, mcsel = selMC, xtitle = 'bdt score'      , norm = True, nbins = 100 , xmin =  0 , xmax = 1 , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'cos_theta_k'           , sel = selData, mcsel = selMC, xtitle = 'cos#theta_{k}'   , norm = True, nbins = 100, xmin =  -1  , xmax =  1   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'cos_theta_l'           , sel = selData, mcsel = selMC, xtitle = 'cos#theta_{l}'   , norm = True, nbins = 100, xmin =  -1  , xmax =  1   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'phi_kst_mumu'          , sel = selData, mcsel = selMC, xtitle = '#phi'            , norm = True, nbins = 100, xmin =  -3.2, xmax =  3.2 , fill = True, pltopt = 'HIST'),
# 
]

def create_histo(name, nbins, xmin, xmax, mybins):
    h_data = ROOT.TH1F(name+'_%d'%i    , '', nbins, xmin, xmax)
    if mybins == 0:
      hr_data  = ROOT.TH1F(name+'_ratio_%d'%i   , '', nbins, xmin, xmax)
    else:
      hr_data  = ROOT.TH1F(name+'_ratio_%d'%i   , '', len(mybins)-1, array('d',mybins))
            
    return h_data, hr_data


ROOT.TH1.SetDefaultSumw2()
if args.channel == 'LMNR':
    colorlist = [ROOT.kBlack, ROOT.kBlack, ROOT.kMagenta+2, ROOT.kBlue-3]
if args.channel == 'Jpsi':
#     colorlist = [ROOT.kBlack, ROOT.kBlack, ROOT.kOrange-3, ROOT.kGreen+3]
    colorlist = [ROOT.kBlack, ROOT.kBlack, ROOT.kOrange-3, ROOT.kRed+2, ROOT.kMagenta+2, ROOT.kBlue-3, ROOT.kAzure+1, ROOT.kAzure+3] 
#     colorlist = [ROOT.kBlack, ROOT.kBlack, ROOT.kGreen+2, ROOT.kGray+2, ROOT.kMagenta+2, ROOT.kBlue-3, ROOT.kAzure+1, ROOT.kAzure+3] 

for i, iplot in enumerate(toplot):
    
    nbins  = iplot.nbins
    if args.channel == 'LMNR':
        nbins = nbins/2
    
    xmin   = iplot.xmin  
    xmax   = iplot.xmax  
    
    sel    = iplot.sel
    mcsel  = iplot.mcsel

    var    = iplot.var
    
    h_data, hr_data = create_histo('hdata', nbins, xmin, xmax, iplot.mybins)
    h_mc, hr_mc = create_histo('hmc', nbins, xmin, xmax, iplot.mybins)

    h_data3, hr_data3 = create_histo('hdata3', nbins, xmin, xmax, iplot.mybins)
    h_mc3, hr_mc3 = create_histo('hmc3', nbins, xmin, xmax, iplot.mybins)
    
    list_to_draw = OrderedDict()
    list_to_draw[h_data]  = 'B^{0}#rightarrow J/#psi K* data'
    list_to_draw[h_mc]    = 'B^{0}#rightarrow J/#psi K* MC'
#     list_to_draw[h_data3] = 'B^{0}#rightarrow J/#psi K* data'
    list_to_draw[h_mc3]   = 'B^{0}#rightarrow J/#psi K* MC, xgbv8'
    
    xtitle = iplot.xtitle
    ytitle = 'a.u.' if iplot.norm else 'events'
    
    setHistoProperties(h_data, colorlist[0], xtitle, ytitle, 1.2, 1.5,  8, 0.6 )
    setHistoProperties(h_mc, colorlist[2], xtitle, ytitle, 1.2, 1.5, 1, 0.6 )
    setHistoProperties(h_data3, colorlist[4], xtitle, ytitle, 1.2, 1.5,  8, 0.6 )
    setHistoProperties(h_mc3, colorlist[3], xtitle, ytitle, 1.2, 1.5, 1, 0.6 )
    
    t0.Draw('%s >> %s' %(var, h_data.GetName()), '(%s)*(nsig_sw)'  %(sel)              , iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, h_mc.GetName()),   '(%s)*(%s) '   %( mcsel,'weight'), iplot.pltopt   )
    t0.Draw('%s >> %s' %(var, h_data3.GetName()), '(%s)*(nsig_sw)'  %(sel)              , iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, h_mc3.GetName()),   '(%s)*(%s*%s) '   %( mcsel,'weight', 'MCw'), iplot.pltopt   )

    t0.Draw('%s >> %s' %(var, hr_data.GetName()),  '(%s)*(nsig_sw)'  %(sel)              , iplot.pltopt   ) # nsig_sw
    t7.Draw('%s >> %s' %(var, hr_mc.GetName()),    '(%s)*(%s) '   %( mcsel,'weight')  , iplot.pltopt   )
    t0.Draw('%s >> %s' %(var, hr_data3.GetName()),  '(%s)*(nsig_sw)'  %(sel)              , iplot.pltopt   ) # nsig_sw
    t7.Draw('%s >> %s' %(var, hr_mc3.GetName()),    '(%s)*(%s*%s) '   %( mcsel,'weight', 'MCw')  , iplot.pltopt   )

    max_y_all = []
    for ih in list_to_draw.keys():  max_y_all.append(ih.DrawNormalized().GetMaximum())
    ymax = max(max_y_all)
    print 'max y value: ' , ymax
    
    p_th1  = ROOT. TH1F("p_th1","",nbins, xmin, xmax)
    p_th1 .GetXaxis().SetLabelSize( 0.0)
    p_th1.GetYaxis().SetRangeUser(0.00001,1.25*ymax)
    p_th1.SetMaximum(1.25*ymax)
    
    c1 = ROOT.TCanvas('c1', 'c1', 700, 700)
    upperPad  = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad  = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.348 )  
    upperPad.Draw()
    lowerPad .Draw()
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)
   
    upperPad.cd()
    p_th1.Draw()
    if 'cos_theta' not in var and 'phi' not in var:
        h_data.DrawNormalized('e same')
#         h_data3.DrawNormalized('e same')
    h_mc.DrawNormalized('SAME e1 hist')
    h_mc3.DrawNormalized('SAME e1 hist')


    l1 = ROOT.TLegend(0.55,0.73,0.88,0.88)
#     l1 = ROOT.TLegend(0.55,0.6,0.88,0.88)
    for ih,jh in list_to_draw.items():
        l1.AddEntry(ih  ,  jh  , 'p'*('data' in jh) + 'l'*('MC' in jh) )
    l1.SetBorderSize(0)
    l1.SetTextSize(0.036)
    l1.Draw('same')

    c1.Update()
    c1.Modified()

    r_data_mc   = ROOT.TGraphAsymmErrors(hr_data.GetNbinsX())
    r_data3_mc3 = ROOT.TGraphAsymmErrors(hr_data3.GetNbinsX())

    if hr_data.Integral() != 0:
      hr_data      .Scale(1./hr_data.Integral())
    hr_mc   .Scale(1./hr_mc  .Integral())

    if hr_data3.Integral() != 0:
      hr_data3      .Scale(1./hr_data3.Integral())
    hr_mc3   .Scale(1./hr_mc3  .Integral())

    h_data_mc = hr_data.Clone()
    h_data_mc.Divide(hr_data, hr_mc, 1, 1, 'B')
    r_data_mc = ROOT.TGraphAsymmErrors(h_data_mc)

    h_data3_mc3 = hr_data3.Clone()
    if 'cos_theta' not in var and 'phi' not in var:
      h_data3_mc3.Divide(hr_data3, hr_mc3, 1, 1, 'B')
    else:
      h_data3_mc3.Divide(hr_mc, hr_mc3, 1, 1, 'B')
    r_data3_mc3 = ROOT.TGraphAsymmErrors(h_data3_mc3)

    r_data_mc.SetLineColor  (colorlist[2])
    r_data_mc.SetMarkerColor  (colorlist[2])

    r_data3_mc3.SetLineColor  (colorlist[3])
    r_data3_mc3.SetMarkerColor  (colorlist[3])

    rp_th1  = ROOT. TH1F("rp_th1","rp_th1",nbins,xmin,xmax)
    set_baseline_histo(rp_th1)
    if 'cos_theta' in var or 'phi'  in var:
      rp_th1 .GetYaxis().SetTitle('MC/MCw');

    lowerPad.cd()
    rp_th1. Draw("")
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(ROOT.kGreen+4)
    line.Draw()
    if 'cos_theta' not in var and 'phi' not in var:
      r_data_mc.Draw('same P')
    r_data3_mc3.Draw('same P')
#     else: 
#     r_data3_mc3.SetLineColor  (colorlist[4])
#     r_data3_mc3.SetMarkerColor  (colorlist[4])
#     r_data3_mc3.Draw('same P')
# 
    upperPad.cd()
    txt = ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.042)
    txt.DrawLatexNDC(.68, .91, '%s (13 TeV)'%args.year)   

    c1.Update()
    c1.Modified()

    c1.SaveAs('plots/MCw_xgbv8/%s_%s_%s_%s_afterMCw_xgbv8.pdf' %(args.year,args.channel, iplot.label, args.l1))
