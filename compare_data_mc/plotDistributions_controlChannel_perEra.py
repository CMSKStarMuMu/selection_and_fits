import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("year" , help = "year")
parser.add_argument("-c", "--channel" , dest = "channel",  help = "Define if Jpsi or LMNR or PsiPrime", default='Jpsi')
parser.add_argument("-e", "--era"     , dest = "era"    ,  help = "Define BH, GH, BF"                 , default='BH')
parser.add_argument("-l", "--l1"      , dest = "l1seed"     , help = "L1 seed: l1_11_4, l1_12_5, l1_10_0, l1_00, l1_00_OS", default = 'all')

args = parser.parse_args()

import ROOT
import sys
import math, pdb
from array import array


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
    h.SetFillColor(color )
    h.SetMarkerColor(color)
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



t0  = ROOT.TChain('ntuple')
t1  = ROOT.TChain('ntuple') ##data

t7  = ROOT.TChain('ntuple') ##mc
 
trig_index = 1
### trig = 0 for LMNR
### trig = 1 for Jpsi
### trig = 2 for PsiPrime


print 'still need to apply cut on l1 for MC!!!!!!!'
if args.channel=='Jpsi':

    if args.year == '2016' and args.l1seed == 'l1_00_OS':
      t0 .Add('../final_ntuples/for_splot/2016data_noIP2D_addxcutvariable_passSPlotCuts_2016BF_L100_mergeSweights.root')
      t1 .Add('../final_ntuples/for_splot/2016data_noIP2D_addxcutvariable_passSPlotCuts_2016GH_L100_mergeSweights.root')
      t7 .Add('../final_ntuples/2016MC_JPSI_noIP2D_addxcutvariable.root')
    elif args.year == '2016' and args.l1seed == 'l1_12_5':
      t0 .Add('../final_ntuples/for_splot/2016data_noIP2D_addxcutvariable_passSPlotCuts_2016BF_L1125_mergeSweights.root')
      t1 .Add('../final_ntuples/for_splot/2016data_noIP2D_addxcutvariable_passSPlotCuts_2016GH_L1125_mergeSweights.root')
      t7 .Add('../final_ntuples/2016MC_JPSI_noIP2D_addxcutvariable.root')
    else:
      t0 .Add('../final_ntuples/for_splot/%sdata_noIP2D%s_addxcutvariable_passSPlotCuts_mergeSweights.root'%(args.year, ''+'_noNan'*(args.year=='2017')))
      t1 .Add('../final_ntuples/for_splot/%sdata_noIP2D%s_addxcutvariable_passSPlotCuts_mergeSweights.root'%(args.year, ''+'_noNan'*(args.year=='2017')))
      t7 .Add('../final_ntuples/%sMC_JPSI_noIP2D%s_addxcutvariable.root'%(args.year, ''+'_noNan'*(args.year=='2017')))

    psiselection = '((tagged_mass > 5 & tagged_mass < 5.6) && (mumuMass*mumuMass > 8.68 && mumuMass*mumuMass < 10.09) && pass_preselection == 1 && passB0Psi_jpsi == 1 && xcut == 0)'
    trig_index = 1

elif args.channel=='LMNR':

   t1.Add('out_distribution_NRKstar_2017%s_2017%s.root'%(args.era[0], args.era[1]))
   t7.Add('../../../angular_analysis/final_ntuples/2017MC_LMNR_bdt0p96.root')
   psiselection = '&& ( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
                        (( mumuMass < {JPSIM} && !( abs(tagged_mass - {B0M} - (mumuMass - {JPSIM})) < 0.18 || abs(tagged_mass - {B0M} - (mumuMass - {PSIM})) < 0.0) ) || \
                         ( mumuMass > {PSIM}  && !( abs(tagged_mass - {B0M} - (mumuMass - {JPSIM})) < 0.0  || abs(tagged_mass - {B0M} - (mumuMass - {PSIM})) < 0.09) ) || \
                         ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( abs(tagged_mass - {B0M} - (mumuMass - {JPSIM})) < 0.08 || abs(tagged_mass - {B0M} - (mumuMass - {PSIM})) < 0.08 ))))'.format(B0M=B0Mass_, JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
   trig_index = 0

print 'data', t1.GetName()  
print 'MC', t7.GetName()  


pt_mu1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(10)]
pt_mu2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(1)] + [i*4+22 for i in range(3)]
pt_tk1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] + [i*2.5+15 for i in range(3)]
pt_tk2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] 

L_bins     = [i*0.01 for i in range(30)] + [i*0.05+0.3 for i in range(14)] + [i*0.1+1. for i in range(6)]
LE_bins    = [0] + [0.002+i*0.0004 for i in range(5)] + [i*0.0002+0.004 for i in range(20)] + [i*0.0004+0.008 for i in range(5)] + [i*0.001+0.01 for i in range(4)] + [i*0.002+0.014 for i in range(4)]
LS_bins    = [i*2 for i in range(50)]    + [i*4+100 for i in range(20)]   + [i*10+180. for i in range(3)]
DCA_bins   = [i*0.05 for i in range(50)] + [i*0.1+2.5 for i in range(10)] + [i*0.5+3.5 for i in range(4)]

cos_bins   = [0.994, 0.997, 0.998, 0.999, 0.9992, 0.9994, 0.9996, 0.9997, 0.9998, 0.9999, 0.99995, 1.]

selMC   = '((truthMatchMum == 1 && truthMatchMup == 1 && truthMatchTrkm == 1 && truthMatchTrkp == 1) &&\
            trig == {TR}) && \
            '.format(M=mc_mass,S=mc_sigma, TR=trig_index) + psiselection

if args.l1seed != 'all':
    selMC = selMC + ' && ({L1} == 1)'.format(L1=args.l1seed)
selData = 'tagged_mass > 0'
              
# import pdb
# pdb.set_trace()   
# selData = '((tagB0==1 && (bMass    > {M}-2.5*{S} && bMass    < {M}+2.5*{S}) ) || \
#             (tagB0==0 && (bBarMass > {M}-2.5*{S} && bBarMass < {M}+2.5*{S}) ))\
#           '.format(M=mc_mass,S=mc_sigma)
#           (kstMass*tagB0 + kstBarMass*(1-tagB0) > 0.806 && kstMass*tagB0 + kstBarMass*(1-tagB0) < 0.986 )  '.format(M=mc_mass,S=mc_sigma)


thebmass     = 'bMass*tagB0 + bBarMass*(1-tagB0)'
thekstmass   = 'kstMass*tagB0 + kstBarMass*(1-tagB0)'


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
#     PlotContainer(var = 'kstTrk1DCABS'          , sel = selData, mcsel = selMC, xtitle = 'leading trk DCA(BS) [cm]'     , norm = True, nbins = 80  , xmin =  -0.4        , xmax =  0.4            ,                    fill = True, pltopt = 'HIST'),
#     PlotContainer(var = 'kstTrk2DCABS'          , sel = selData, mcsel = selMC, xtitle = 'trailing trk DCA(BS) [cm]'    , norm = True, nbins = 80  , xmin =  -0.4         , xmax =  0.4            ,                    fill = True, pltopt = 'HIST'),
]



ROOT.TH1.SetDefaultSumw2()
if args.channel == 'LMNR':
    colorlist = [ROOT.kBlack, ROOT.kBlack, ROOT.kMagenta+2, ROOT.kBlue-3]
if args.channel == 'Jpsi':
    if args.year == '2016':
        colorlist = [ROOT.kBlack, ROOT.kGray, ROOT.kOrange-3, ROOT.kGreen+3,ROOT.kCyan-2, ROOT.kBlue-3, ROOT.kViolet+5, ROOT.kMagenta+1, ROOT.kSpring-5]
    elif args.year == '2017':
        colorlist = [ROOT.kBlack, ROOT.kGray, ROOT.kAzure+1, ROOT.kAzure+1, ROOT.kViolet+5, ROOT.kMagenta+1, ROOT.kSpring-5]
    elif args.year == '2018':
        colorlist = [ROOT.kBlack, ROOT.kGray, ROOT.kGreen+3, ROOT.kGreen+3]



for i, iplot in enumerate(toplot):
    
    nbins  = iplot.nbins
    if args.channel == 'LMNR':
        nbins = nbins/2
    
    xmin   = iplot.xmin  
    xmax   = iplot.xmax  
    
    sel    = iplot.sel
    mcsel  = iplot.mcsel

    var    = iplot.var
    
    h_b       = ROOT.TH1F('h_b_%d'%i     , '', nbins, xmin, xmax)
    h_c       = ROOT.TH1F('h_c_%d'%i     , '', nbins, xmin, xmax)
    h_mc_bf   = ROOT.TH1F('h_mc_bf_%d'%i , '', nbins, xmin, xmax)
    h_mc_gh   = ROOT.TH1F('h_mc_gh_%d'%i , '', nbins, xmin, xmax)
    if iplot.mybins == 0:
        hr_b       = ROOT.TH1F('hr_b_%d'%i    , '', nbins, xmin, xmax)
        hr_c       = ROOT.TH1F('hr_c_%d'%i    , '', nbins, xmin, xmax)
        hr_mc_bf   = ROOT.TH1F('hr_mc_bf_%d'%i, '', nbins, xmin, xmax)
        hr_mc_gh   = ROOT.TH1F('hr_mc_gh_%d'%i, '', nbins, xmin, xmax)
    else:
        hr_b       = ROOT.TH1F('hr_b_%d'%i    , '', len(iplot.mybins)-1, array('d',iplot.mybins))
        hr_c       = ROOT.TH1F('hr_c_%d'%i    , '', len(iplot.mybins)-1, array('d',iplot.mybins))
        hr_mc_bf   = ROOT.TH1F('hr_mc_bf_%d'%i, '', len(iplot.mybins)-1, array('d',iplot.mybins))
        hr_mc_gh   = ROOT.TH1F('hr_mc_gh_%d'%i, '', len(iplot.mybins)-1, array('d',iplot.mybins))

    xtitle = iplot.xtitle
    ytitle = 'a.u.' if iplot.norm else 'events'
    
    
    p_th1  = ROOT. TH1F("rp_th1","rp_th1",nbins,xmin,xmax)
    setHistoProperties(h_b    , colorlist[0], xtitle, ytitle, 1.2, 1.5, 20, 0.6 )
    setHistoProperties(h_c    , colorlist[1], xtitle, ytitle, 1.2, 1.5, 24, 0.6 )
    setHistoProperties(h_mc_bf, colorlist[2], xtitle, ytitle, 1.2, 1.5, 1, 0.06 )
    setHistoProperties(h_mc_gh, colorlist[3], xtitle, ytitle, 1.2, 1.5, 1, 0.06 )
    
    era_string_b = ''
    era_string_c = ''
    era_string_b_mc = ''
    era_string_c_mc = ''
    weight_string_bf = 'weight'
    weight_string_gh = 'weight'
    if args.year == '2016':
        era_string_b_mc = ', as BF'
        era_string_c_mc = ', asGH'
        era_string_b = ', 2016BF'
        era_string_c = ', 2016GH'
        weight_string_bf = 'weightBF'
        weight_string_gh = 'weightGH'
    
#     t1.Draw('%s >> %s' %(var, h_b.GetName()), '(%s)' %(sel)   , iplot.pltopt   )
#     t2.Draw('%s >> %s' %(var, h_mc  .GetName()), '(%s)' %( mcsel), iplot.pltopt + 'SAME')
    t0.Draw('%s >> %s' %(var, h_b.GetName())    , '(%s)*(nsig_sw)'   %(sel)              , iplot.pltopt   )
    t1.Draw('%s >> %s' %(var, h_c.GetName())    , '(%s)*(nsig_sw)'   %(sel)              , iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, h_mc_bf.GetName()), '(%s )*(%s) '      %( mcsel,'weight'), iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, h_mc_gh.GetName()), '(%s )*(%s) '      %( mcsel,'weight'), iplot.pltopt   )

    t0.Draw('%s >> %s' %(var, hr_b.GetName())    , '(%s)*(nsig_sw)'   %(sel)              , iplot.pltopt   )
    t1.Draw('%s >> %s' %(var, hr_c.GetName())    , '(%s)*(nsig_sw)'   %(sel)              , iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, hr_mc_bf.GetName()), '(%s )*(%s) '      %( mcsel,'weight'), iplot.pltopt   )
    t7.Draw('%s >> %s' %(var, hr_mc_gh.GetName()), '(%s )*(%s) '      %( mcsel,'weight'), iplot.pltopt   )

#     import pdb; pdb.set_trace()
    ymax = max(h_b.DrawNormalized().GetMaximum(), h_c.DrawNormalized().GetMaximum(), h_mc_gh.DrawNormalized().GetMaximum(), h_mc_bf.DrawNormalized().GetMaximum())
    print 'max y value: ' , ymax
    p_th1 .SetTitle("");
    p_th1.GetYaxis().SetRangeUser(0.00001,1.15*ymax)
    p_th1.SetMaximum(1.2*ymax)
    p_th1.GetXaxis().SetLabelColor(ROOT.kWhite)
#     h_b.SetMinimum(0.)
#     h_b.SetMaximum(1.1*ymax)
#     h_mc  .SetMaximum(1.1*ymax)
    
       
    c1 = ROOT.TCanvas('c1', 'c1', 700, 700)
    upperPad  = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad  = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.Draw()
    lowerPad .Draw()
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)
   
    upperPad.cd()
    p_th1.Draw()
#     h_b.GetXaxis().SetLabelColor(ROOT.kWhite)
#     h_mc.GetXaxis().SetLabelColor(ROOT.kWhite)
#     h_b.GetXaxis().SetLabelSize(0.)
#     h_mc.GetXaxis().SetLabelSize(0)
# 

    h_b .DrawNormalized('e same')
    if args.year == '2016':
        h_c .DrawNormalized('same e')

    h_mc_bf.SetFillStyle(3002)
    h_mc_bf.DrawNormalized('SAME e2')
    h_mc_bf.SetFillColor(ROOT.kWhite  )
    h_mc_bf.DrawNormalized('SAME hist')

    h_mc_gh.SetFillStyle(3004)
    h_mc_gh.DrawNormalized('SAME e2')
    h_mc_gh.SetFillColor(ROOT.kWhite  )
    if args.year == '2016':
        h_c .DrawNormalized('same e')
    h_mc_gh.DrawNormalized('SAME hist')
#     h_c .DrawNormalized('SAME hist')

    l1 = ROOT.TLegend(0.55,0.73,0.88,0.88)
    if args.year == '2016':
        l1 = ROOT.TLegend(0.55,0.68,0.88,0.88)
    if args.channel=='LMNR':
        l1.AddEntry(h_b     , 'B^{0}#rightarrow J/#psi K* data' + era_string_b  , 'p')
        l1.AddEntry(h_c     , 'B^{0}#rightarrow J/#psi K* data' + era_string_c  , 'p')
        l1.AddEntry(h_mc_bf , 'B^{0}#rightarrow #mu^{+}#mu^{-} K* MC' + era_string_b_mc , 'f')
        l1.AddEntry(h_mc_gh , 'B^{0}#rightarrow #mu^{+}#mu^{-} K* MC' + era_string_c_mc , 'f')
    elif args.channel=='Jpsi':
        l1.AddEntry(h_b     , 'B^{0}#rightarrow J/#psi K* data' + era_string_b  , 'p')
        if args.year == '2016':
            l1.AddEntry(h_c     , 'B^{0}#rightarrow J/#psi K* data' + era_string_c  , 'p')
        l1.AddEntry(h_mc_bf , 'B^{0}#rightarrow J/#psi K* MC ' + era_string_b_mc, 'f')
        if args.year == '2016':
            l1.AddEntry(h_mc_gh , 'B^{0}#rightarrow J/#psi K* MC '+ era_string_c_mc, 'f')

    l1.SetBorderSize(0)
    l1.SetTextSize(0.036)
    l1.Draw('same')

# #     ROOT.gPad.SetLogx(iplot.logx)
# #     ROOT.gPad.SetLogy(iplot.logy)
# #     ROOT.gPad.Update()

    c1.Update()
    c1.Modified()

    r_b  = ROOT.TGraphAsymmErrors(hr_b.GetNbinsX())
    r_c  = ROOT.TGraphAsymmErrors(hr_c.GetNbinsX())
    hr_b    .Scale(1./hr_b.Integral())
    hr_c    .Scale(1./hr_c.Integral())
    hr_mc_bf .Scale(1./hr_mc_bf  .Integral())
    hr_mc_gh .Scale(1./hr_mc_gh  .Integral())

    hr_data_mc_b = hr_b.Clone()
    hr_data_mc_b.Divide(hr_b, hr_mc_bf, 1, 1, 'B')
    r_b = ROOT.TGraphAsymmErrors(hr_data_mc_b)

    hr_data_mc_c = hr_b.Clone()
    hr_data_mc_c.Divide(hr_c, hr_mc_gh, 1, 1, 'B')
    r_c = ROOT.TGraphAsymmErrors(hr_data_mc_c)
    
#     r_b .Divide( hr_b , hr_mc_bf, 'pois v' )
#     r_c .Divide( hr_c , hr_mc_gh, 'pois v' )
    
    r_b .SetMarkerStyle(8)
    r_c .SetMarkerStyle(4)
    r_b .SetMarkerColor(colorlist[2])
    r_c .SetMarkerColor(colorlist[3])

    r_b .SetLineColor  (colorlist[2])
    r_c .SetLineColor  (colorlist[3])

    rp_th1  = ROOT. TH1F("rp_th1","rp_th1",nbins,xmin,xmax)
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

#     rp . SetHistogram(rp_th1);   
    rp_th1.GetYaxis().SetRangeUser(0.05,2.1)
    rp_th1.SetMaximum(2.1)
    rp_th1.SetMinimum(0.05)
#     rp_th1.GetYaxis().SetMaximum(0.5,2)

# # #     rp.GetLowYaxis().SetNdivisions(505)
    lowerPad.cd()
    rp_th1. Draw("")
    line =ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(ROOT.kGreen+4)
    line.Draw()
    r_b .Draw('same P')
    if args.year == '2016':
        r_c .Draw('same P')

    upperPad.cd()
    txt = ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.042)
    txt.DrawLatexNDC(.73, .91, '%s (13 TeV)'%args.year)   

    c1.Update()
    c1.Modified()

#     c1.SaveAs('plots/Mar21ntuples_BDTSelection_%s/%s_%s_perSingleEra_8b4e.pdf' %(args.channel, iplot.label, args.l1))
    c1.SaveAs('plots/compare_for_AN/Jpsi_%s_%s_%s.pdf' %(iplot.label, args.l1seed, args.year))
#     c1.SaveAs('plots/compare_for_AN/%s_%s_perDCDC.png' %(iplot.label, args.l1))
