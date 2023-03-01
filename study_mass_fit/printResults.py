import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("bin", help = "choose q2 bin range", default = -1, type = float)
parser.add_argument("year", help = "choose year [format:2016, 20162017]", default = '2016')
args = parser.parse_args()

import ROOT
from ROOT import RooFit
from math import sqrt, sin, cos
import math
import itertools, pdb
from pdb import set_trace
from uncertainties import ufloat
from uncertainties.umath import sqrt as usqrt
from uncertainties.umath import sin as usin
from uncertainties.umath import cos as ucos
import re
import pandas as pd
import numpy as np

from collections import OrderedDict
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('../utils/func_roofit/libRooDoubleCBFast')

class fit_pars(object):
    '''
    '''
    def __init__(self):
        self.Reset()

    def Reset(self):
        self.val          = ufloat(-10,0)
        self.year         = 99
        self.bin          = 99
        self.name         = 'parname'

#     def __str__(self):
#         toWrite = ' fRT MC  = %f +/- %f \n \
# fRT fit  = %f +/- %f \n \
# deltaM   = %f +/- %f \n '%(self.fRT_mc.n, self.fRT.s, self.fRT.n - self.fRT_mc)
#         return toWrite 


my_year = args.year
bin_range = [0,1,2,3,5,7]
if args.bin != -1:
    bin_range = [int(args.bin)]  

from collections import OrderedDict
fnames = OrderedDict({
  '/eos/cms/store/user/fiorendi/p5prime/massFits/noIP2D/results_fits_%s_fM_Psi_noIP2D.root'%my_year : 'nominal',
  '/eos/cms/store/user/fiorendi/p5prime/massFits/noIP2D/xgbv8/results_fits_%s_fM_Psi_noIP2D_MCw_xgbv8.root'%my_year : 'mcw',
})
# fname_data = 'results_data_fits_2018_newSigmaFRT.root'

fr = []
tags = []
for i,j in fnames.items():
#     pdb.set_trace()
    year = int(re.findall(r'-?\d+', i)[1]) ## was 0
    print 'opening ', i
    fo = ROOT.TFile(i,'open')
    keys = fo.GetListOfKeys()
    for k in keys:  
        if 'results_' in k.GetName():   
            fr.append(fo.Get(k.GetName()))
            fr[-1].SetName(k.GetName() + '_' + j + '_' + str(year) ) 
            tags.append(j)

parlist = [
#            'Yield', 
#            'nbkg', 
           'mean_{RT}',
           'mean_{WT}',
           '#sigma_{RT1}',
           '#sigma_{RT2}',
           '#alpha_{RT1}',
           'n_{RT1}',
           'f^{RT}',
           '#sigma_{WT1}',
           '#alpha_{WT1}',
           '#alpha_{WT2}',
           'n_{WT1}',
           'n_{WT2}',
           
           ]

# format_list = [parname, year, bin, value, error]


value_list = []
# par_dict = OrderedDict()
for ipar in parlist:
    for i,ifr in enumerate(fr):
#         pdb.set_trace()
        if 'Yield' not in ipar and 'nbkg' not in ipar:
          parname = ipar + '^{%s}'%ifr.GetName().split('_')[2]
          if 'mean_{WT}' in ipar and 'deltaPeak' in tags[i]:
            parname = "deltaPeakVar%s"%ifr.GetName().split('_')[2]

        else:
          parname = ipar  
#         print parname  
        
        getpar = ifr.floatParsFinal().find(parname)
        if getpar == None:
          continue

#         print int(ifr.GetName().split('_')[-1]), int(ifr.GetName().split('_')[-2]), '\t'  , ifr.floatParsFinal().find(parname).Print() 
        tmp_list = [ ipar, 
                     int(ifr.GetName().split('_')[-1]), 
                     int(ifr.GetName().split('_')[2]),
                     ufloat(ifr.floatParsFinal().find( parname).getVal(),
                     ifr.floatParsFinal().find( parname).getError()),
                     tags[i]
                     ] 
        if 'mean_{WT}' in ipar and 'deltaPeak' in tags[i]:
          tmp_list = [ ipar, 
                       int(ifr.GetName().split('_')[-1]), 
                       int(ifr.GetName().split('_')[2]),
                       ufloat(-ifr.floatParsFinal().find( parname).getVal()+ifr.floatParsFinal().find('mean_{RT}'+ '^{%s}'%ifr.GetName().split('_')[2]).getVal(),
                       ifr.floatParsFinal().find( parname).getError()),
                       tags[i]
                       ] 
        value_list.append(tmp_list)             
        
df = pd.DataFrame(value_list)
df.rename(columns={0 : 'parname', 
                   1 : 'year',
                   2 : 'bin',
                   3 : 'val',
#                    4 : 'error',
                   4 : 'tag',
                   }, inplace=True)

        

this_tag = 'nominal'
mask_tag  = (df['tag'] == this_tag)
# df2 = df[mask_tag][['parname', 'year', 'bin', 'val', 'error']]
df2 = df[mask_tag][['parname', 'year', 'bin', 'val']]
df2 = df2.rename(columns={ "val": "val_%s"%this_tag,\
#                            "error": "err_%s"%this_tag\
                         })

next_tag = 'mcw'
mask_tag  = (df['tag'] == next_tag)
# df3 = df[mask_tag][['parname', 'year', 'bin', 'val', 'error']]
df3 = df[mask_tag][['parname', 'year', 'bin', 'val']]
df3 = df3.rename(columns={ "val": "val_%s"%next_tag,\
#                          "error": "err_%s"%next_tag\
                         })

df4 = df3.merge(df2, how='left',on=['parname','year', 'bin'], suffixes=('', '_remove'))
df4['difference'] = df4['val_%s'%this_tag] -  df4['val_%s'%next_tag]
print (df)
exit(0)
pdb.set_trace()

## now examine bin2 

df['difference'] = np.nan
for ipar in parlist:
  for ibin in bin_range:
    mask = (df['bin'] == ibin) & (df['year'] == my_year) & (df['parname'] == ipar)
    df_mask = df[mask]
    df_mask['difference'] = df_mask['val'] - df_mask['val'].shift(-1)
    df.update(df_mask)


for ipar in parlist:
    mask = (df['year'] == my_year) & (df['parname'] == ipar) & (df['difference']>-99999999)
    df_mask = df[mask]
    print df_mask[['parname', 'bin', 'val', 'difference']]

pdb.set_trace()    
# for ipar in ['Yield']:
#     mask = (df['year'] == my_year) & (df['parname'] == ipar) 
#     df_mask = df[mask]
#     print df_mask[df_mask.bin==0]
#     print df_mask[['parname', 'bin', 'val']]
# 
# print 'now printing significances'
# for ibin in bin_range:
#   for itag in fnames.values():
#     ipar = 'Yield'
#     mask_s = (df['year'] == my_year) & (df['parname'] == ipar) & (df['tag'] == itag) & (df['bin'] == ibin)
#     df_s = df[mask_s]
#     jpar = 'nbkg'
#     mask_b = (df['year'] == my_year) & (df['parname'] == jpar) & (df['tag'] == itag) & (df['bin'] == ibin)
#     df_b = df[mask_b]
#     pdb.set_trace()    
    
#     print itag, ibin, df_s['val'].iloc[0] / math.sqrt(df_b['val'].iloc[0] + df_s['val'].iloc[0])
# df[df.parname==ipar]    

# pdb.set_trace()    
#     sara['diff'] = sara['val'].diff()
#     if firsttime==0:
#       df = df.merge(sara, how='left',on=['parname','year', 'bin'], suffixes=('', '_remove'))
#       df.drop([icol for icol in df.columns if 'remove' in icol], axis=1, inplace=True)
#       df.drop_duplicates(keep='first', subset=['val', 'error'] ,inplace=True)
#       firsttime = 1
#     df.update(sara)
# #     df = df.merge(sara,how='left',on=['parname','year', 'bin'], suffixes=('', '_remove'))
#     pdb.set_trace()
# #     df.drop([icol for icol in df.columns if 'remove' in icol], axis=1, inplace=True)
#   
# #   tmp_df = df.loc[(df['bin'] == 2) & (df['year'] == 2018) & (df['parname'] == ipar)]        
# #   tmp_df['val_dif'] = tmp_df['val'].diff()
# #   df = df.merge(tmp_df,how='left',on=['parname','year', 'bin'], suffixes=('', '_remove'))
# #   pdb.set_trace()
# # 
# # df.drop([icol for icol in df.columns if 'remove' in icol], axis=1, inplace=True)
# 
# #   df['val_dif'] = np.where((df['bin'] == 2) & (df['year'] == 2018) & (df['parname'] == ipar), df['val'] - df['val'].shift(), np.nan)
# 
# #   mask = (df['bin'] == 2) & (df['year'] == 2018) & (df['parname'] == ipar)
# #   print ipar
# #   df[mask]['val'], df[mask]['val'].diff()
# # 
# # #             thewin  = thedupl.loc[thedupl.groupby(['eventN', 'lumi', 'runN'])['bdt_prob'].idxmax()]
# # 
# # # prova groupby
# # df.loc[df.groupby(['tag'])['val']]
# # 
# # sara=df[mask]
# # sara['diff'] = sara['val'].diff()
# # df = df.merge(sara,how='left',on=['parname','year', 'bin'], suffixes=('', '_remove'))
# # df.drop([icol for icol in df.columns if 'remove' in icol], axis=1, inplace=True)
# #   df['new_dif'] = df[mask]['val'].diff()
#      
# #         dataset = dataset.merge(dataset_support,how='left',on=['bEta', 'bPt'], suffixes=('', '_remove'))
# 
#         
# #     this_fit = fit_pars()
# #         setattr(this_fit, 'val', ufloat(ifr.floatParsFinal().find( ipar).getVal(), ifr.floatParsFinal().find( ipar).getError()))
# #         setattr(this_fit, 'year', int(ifr.GetName().split('_')[-1]))
# #         setattr(this_fit, 'bin', int(ifr.GetName().split('_')[-2]))
# 
# #         par_dict[parList]
# # 
# #     re.findall(r'-?\d+', i)
#         
# #     year bin yield uncertainty
#     
# #     fr.append(f.Get('w')) 
# #     f.Close()
# #   except:
# #     print ('file %s or workspace not found'%(f))
# # 
# # 
# # # study mass
# # for ii, ibin in enumerate(bin_range):
# #     
# #     keys = fo.GetListOfKeys()
# #     pars_o = fit_pars_mass()
# #     for i in keys:
# #         if 'results_RT_%s'%ibin in i.GetName():   
# #             fitres = fo.Get(i.GetName()) 
# #             setattr(pars_o, 'mRT_mc', ufloat (fitres.floatParsFinal().find( 'mean^{RT%s}'%ibin).getVal(), 
# #                                               fitres.floatParsFinal().find( 'mean^{RT%s}'%ibin).getError()))
# #         if 'results_WT_%s'%ibin in i.GetName():   
# #             fitres = fo.Get(i.GetName()) 
# #             setattr(pars_o, 'mWT_mc', ufloat (fitres.floatParsFinal().find( 'mean^{WT%s}'%ibin).getVal(), 
# #                                               fitres.floatParsFinal().find( 'mean^{WT%s}'%ibin).getError()))
# # 
# #         if 'results_data_%s'%ibin in i.GetName():   
# #             fitres = fo.Get(i.GetName()) 
# #             setattr(pars_o, 'mRT', ufloat (fitres.floatParsFinal().find( 'mean^{RT%s}'%ibin).getVal(), 
# #                                            fitres.floatParsFinal().find( 'mean^{RT%s}'%ibin).getError()))
# #             setattr(pars_o, 'mWT', ufloat (fitres.floatParsFinal().find( 'mean^{WT%s}'%ibin).getVal(), 
# #                                            fitres.floatParsFinal().find( 'mean^{WT%s}'%ibin).getError()))
# #             
# # #             for ipar in par_list_mass.keys():
# # #                 setattr(pars_o, ipar, ufloat (fitres.floatParsFinal().find( ipar).getVal(), 
# # #                                               fitres.floatParsFinal().find( ipar).getError()))
# # # #             print 'bin: ', ibin
# # #             print pars_o
# #     
# # 
# #     print '\nbin: ', ibin
# #     for ipar in par_list_mass.keys():
# #         print ipar, '\t', getattr(pars_o, ipar)
# #     print 'deltaM MC: \t'  , pars_o.mRT_mc - pars_o.mWT_mc
# #     print 'deltaM data: \t', pars_o.mRT    - pars_o.mWT  
# # 
# 
# 
# # fo.Close()
# # fn.Close()
# 
# 
