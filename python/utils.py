import os
import glob
from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors, TH1
from Analysis.Tools.utils import Process

TH1.SetDefaultSumw2()

def data_info(name):
   # sample info from directories with name like mssmhbb_fh_2017.cfg
   name = name.split('.')[0]
   info = {}
   info['cfg']   = '_'.join(name.split('_')[:-2])
   info['year']  = name.split('_')[-1]
   info['mode']  = name.split('_')[-2]
   return info

      # this class has to prepare pt, eta, phi and possibly m12 distributions, flavour dependent?
class MssmHbb:
   def __init__(self,process,histo_dir='',flavours=['b','bb','c','cc','udsg'],variables=['pt','eta','phi']):
      self.m_proc = process
      self.m_histo_dir = histo_dir
      self.m_flvs = flavours
      self.m_vars = variables
      self.m_jets = []
      self.m_dijets = []
      self.m_histos = {}
      self.m_histos_dijet = {}
      self.m_histos_proc = {}
      self.m_histos_dijet_proc = {}
      
   def histogram_directory(self):
      return self.m_histo_dir
      
   def flavours(self):
      return self.m_flvs
      
   def histograms(self,process=False):
      if process:
         return self.m_histos_proc
      return self.m_histos
      
   def histograms_dijet(self,process=False):
      if process:
         return self.m_histos_dijet_proc
      return self.m_histos_dijet
      
   def jets(self):
      return self.m_jets
      
   def dijets(self):
      return self.m_dijets
      
   # read histograms from the files   
   def readTH1(self,lumi_scale=False):
      h1 = {}
      pbins = []
      jets = []
      njets = -1
      h_types = self.m_vars
      flv_all = self.flavours()
      flv_all.append('all')
      for pbin, pfile in self.m_proc.rootFiles().items():  # need to think better these loops to avoid more loops below
         ls = self.m_proc.luminosityScale()[pbin]
         print('anal_mssm: lumi scale: ',pbin,ls,self.m_proc.crossSections()[pbin],self.m_proc.alias())
         h1[pbin] = {}
         pbins.append(pbin)
         tfile = TFile( pfile, 'OLD' )
         # get number of jets
         h_dir = tfile.GetDirectory(self.m_histo_dir)
         h_dir_keys = h_dir.GetListOfKeys()
         histo_names = {}
         for typ in h_types:
            histo_names[typ] = [x.GetName() for x in h_dir_keys if typ+'_jet' in x.GetName() and x.GetName().count('_')==1 and len(x.GetName())==(len(typ)+5)]
         if not jets:  # number of jets is the same for all process bins
            njets = len(histo_names['pt'])
            for j in range(1,njets+1):
               jet = 'jet{}'.format(j)
               jets.append(jet)
         for flv in flv_all:
            h1[pbin][flv]={}
            for jet in jets:
               h1[pbin][flv][jet]={}
               for typ in h_types:
                  h_name = '{}_{}_{}'.format(typ,jet,flv)
                  if flv == 'all':
                     h_name = '{}_{}'.format(typ,jet)
                  h1[pbin][flv][jet][typ] = h_dir.Get(h_name)
                  h1[pbin][flv][jet][typ].SetName('{}_{}'.format(h_name,pbin))
                  h1[pbin][flv][jet][typ].SetDirectory(0) # Important!
                  if lumi_scale:
                     h1[pbin][flv][jet][typ].Scale(ls)
         tfile.Close()
         
      self.m_histos_proc = h1
      self.m_jets = jets
         
      h1_add_proc = {}
      for f,flv in enumerate(flv_all):
         h1_add_proc[flv] = {}
         for jet in jets:
            h1_add_proc[flv][jet] = {}
            for typ in h_types:
               h_name = typ+'_'+jet+'_'+flv
               for b,pbin in enumerate(pbins):
                  if b == 0:
                     h1_add_proc[flv][jet][typ] = h1[pbin][flv][jet][typ].Clone(h_name)
                  else:
                     h1_add_proc[flv][jet][typ].Add(h1[pbin][flv][jet][typ])
                  h1_add_proc[flv][jet][typ].SetDirectory(0) # Important!
               
               
      self.m_histos = h1_add_proc
      
   # read histograms from the files   
   def readTH1Dijet(self,lumi_scale=False):
      h1 = {}
      pbins = []
      jets = []
      njets = -1
      h_types = self.m_vars
      flv_all = self.flavours()
      flv_all = ['all']
      for pbin, pfile in self.m_proc.rootFiles().items():  # need to think better these loops to avoid more loops below
         ls = self.m_proc.luminosityScale()[pbin]
         h1[pbin] = {}
         pbins.append(pbin)
         tfile = TFile( pfile, 'OLD' )
         # get number of jets
         h_dir = tfile.GetDirectory(self.m_histo_dir)
         h_dir_keys = h_dir.GetListOfKeys()
         histo_names = {}
         for typ in h_types:
            histo_names[typ] = [x.GetName() for x in h_dir_keys if typ+'_jet' in x.GetName() and x.GetName().count('_')==1 and len(x.GetName())==(len(typ)+6)]
         if not jets:  # number of jets is the same for all process bins
            for j2 in histo_names['m']:
               jets.append(j2[2:])
            njets = len(jets)
         for flv in flv_all:
            h1[pbin][flv]={}
            for jet in jets:
               h1[pbin][flv][jet]={}
               for typ in h_types:
                  h_name = '{}_{}_{}'.format(typ,jet,flv)
                  if flv == 'all':
                     h_name = '{}_{}'.format(typ,jet)
                  h1[pbin][flv][jet][typ] = h_dir.Get(h_name)
                  h1[pbin][flv][jet][typ].SetName('{}_{}'.format(h_name,pbin))
                  h1[pbin][flv][jet][typ].SetDirectory(0) # Important!
                  if lumi_scale:
                     h1[pbin][flv][jet][typ].Scale(ls)
         tfile.Close()
         
      self.m_histos_dijet_proc = h1
      self.m_dijets = jets
         
      h1_add_proc = {}
      for f,flv in enumerate(flv_all):
         h1_add_proc[flv] = {}
         for jet in jets:
            h1_add_proc[flv][jet] = {}
            for typ in h_types:
               h_name = typ+'_'+jet+'_'+flv
               for b,pbin in enumerate(pbins):
                  if b == 0:
                     h1_add_proc[flv][jet][typ] = h1[pbin][flv][jet][typ].Clone(h_name)
                  else:
                     h1_add_proc[flv][jet][typ].Add(h1[pbin][flv][jet][typ])
                  h1_add_proc[flv][jet][typ].SetDirectory(0) # Important!
               
               
      self.m_histos_dijet = h1_add_proc
      
############################################################################      
