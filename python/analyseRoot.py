import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
# load defined functions
from utils import *

r.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--removePoles', action='store_true', help='Remove of polar regions, by cut on latidute <50.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

# read tree
inRoot = r.TFile( filename , 'update' )

histList1D = []
histList1D_ln = []
histList1D_rate = []
hist2D = []
hist2D_en = []

energyBins = 12
energies = [2.0, 6.5, 9.9, 12.2, 14.9, 18.2, 22.3, 27.2, 33.3, 40.7, 49.7, 60.8]
energyMax = 70

# geometrical factors
ele_GF = [ 0.76, 188.26, 326.64, 339.65, 344.99, 331.83, 304.73, 263.56, 217.33, 169.48, 117.31, 71.45 ]
ele_corr_GF = [ 131.128, 545.639, 560.297, 530.937, 477.827, 413.133, 334.176, 252.3, 204.52, 103.216, 77.5552, 61.1536 ]

# L bins
l_x_bins = []
for x in range(0,5):
      l_x_bins.append(1+0.2*x)
for x in range(0,9):
      l_x_bins.append(2+x)
l_bins=len(l_x_bins)-1

# pick color map
colors = [1,   920,  632,  416,
          600,  400, 616,   432,  800,
          820,  840,  860,  880,  900]

# define cuts
cut = 'field>25000'
cutFlux = 'field>25000&&flux!=0'
treecut = 'ev.field>25000'
if args.removePoles:
      cutFlux+='&&abs(Lat)<50'

hist_tot = r.TH1D('hist_en_tot', 'hist_en_tot', 1000*200,0,1000)
if filename.find('rate'):
      hist_tot = r.TH1D('hist_en_tot', 'hist_en_tot', 1000*200,0,1000)
inRoot.tree.Draw('flux>>hist_en_tot',cutFlux,'')
if filename.find('rate'):
      hist_tot.GetXaxis().SetTitle("rate [Hz]")
else:
      hist_tot.GetXaxis().SetTitle("flux [Hz/cm^{2}#upoint sr]")
hist_tot.SetLineColor(1)
hist_tot.SetLineWidth(2)
hist_tot.SetDirectory(0)
hist_tot.Write('',r.TObject.kOverwrite)
hist_tot_log = r.TH1D('hist_en_tot_log', 'hist_en_tot_log', 40,-3,3)
inRoot.tree.Draw('log10(flux)>>hist_en_tot_log',cutFlux,'')
if filename.find('rate'):
      hist_tot_log.GetXaxis().SetTitle("log10(rate) [Hz]")
else:
      hist_tot_log.GetXaxis().SetTitle("log10(rate) [Hz]")
hist_tot_log.SetLineColor(1)
hist_tot_log.SetLineWidth(2)
hist_tot_log.SetDirectory(0)
hist_tot_log.Write('',r.TObject.kOverwrite)

for e in range(energyBins):
      
      diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
      if not diren:
            newdir = inRoot.mkdir("energy_"+str(energies[e])+"MeV")
            newdir.cd()
      else:
            diren.cd()

      # rate/flux 1D histograms
      histList1D.append( r.TH1D('hist_en_'+str(energies[e])+'MeV', 'hist_en_'+str(energies[e])+'MeV', 1000*200,0,1000) )
      histList1D_ln.append( r.TH1D('hist_en_log10_'+str(energies[e])+'MeV', 'hist_en_log10_'+str(energies[e])+'MeV', 40,-3,3) )
      histList1D_rate.append( r.TH1D('hist_en_rate_'+str(energies[e])+'MeV', 'hist_en_rate_'+str(energies[e])+'MeV', 1000*200,0,1000) )
      # L-pitch maps
      hist2D_en.append( r.TH2D('hist2D_en_'+str(energies[e])+'MeV', 'hist2D_en_'+str(energies[e])+'MeV', l_bins,np.array(l_x_bins), 9,0,180) ) 
      hist2D.append( r.TH2D('hist2D_flux_'+str(energies[e])+'MeV', 'hist2D_flux_'+str(energies[e])+'MeV', l_bins,np.array(l_x_bins), 9,0,180) )

      cutFluxEn = 'field>25000&&flux_en['+str(e)+']!=0'
      if args.removePoles:
            cutFluxEn += '&&abs(Lat)<50'

      r.gROOT.SetBatch(True)
      inRoot.tree.Draw('flux_en['+str(e)+']*'+str(ele_GF[e])+'/'+str(ele_corr_GF[e])+'>>hist_en_'+str(energies[e])+'MeV',cutFluxEn,'')
      inRoot.tree.Draw('log10(flux_en['+str(e)+']*'+str(ele_GF[e])+'/'+str(ele_corr_GF[e])+')>>hist_en_log10_'+str(energies[e])+'MeV',cutFluxEn,'')
      inRoot.tree.Draw('flux_en['+str(e)+']*'+str(ele_GF[e])+'>>hist_en_rate_'+str(energies[e])+'MeV',cutFluxEn,'')

      histList1D[e].GetXaxis().SetTitle("flux [Hz/cm^{2}#upoint sr]")
      histList1D[e].SetLineColor(colors[e])
      histList1D[e].SetLineWidth(2)
      histList1D_ln[e].GetXaxis().SetTitle("log10(flux) [Hz/cm^{2}#upoint sr]")
      histList1D_ln[e].SetLineColor(colors[e])
      histList1D_ln[e].SetLineWidth(2)
      histList1D_rate[e].GetXaxis().SetTitle("rate [Hz]")
      histList1D_rate[e].SetLineColor(colors[e])
      histList1D_rate[e].SetLineWidth(2)

      diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
      histList1D[e].SetDirectory(0)
      histList1D[e].SetDirectory(diren)
      histList1D[e].Write('',r.TObject.kOverwrite)
      histList1D_ln[e].SetDirectory(0)
      histList1D_ln[e].SetDirectory(diren)
      histList1D_ln[e].Write('',r.TObject.kOverwrite)
      histList1D_rate[e].SetDirectory(0)
      histList1D_rate[e].SetDirectory(diren)
      histList1D_rate[e].Write('',r.TObject.kOverwrite)

      inRoot.cd()
e=0


# loop through ev to fill maps
for ev in inRoot.tree:
      ev.L
      ev.pitch
      for enbin in range(energyBins):
            if ev.flux_en[enbin]!=0 and eval(treecut):
                  oldbin = hist2D[enbin].FindBin(ev.L,float(ev.pitch[0]))
                  oldCount = hist2D[enbin].GetBinContent(oldbin)
                  if oldCount != 0.:
                        # fill with flux using new geom factor
                        hist2D[enbin].SetBinContent(oldbin, oldCount + ev.flux_en[enbin]*ele_GF[enbin]/ele_corr_GF[enbin])
                        hist2D_en[enbin].SetBinContent(oldbin, hist2D_en[enbin].GetBinContent(hist2D_en[enbin].FindBin(ev.L,float(ev.pitch[0]))) + 1)
                  else:
                        hist2D[enbin].SetBinContent(oldbin, ev.flux_en[enbin]*ele_GF[enbin]/ele_corr_GF[enbin])
                        hist2D_en[enbin].SetBinContent(oldbin, 1)

for enbin in range(energyBins):
      diren = inRoot.GetDirectory("energy_"+str(energies[enbin])+"MeV")
      diren.cd()
      prep2D(hist2D[enbin], 'L value', 'pitch [deg]', '#sum flux', False)
      prep2D(hist2D_en[enbin], 'L value', 'pitch [deg]', '#entries', False)
      hist2D[enbin].SetDirectory(0)
      hist2D[enbin].SetDirectory(diren)
      hist2D[enbin].Write('',r.TObject.kOverwrite)
      hist2D_en[enbin].SetDirectory(0)
      hist2D_en[enbin].SetDirectory(diren)
      hist2D_en[enbin].Write('',r.TObject.kOverwrite)
      inRoot.cd()

# SAA cut:
inRoot.tree.Draw('pitch:L>>hist2D_l_pitch_en_noSAA('+str(l_bins)+','+str(np.array(l_x_bins))+',9,10,190)', cut,'colz')
hist2D_l_pitch_en_noSAA = r.gDirectory.Get("hist2D_l_pitch_en_noSAA")
inRoot.tree.Draw('L:field>>hist2D_B_l_en_noSAA(100,17000,55000,'+str(l_bins)+','+str(np.array(l_x_bins))+')', cut,'colz')
hist2D_B_l_en_noSAA = r.gDirectory.Get("hist2D_B_l_en_noSAA")
inRoot.tree.Draw('pitch:energy>>hist2D_energy_pitch_en_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', cut,'colz')
hist2D_energy_pitch_en_noSAA = r.gDirectory.Get("hist2D_energy_pitch_en_noSAA")

prep2D(hist2D_l_pitch_en_noSAA, 'L value', 'pitch [deg]', '#entries', False)
prep2D(hist2D_B_l_en_noSAA, 'field [nT]', 'L value', '#entries', False)
prep2D(hist2D_energy_pitch_en_noSAA, 'energy [MeV]', 'pitch [deg]', '#entries', False)

inRoot.WriteObject(hist2D_l_pitch_en_noSAA,"hist2D_l_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_en_noSAA,"hist2D_energy_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_en_noSAA,"hist2D_B_l_en_noSAA",'kOverwrite')

inRoot.Close()
