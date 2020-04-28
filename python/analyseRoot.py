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
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

# read tree
inRoot = r.TFile( filename , 'update' )

l_bins = 50

## write 2d histograms                                                                                                                                                      
##entries
#hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,1,5,9,10,190)
#hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,1,5,9,10,190)
#hist2D_l_pitch_en_noSAA=r.TH2D("hist2D_l_pitch_en_noSAA","hist2D_l_pitch_en_noSAA",l_bins,1,5,9,10,190)
#hist2D_B_l_en=r.TH2D("hist2D_B_l_en","hist2D_B_l_en",100,17000,55000,l_bins,1,5)
#hist2D_B_l_en_noSAA=r.TH2D("hist2D_B_l_en_noSAA","hist2D_B_l_en_noSAA",100,17000,55000,l_bins,1,5)
#hist2D_energy_pitch_en=r.TH2D("hist2D_energy_pitch_en","hist2D_energy_pitch_en",12,0,70,9,10,190)
#hist2D_energy_pitch_en_noSAA=r.TH2D("hist2D_energy_pitch_en_noSAA","hist2D_energy_pitch_en_noSAA",12,0,70,9,10,190)
#
##flux
#hist2D_l_pitch_noSAA=r.TH2D("hist2D_l_pitch_noSAA","hist2D_l_pitch_avFlux_noSAA",l_bins,1,5,9,10,190)
#hist2D_l_pitch_rms_noSAA=r.TH2D("hist2D_l_pitch_rms_noSAA","hist2D_l_pitch_rms_noSAA",l_bins,1,5,9,10,190)
#hist2D_l_pitch_rms=r.TH2D("hist2D_l_pitch_rms","hist2D_l_pitch_rms",l_bins,1,5,9,10,190)
#hist2D_energy_pitch=r.TH2D("hist2D_energy_pitch","hist2D_energy_pitch",12,0,70,9,10,190)
#hist2D_energy_pitch_noSAA=r.TH2D("hist2D_energy_pitch_noSAA","hist2D_energy_pitch_noSAA",12,0,70,9,10,190)
#hist2D_energy_pitch_rms=r.TH2D("hist2D_energy_pitch_rms","hist2D_energy_pitch_rms",12,0,70,9,10,190)
#hist2D_energy_pitch_rms_noSAA=r.TH2D("hist2D_energy_pitch_rms_noSAA","hist2D_energy_pitch_rms_noSAA",12,0,70,9,10,190)
#hist2D_B_l=r.TH2D("hist2D_B_l","hist2D_B_l",100,17000,55000,l_bins,1,5)
#hist2D_B_l_noSAA=r.TH2D("hist2D_B_l_noSAA","hist2D_B_l_noSAA",100,17000,55000,l_bins,1,5)
#hist2D_B_l_rms=r.TH2D("hist2D_B_l_rms","hist2D_B_l_rms",100,17000,55000,l_bins,1,5)
#hist2D_B_l_rms_noSAA=r.TH2D("hist2D_B_l_rms_noSAA","hist2D_B_l_rms_noSAA",100,17000,55000,l_bins,1,5)

histList1D = []
histList1D_rate = []
energyBins = 12
energies = [4.0, 8.96811, 10.9642, 13.4045, 16.388, 20.0355, 24.4949,  29.9468, 36.6122, 44.7611, 54.7237,  66.9037]
energyMax = 70
# geometrical factors
ele_GF = [ 0.76, 188.26, 326.64, 339.65, 344.99, 331.83, 304.73, 263.56, 217.33, 169.48, 117.31, 71.45 ]

# pick color map
# colormap = plt.cm.nipy_spectral
colors = [1,   920,  632,  416,
          600,  400, 616,   432,  800,
          820,  840,  860,  880,  900]
#[colormap(i) for i in np.linspace(0, 0.9,energyBins)]

for e in range(energyBins):
      
      diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
      if not diren:
            newdir = inRoot.mkdir("energy_"+str(energies[e])+"MeV")
            newdir.cd()
      else:
            print(diren)
            diren.cd()

      histList1D.append( r.TH1D('hist_en_'+str(energies[e])+'MeV', 'hist_en_'+str(energies[e])+'MeV', 1000000,0,1000) )
      histList1D_rate.append( r.TH1D('hist_en_rate_'+str(energies[e])+'MeV', 'hist_en_rate_'+str(energies[e])+'MeV', 1000000,0,1000) )

      r.gROOT.SetBatch(True)
      inRoot.tree.Draw('flux_en['+str(e)+']>>hist_en_'+str(energies[e])+'MeV','field>25000','')
      inRoot.tree.Draw('flux_en['+str(e)+']*'+str(ele_GF[e])+'>>hist_en_rate_'+str(energies[e])+'MeV','field>25000','')

      histList1D[e].GetXaxis().SetTitle("flux [Hz/cm^{2}#upoint sr]")
      histList1D[e].SetLineColor(colors[e])
      histList1D[e].SetLineWidth(2)
      histList1D_rate[e].GetXaxis().SetTitle("rate [Hz]")
      histList1D_rate[e].SetLineColor(colors[e])
      histList1D_rate[e].SetLineWidth(2)

      diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
      histList1D[e].SetDirectory(0)
      histList1D[e].SetDirectory(diren)
      histList1D[e].Write('',r.TObject.kOverwrite)
      histList1D_rate[e].SetDirectory(0)
      histList1D_rate[e].SetDirectory(diren)
      histList1D_rate[e].Write('',r.TObject.kOverwrite)

      inRoot.cd()
e=0

# SAA cut:
inRoot.tree.Draw('pitch:L>>hist2D_l_pitch_en_noSAA('+str(l_bins)+',1,5,9,10,190)', 'field>25000','colz')
hist2D_l_pitch_en_noSAA = r.gDirectory.Get("hist2D_l_pitch_en_noSAA")
inRoot.tree.Draw('L:field>>hist2D_B_l_en_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_en_noSAA = r.gDirectory.Get("hist2D_B_l_en_noSAA")
inRoot.tree.Draw('pitch:energy>>hist2D_energy_pitch_en_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_en_noSAA = r.gDirectory.Get("hist2D_energy_pitch_en_noSAA")
inRoot.tree.Draw('pitch:L:flux_en.flux_en>>hist2D_l_pitch_noSAA('+str(l_bins)+',1,5,9,10,190)', 'field>25000','colz')
hist2D_l_pitch_noSAA = r.gDirectory.Get("hist2D_l_pitch_noSAA")
inRoot.tree.Draw('pitch:L:pow(flux_en.flux_en,2)>>hist2D_l_pitch_rms_noSAA('+str(l_bins)+',1,5,9,10,190)', 'field>25000','colz')
hist2D_l_pitch_rms_noSAA = r.gDirectory.Get("hist2D_l_pitch_rms_noSAA")
inRoot.tree.Draw('pitch:energy:flux_en.flux_en>>hist2D_energy_pitch_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_noSAA = r.gDirectory.Get("hist2D_energy_pitch_noSAA")
inRoot.tree.Draw('pitch:energy:pow(flux_en.flux_en,2)>>hist2D_energy_pitch_rms_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_rms_noSAA =  r.gDirectory.Get("hist2D_energy_pitch_rms_noSAA")
inRoot.tree.Draw('L:field:flux_en.flux_en>>hist2D_B_l_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_noSAA=r.gDirectory.Get("hist2D_B_l_noSAA")
inRoot.tree.Draw('L:field:pow(flux_en.flux_en,2)>>hist2D_B_l_rms_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_rms_noSAA=r.gDirectory.Get("hist2D_B_l_rms_noSAA")

inRoot.tree.Draw('pitch:L:flux_en.flux_en>>hist2D_l_pitch_flux('+str(l_bins)+',1,5,9,10,190)','','colz')
hist2D_l_pitch_flux = r.gDirectory.Get('hist2D_l_pitch_flux')
inRoot.tree.Draw('pitch:L>>hist2D_l_pitch_en('+str(l_bins)+',1,5,9,10,190)')
hist2D_l_pitch_en = r.gDirectory.Get("hist2D_l_pitch_en")
inRoot.tree.Draw('L:field>>hist2D_B_l_en(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l_en= r.gDirectory.Get("hist2D_B_l_en")
inRoot.tree.Draw('pitch:L:pow(flux_en.flux_en,2)>>hist2D_l_pitch_rms('+str(l_bins)+',1,5,9,10,190)','','colz')
hist2D_l_pitch_rms = r.gDirectory.Get("hist2D_l_pitch_rms")
inRoot.tree.Draw('pitch:energy>>hist2D_energy_pitch_en('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch_en =  r.gDirectory.Get("hist2D_energy_pitch_en")
inRoot.tree.Draw('pitch:energy:flux_en.flux_en>>hist2D_energy_pitch('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch = r.gDirectory.Get("hist2D_energy_pitch")
inRoot.tree.Draw('pitch:energy:pow(flux_en.flux_en,2)>>hist2D_energy_pitch_rms('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch_rms = r.gDirectory.Get("hist2D_energy_pitch_rms")
inRoot.tree.Draw('L:field:flux_en.flux_en>>hist2D_B_l(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l=r.gDirectory.Get("hist2D_B_l")
inRoot.tree.Draw('L:field:pow(flux_en.flux_en,2)>>hist2D_B_l_rms(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l_rms=r.gDirectory.Get("hist2D_B_l_rms")


print ("hist entries: ", hist2D_l_pitch_rms.GetEntries())

# scale and normalise
# for rms, devide by n and take sqrt of entries
hist2D_l_pitch_rms_norm            = takeSqrt(hist2D_l_pitch_rms, hist2D_l_pitch_en)
hist2D_energy_pitch_rms_norm       = takeSqrt(hist2D_energy_pitch_rms, hist2D_energy_pitch_en)
hist2D_B_l_rms_norm                = takeSqrt(hist2D_B_l_rms, hist2D_B_l_en)
hist2D_l_pitch_rms_noSAA_norm      = takeSqrt(hist2D_l_pitch_rms_noSAA, hist2D_l_pitch_en_noSAA)
hist2D_energy_pitch_rms_noSAA_norm = takeSqrt(hist2D_energy_pitch_rms_noSAA, hist2D_energy_pitch_en_noSAA)
hist2D_B_l_rms_noSAA_norm          = takeSqrt(hist2D_B_l_rms_noSAA, hist2D_B_l_en_noSAA)

# for the average, devide by n per cell
hist2D_l_pitch_norm            = divideHists(hist2D_l_pitch_flux, hist2D_l_pitch_en)
hist2D_energy_pitch_norm       = divideHists(hist2D_energy_pitch, hist2D_energy_pitch_en)
hist2D_B_l_norm                = divideHists(hist2D_B_l, hist2D_B_l_en)
hist2D_l_pitch_noSAA_norm      = divideHists(hist2D_l_pitch_noSAA, hist2D_l_pitch_en_noSAA)
hist2D_energy_pitch_noSAA_norm = divideHists(hist2D_energy_pitch_noSAA, hist2D_energy_pitch_en_noSAA)
hist2D_B_l_noSAA_norm          = divideHists(hist2D_B_l_noSAA, hist2D_B_l_en_noSAA)

prep2D(hist2D_l_pitch_flux, 'L value', 'pitch [deg]', '#LT flux #GT', False)
prep2D(hist2D_l_pitch_en, 'L value', 'pitch [deg]', '#entries', False)
prep2D(hist2D_l_pitch_en_noSAA, 'L value', 'pitch [deg]', '#entries', False)
prep2D(hist2D_B_l_en_noSAA, 'field [nT]', 'L value', '#entries', False)
prep2D(hist2D_l_pitch_noSAA, 'L value', 'pitch [deg]', '#LT flux #GT', False)
prep2D(hist2D_l_pitch_rms_noSAA, 'L value', 'pitch [deg]', 'rms(flux)', False)
prep2D(hist2D_energy_pitch_en, 'energy [MeV]', 'pitch [deg]', '#entries', False)
prep2D(hist2D_energy_pitch_en_noSAA, 'energy [MeV]', 'pitch [deg]', '#entries', False)
prep2D(hist2D_energy_pitch_noSAA, 'energy [MeV]', 'pitch [deg]', '#LT flux #GT', False)
prep2D(hist2D_energy_pitch_rms_noSAA, 'energy [MeV]', 'pitch [deg]', 'rms(flux)', False)
prep2D(hist2D_B_l_noSAA, 'field [nT]', 'L value', '#LT flux #GT', False)
prep2D(hist2D_B_l_rms_noSAA,  'field [nT]', 'L value', 'rms(flux)', False)

prep2D(hist2D_B_l_en, 'field [nT]', 'L value', '#entries', False)
prep2D(hist2D_l_pitch_rms, 'L value', 'pitch [deg]', 'rms(flux)', False)
prep2D(hist2D_energy_pitch, 'energy [MeV]', 'pitch [deg]', '#LT flux #GT', False)
prep2D(hist2D_energy_pitch_rms, 'energy [MeV]', 'pitch [deg]', 'rms(flux)', False)
prep2D(hist2D_B_l, 'field [nT]', 'L value', '#LT flux #GT', False)
prep2D(hist2D_B_l_rms,  'field [nT]', 'L value', 'rms(flux)', False)

inRoot.WriteObject(hist2D_l_pitch_en,"hist2D_l_pitch_en",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_norm,"hist2D_l_pitch_flux",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_rms_norm,"hist2D_l_pitch_rms",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_en_noSAA,"hist2D_l_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_noSAA_norm,"hist2D_l_pitch_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_rms_noSAA_norm,"hist2D_l_pitch_rms_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_en,"hist2D_energy_pitch_en",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_norm,"hist2D_energy_pitch",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_rms_norm,"hist2D_energy_pitch_rms",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_en_noSAA,"hist2D_energy_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_noSAA_norm,"hist2D_energy_pitch_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_rms_noSAA_norm,"hist2D_energy_pitch_rms_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_en,"hist2D_B_l_en",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_norm,"hist2D_B_l",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_rms_norm,"hist2D_B_l_rms",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_en_noSAA,"hist2D_B_l_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_noSAA_norm,"hist2D_B_l_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_rms_noSAA_norm,"hist2D_B_l_rms_noSAA",'kOverwrite')

inRoot.Close()
