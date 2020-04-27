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
energyBins = 12
energies = [4.0, 8.96811, 10.9642, 13.4045, 16.388, 20.0355, 24.4949,  29.9468, 36.6122, 44.7611, 54.7237,  66.9037]
energyMax = 70
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

      histList1D.append( r.TH1D('hist_en_'+str(e), 'hist_en_'+str(e), 1000000,0,1000) )
      r.gROOT.SetBatch(True)
      inRoot.tree.Draw('flux_en['+str(e)+']>>hist_en_'+str(e),'field>25000','')

      histList1D[e].GetXaxis().SetTitle("flux [Hz/cm^{2}#upoint sr]")
      histList1D[e].SetLineColor(colors[e])
      histList1D[e].SetLineWidth(2)
      
      diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
      print("Write histo in: ", diren)
      histList1D[e].SetDirectory(0)
      histList1D[e].SetDirectory(diren)
      #      inRoot.WriteObject(histList1D[e], 'hist1D_flux_en_'+str(e),'kOverwrite')
      histList1D[e].Write('kOverwrite')

      inRoot.cd()
e=0

# SAA cut:
inRoot.tree.Draw('pitch:L>>hist2D_l_pitch_en_noSAA('+str(l_bins)+',1,5,9,10,190)', 'field>25000','colz')
hist2D_l_pitch_en_noSAA = r.gDirectory.Get("hist2D_l_pitch_en_noSAA")
inRoot.tree.Draw('L:field>>hist2D_B_l_en_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_en_noSAA = r.gDirectory.Get("hist2D_B_l_en_noSAA")
inRoot.tree.Draw('pitch:energy>>hist2D_energy_pitch_en_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_en_noSAA = r.gDirectory.Get("hist2D_energy_pitch_en_noSAA")
inRoot.tree.Draw('pitch:L:flux_en['+str(e)+']>>hist2D_l_pitch_noSAA('+str(l_bins)+',1,5,9,10,190,100,0,1)', 'field>25000','colz')
hist2D_l_pitch_noSAA = r.gDirectory.Get("hist2D_l_pitch_noSAA")
inRoot.tree.Draw('pitch:L:pow(flux_en['+str(e)+'],2)>>hist2D_l_pitch_rms_noSAA('+str(l_bins)+',1,5,9,10,190)', 'field>25000','colz')
hist2D_l_pitch_rms_noSAA = r.gDirectory.Get("hist2D_l_pitch_rms_noSAA")
inRoot.tree.Draw('pitch:energy:flux_en['+str(e)+']>>hist2D_energy_pitch_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_noSAA = r.gDirectory.Get("hist2D_energy_pitch_noSAA")
inRoot.tree.Draw('pitch:energy:pow(flux_en['+str(e)+'],2)>>hist2D_energy_pitch_rms_noSAA('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)', 'field>25000','colz')
hist2D_energy_pitch_rms_noSAA =  r.gDirectory.Get("hist2D_energy_pitch_rms_noSAA")
inRoot.tree.Draw('L:field:flux_en['+str(e)+']>>hist2D_B_l_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_noSAA=r.gDirectory.Get("hist2D_B_l_noSAA")
inRoot.tree.Draw('L:field:pow(flux_en['+str(e)+'],2)>>hist2D_B_l_rms_noSAA(100,17000,55000,'+str(l_bins)+',1,5)', 'field>25000','colz')
hist2D_B_l_rms_noSAA=r.gDirectory.Get("hist2D_B_l_rms_noSAA")
print("No error until here.")
inRoot.tree.Draw('pitch:L:flux_en['+str(e)+']>>hist2D_l_pitch('+str(l_bins)+',1,5,9,10,190,100,0,1)','','colz')
hist2D_l_pitch = r.gDirectory.Get('hist2D_l_pitch')
inRoot.tree.Draw('pitch:L>>hist2D_l_pitch_en('+str(l_bins)+',1,5,9,10,190)')
hist2D_l_pitch_en = r.gDirectory.Get("hist2D_l_pitch_en")
inRoot.tree.Draw('L:field>>hist2D_B_l_en(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l_en= r.gDirectory.Get("hist2D_B_l_en")
inRoot.tree.Draw('pitch:L:pow(flux_en['+str(e)+'],2)>>hist2D_l_pitch_rms('+str(l_bins)+',1,5,9,10,190)','','colz')
hist2D_l_pitch_rms = r.gDirectory.Get("hist2D_l_pitch_rms")
inRoot.tree.Draw('pitch:energy>>hist2D_energy_pitch_en('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch_en =  r.gDirectory.Get("hist2D_energy_pitch_en")
inRoot.tree.Draw('pitch:energy:flux_en['+str(e)+']>>hist2D_energy_pitch('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch = r.gDirectory.Get("hist2D_energy_pitch")
inRoot.tree.Draw('pitch:energy:pow(flux_en['+str(e)+'],2)>>hist2D_energy_pitch_rms('+str(energyBins)+',0,'+str(energyMax)+',9,10,190)','','colz')
hist2D_energy_pitch_rms = r.gDirectory.Get("hist2D_energy_pitch_rms")
inRoot.tree.Draw('L:field:flux_en['+str(e)+']>>hist2D_B_l(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l=r.gDirectory.Get("hist2D_B_l")
inRoot.tree.Draw('L:field:pow(flux_en['+str(e)+'],2)>>hist2D_B_l_rms(100,17000,55000,'+str(l_bins)+',1,5)','','colz')
hist2D_B_l_rms=r.gDirectory.Get("hist2D_B_l_rms")

# scale and normalise
# for rms, devide by n and take sqrt of entries
#takeSqrt(hist2D_l_pitch_rms, hist2D_l_pitch_en)
#takeSqrt(hist2D_energy_pitch_rms, hist2D_energy_pitch_en)
#takeSqrt(hist2D_B_l_rms, hist2D_B_l_en)
#takeSqrt(hist2D_l_pitch_rms_noSAA, hist2D_l_pitch_en_noSAA)
#takeSqrt(hist2D_energy_pitch_rms_noSAA, hist2D_energy_pitch_en_noSAA)
#takeSqrt(hist2D_B_l_rms_noSAA, hist2D_B_l_en_noSAA)

# for the average, devide by n per cell
#hist2D_l_pitch.Divide(hist2D_l_pitch_en)
#hist2D_energy_pitch.Divide(hist2D_energy_pitch_en)
#hist2D_B_l.Divide(hist2D_B_l_en)
#hist2D_l_pitch_noSAA.Divide(hist2D_l_pitch_en_noSAA)
#hist2D_energy_pitch_noSAA.Divide(hist2D_energy_pitch_en_noSAA)
#hist2D_B_l_noSAA.Divide(hist2D_B_l_en_noSAA)

#prep2D(hist2D_l_pitch, 'L value', 'pitch [deg]', '#LT flux #GT', False)
#prep2D(hist2D_l_pitch_en, 'L value', 'pitch [deg]', '#entries', False)
#prep2D(hist2D_l_pitch_en_noSAA, 'L value', 'pitch [deg]', '#entries', False)
#prep2D(hist2D_B_l_en_noSAA, 'field [nT]', 'L value', '#entries', False)
#prep2D(hist2D_l_pitch_noSAA, 'L value', 'pitch [deg]', '#LT flux #GT', False)
#prep2D(hist2D_l_pitch_rms_noSAA, 'L value', 'pitch [deg]', 'rms(flux)', False)
#prep2D(hist2D_energy_pitch_en, 'energy [MeV]', 'pitch [deg]', '#entries', False)
#prep2D(hist2D_energy_pitch_en_noSAA, 'energy [MeV]', 'pitch [deg]', '#entries', False)
#prep2D(hist2D_energy_pitch_noSAA, 'energy [MeV]', 'pitch [deg]', '#LT flux #GT', False)
#prep2D(hist2D_energy_pitch_rms_noSAA, 'energy [MeV]', 'pitch [deg]', 'rms(flux)', False)
#prep2D(hist2D_B_l_noSAA, 'field [nT]', 'L value', '#LT flux #GT', False)
#prep2D(hist2D_B_l_rms_noSAA,  'field [nT]', 'L value', 'rms(flux)', False)
#
#prep2D(hist2D_B_l_en, 'field [nT]', 'L value', '#entries', False)
#prep2D(hist2D_l_pitch_rms, 'L value', 'pitch [deg]', 'rms(flux)', False)
#prep2D(hist2D_energy_pitch, 'energy [MeV]', 'pitch [deg]', '#LT flux #GT', False)
#prep2D(hist2D_energy_pitch_rms, 'energy [MeV]', 'pitch [deg]', 'rms(flux)', False)
#prep2D(hist2D_B_l, 'field [nT]', 'L value', '#LT flux #GT', False)
#prep2D(hist2D_B_l_rms,  'field [nT]', 'L value', 'rms(flux)', False)


inRoot.WriteObject(hist2D_l_pitch_en,"hist2D_l_pitch_en",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch,"hist2D_l_pitch",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_rms,"hist2D_l_pitch_rms",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_en_noSAA,"hist2D_l_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_noSAA,"hist2D_l_pitch_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_l_pitch_rms_noSAA,"hist2D_l_pitch_rms_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_en,"hist2D_energy_pitch_en",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch,"hist2D_energy_pitch",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_rms,"hist2D_energy_pitch_rms",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_en_noSAA,"hist2D_energy_pitch_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_noSAA,"hist2D_energy_pitch_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_energy_pitch_rms_noSAA,"hist2D_energy_pitch_rms_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_en,"hist2D_B_l_en",'kOverwrite')
inRoot.WriteObject(hist2D_B_l,"hist2D_B_l",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_rms,"hist2D_B_l_rms",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_en_noSAA,"hist2D_B_l_en_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_noSAA,"hist2D_B_l_noSAA",'kOverwrite')
inRoot.WriteObject(hist2D_B_l_rms_noSAA,"hist2D_B_l_rms_noSAA",'kOverwrite')


##print(colors)
#stack = r.THStack('stack','stack')
#for ene in range(energyBins):
#      histList1D[ene].GetXaxis().SetTitle("flux [Hz/cm^{2}#upoint sr]")
#      histList1D[ene].SetLineColor(colors[ene])
#      histList1D[ene].SetLineWidth(2)
#      stack.Add(histList1D[ene],"hist")
#      inRoot.WriteObject(histList1D[ene], 'hist1D_flux_en_'+str(ene),'kOverwrite')
#
#r.gROOT.SetBatch(False)
#canstack = r.TCanvas("can_stack","can_stack")
#stack.Draw("stack hist")
#r.gPad.Modified()
#r.gPad.Update()
#stack.GetXaxis().SetTitle("flux [Hz/cm^{2} sr]")
#stack.GetXaxis().SetRangeUser(0,1)
#stack.Draw("stack hist")
#r.gPad.Modified()
#r.gPad.Update()
#
#inRoot.WriteObject(canstack, 'flux_stack', 'kOverwrite')

inRoot.Close()
