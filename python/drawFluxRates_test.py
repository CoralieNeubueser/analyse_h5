import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
# load defined functions
from utils import *

r.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
#parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filenames = ['/home/LIMADOU/cneubueser/analyse_h5/root/L3_test/L3h5_orig/all.root',
             '/home/LIMADOU/cneubueser/analyse_h5/root/L3_test/L3h5_rate/all.root',
             '/home/LIMADOU/cneubueser/analyse_h5/root/L3_test/L3h5_05_95/all.root',
             '/home/LIMADOU/cneubueser/analyse_h5/root/L3_test/L3h5_rate_05_95/all.root'
       ] #args.inputFile

l_bins = 50

## write 2d histograms
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
          600,  400, 616,   432,  
          800, 820,  840,  880,  
          860,  900]


for filename in filenames:
      # read tree
      inRoot = r.TFile( filename , 'update' )
      
      rates = ['hist_en_' ,'hist_en_log10_', 'hist_en_rate_']
      rebins = [2, 2, 1000]
      rangesMaxX = [15, 3, 15]
      rangesMinX = [0, -3, 0]
      if filename.find('rate')>0:
            print(filename)
            rebins = [1000, 2, 100]


      r.gStyle.SetPadRightMargin(0.1)
      i=0
      for rate in rates:
            
            can=r.TCanvas('can','can')
            can.SetLogy()
            legend=r.TLegend(0.7,0.5,0.9,0.9)

            for e in range(energyBins):
                  
                  diren = inRoot.GetDirectory("energy_"+str(energies[e])+"MeV")
                  if not diren:
                        print('Energy directory not found!')
                        break
                  else:
                        print(diren)
                        diren.cd()

                  readThis = r.TH1D('readThis','',100,0,100)
                  inRoot.GetObject('energy_'+str(energies[e])+'MeV/'+str(rate)+str(energies[e])+'MeV;1', readThis)
                  readThis.Rebin(rebins[i])
                  readThis.SetLineColor(colors[e])
                  readThis.SetMarkerColor(colors[e])
                  readThis.SetMinimum(1)
                  readThis.SetMaximum(1000*rebins[i])

                  if str(filename).find('_rate')>0:
                        readThis.GetXaxis().SetTitle('rate [Hz]')
                        readThis.GetXaxis().SetRangeUser(0,50)

                  readThis.GetXaxis().SetRangeUser(rangesMinX[i],rangesMaxX[i])

                  legend.AddEntry(readThis, str(energies[e])+"MeV")
                  if e==0:
                        readThis.Draw("hist")
                  else:
                        readThis.Draw("hist same")
            
                  r.gDirectory.cd()

            i+=1
            legend.Draw()
            outname=filename.replace("root","pdf").replace("all",rate+'all')
            print(outname)
            can.Print(outname)

      inRoot.Close()
