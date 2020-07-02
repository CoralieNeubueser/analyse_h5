import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
# load defined functions
from utils import *
from drawFunctions import *
from rootpy import stl

r.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

# retrieve energy bins for either hepd: True, hepp: False
energyBins, energies, energyMax = getEnergyBins(True, False)

# get L/alpha bins
l_bins, l_x_bins = getLbins()
numPbin, Pbins = getPitchBins()

#for enbin of second energyBin                                                                                                                                                                    
enbin = 1
energy_selected = 6.484055042266846

# select a day
day = 20180801

# no SAA
treecut = 'ev.field>25000'

hist1D_L = r.TH1D('hist1D_L_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist1D_L_'+str(day)+'_'+str(energies[enbin])+'MeV', l_bins, np.array(l_x_bins))
hist1D_alpha = r.TH1D('hist1D_alpha_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist1D_alpha_'+str(day)+'_'+str(energies[enbin])+'MeV', numPbin-1, np.array(Pbins,dtype=float))
hist2D = r.TH2D('hist2D_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist2D_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float))
hist2D_en = r.TH2D('hist2D_flux_en_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist2D_flux_en_'+str(day)+'_'+str(energies[enbin])+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float))
hist2D_loc = r.TH2D('hist2D_loc_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist2D_loc_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', 180, -180,180, 90,-90,90)
hist2D_time = r.TH2D('hist2D_time_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', 'hist2D_time_flux_'+str(day)+'_'+str(energies[enbin])+'MeV', (60*24), 0, 24, 100, 0, 0.1)

# read in txt files with averages
path = 'averages/hepd/'
file = open(path+str(day)+".txt", "r")
av_Lalpha = dict()

for line in file: 
    # print(line)
    columns = np.array(line.split())
    energy = columns[0]
    if energy == str(energy_selected):
        # filll dictionary from (L, alpha) -> (mean, rms)
        av_Lalpha.update( {(columns[1],columns[2]):(columns[4],columns[6])} )

# read tree
inRoot = r.TFile( filename , 'update' )
tree = inRoot.tree

for ev in tree:
    L = ev.L
    energy = ev.energy
    # select a day
    if day!=ev.day:
        continue
    # reject SAA
    if ev.field<25000:
        continue
    # select an energy bin
    for en in energy:
        if en!=energy_selected:
            continue
        
        for ia,alpha in enumerate(ev.alpha):
            # match L to L bin
            L_bin = hist1D_L.GetBinLowEdge( hist1D_L.FindBin( L ) )
            # match alpha to pitch bin
            alpha_bin = int(hist1D_alpha.GetBinLowEdge( hist1D_alpha.FindBin( alpha ) )) 
        
            # test if keys exist in dict
            if (str(L_bin), str(alpha_bin)) in av_Lalpha:
            
                average = float(av_Lalpha[(str(L_bin), str(alpha_bin))][0])
                rms = float(av_Lalpha[(str(L_bin), str(alpha_bin))][1])            
                # get daily average in L-alpha cell
                fiveSigma = average + 5*rms
                
                flux = getattr(tree,"flux_"+str(L_bin)+"_"+str(alpha_bin))
                if len(flux) > ia and flux[ia]>fiveSigma:

                    print("L-alpha bins : ", L_bin, alpha_bin)
                    print("average :      ",  average)
                    print("rms :          ",  rms)
                    print('Five sigma: ', fiveSigma)
                    print('Found flux: ',flux[ia])

                    hist2D.Fill(L, alpha, flux[ia])
                    hist2D_en.Fill(L, alpha)
                    hist2D_loc.Fill(ev.Long, ev.Lat, flux[ia])
                    hist2D_time.Fill(ev.time,flux[ia])

prep2D(hist2D, 'L value', '#alpha_eq [deg]', '#sum#Phi', False)
prep2D(hist2D_en, 'L value', '#alpha_eq [deg]', '#entries', False)
prep2D(hist2D_loc, 'Longitude', 'Latitude', '#sum#Phi', False)
prep2D(hist2D_time, 't [h]', '#Phi', '#entries', False)

inRoot.WriteObject(hist2D,"hist2D_highFlux_energyBin_1",'kOverwrite')
inRoot.WriteObject(hist2D_en,"hist2D_highFlux_entries_energyBin_1",'kOverwrite')
inRoot.WriteObject(hist2D_loc,"hist2D_highFlux_location_energyBin_1",'kOverwrite')
inRoot.WriteObject(hist2D_time,"hist2D_highFlux_time_energyBin_1",'kOverwrite')

inRoot.Close()
