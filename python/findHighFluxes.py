import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
# load defined functions
from utils import *
from drawFunctions import *

r.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.', required=True)
parser.add_argument('--data', type=str, help='Define input data.', required=True)
parser.add_argument('--thr', type=int, default=100, help='Define minimum statistics used.')
parser.add_argument('--sigma', type=int, default=1, help='Define minimum sigma for flux values > <phi>+sigma*phi_rms.')
parser.add_argument('--fitted', action='store_true', help='Use exponential fit tau value for threshold.')
parser.add_argument('--useVersion', type=str, default='v2', help='Set the version.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

# get L/alpha bins
l_bins, l_x_bins = getLbins()
numPbin, Pbins = getPitchBins()
storedEn = []
# retrieve energy bins for either hepd: True, hepp: False
if args.data=='hepd':
    en_bins, energies, en_max = getEnergyBins(True, False)
    storedEn.append(energies)
else:
    en_bins, energies, en_max = getEnergyBins(False, True)
    for ien,en in enumerate(energies):
        storedEn.append(round(en,1))

hist1D_L = r.TH1D('hist1D_L', 'hist1D_L', l_bins, np.array(l_x_bins))
hist1D_alpha =  r.TH1D('hist1D_alpha', 'hist1D_alpha', numPbin-1, np.array(Pbins,dtype=float))
hist2D = []
hist2D_en = []
hist2D_loc = []
hist2D_time = []
av_Lalpha = []

# read tree
inRoot = r.TFile( filename , 'update' )
tree = inRoot.tree

# check if energies in list                                                                                                        
test_energies = set()
for ev in inRoot.tree:
      if len(test_energies) < en_bins:
            for e in ev.energy:
                  test_energies.add(e)
      else:
            break
test_energies = sorted(test_energies, key=float)
if energies!=test_energies:
    energies=test_energies
print("Energy bins: ", energies)

# output tree
outfilename = filename.replace('all','all_highFluxes')
if args.fitted:
    outfilename = filename.replace('all','all_highFluxes_fittedExp')
outRoot = r.TFile( outfilename, 'recreate' )
out_tree = r.TTree( 'events', 'tree of fluxes' )

Ev = array('i', [0])
Flux = array( 'f', [0.] )
Signal = array( 'f', [0.] )
Day = array('i', [0])
Time = array( 'f', [0.] )
Longitude = array('i', [0])
Latitude = array('i', [0])
Longitude_geom = array('i', [0])
Latitude_geom = array('i', [0])
Lshell = array('f', [0.])
Alpha_eq = array('f', [0.])
Energy = array('f', [0.])
GeomInd = array('i', [0])

out_tree.Branch( 'event', Ev, 'event/I' )
out_tree.Branch( 'flux', Flux, 'flux/F' )
out_tree.Branch( 'counts', Signal, 'counts/F' )
out_tree.Branch( 'day', Day, 'day/I' )
out_tree.Branch( 'time', Time, 'time/F' )
out_tree.Branch( 'lon', Longitude, 'lon/I' )
out_tree.Branch( 'lat', Latitude, 'lat/I' )
out_tree.Branch( 'geom_lon', Longitude_geom, 'geom_lon/I' )
out_tree.Branch( 'geom_lat', Latitude_geom, 'geom_lat/I' )
out_tree.Branch( 'L', Lshell, 'L/F' )
out_tree.Branch( 'alpha', Alpha_eq, 'alpha/F' )
out_tree.Branch( 'energy', Energy, 'energy/F' )
out_tree.Branch( 'geomIndex', GeomInd, 'geomIndex/I' )

# get a list of all days, for which data was taken
days = getDays(tree)
count = int(0)
# determine max counts for histogram later
maxCounts=0

# read Data file with geom Indices
dataDict = readGeomIndex()

for day in days:
    
    # prepare histograms
    hist1D_alpha = r.TH1D('hist1D_alpha','hist1D_alpha', numPbin-1, np.array(Pbins,dtype=float))
    hist1D_L = r.TH1D('hist1D_L','hist1D_L', l_bins, np.array(l_x_bins))
    for en in energies:
        hist2D.append( r.TH2D('hist2D_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_flux_'+str(day)+'_'+str(en)+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float)) )
        hist2D_en.append( r.TH2D('hist2D_flux_en_'+str(day)+'_'+str(en)+'MeV', 'hist2D_flux_en_'+str(day)+'_'+str(en)+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float)) )
        hist2D_loc.append( r.TH2D('hist2D_loc_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_loc_flux_'+str(day)+'_'+str(en)+'MeV', 180, -180,180, 90,-90,90) )
        hist2D_time.append( r.TH2D('hist2D_time_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_time_flux_'+str(day)+'_'+str(en)+'MeV', (60*24), 0, 24, 100, 0, 0.1) )
        av_Lalpha.append( dict() )

    # read in txt files with averages
    path = sharedOutPath()+'/data/averages/'+args.useVersion+'/'+args.data+'/'
    if args.integral:
        path += str(args.integral)+'s/'
    if args.fitted:
        path += 'fittedExp/'

    print("Average/RMS read from file: ", path+str(day)+'_min_'+str(args.thr)+'ev.txt')
    file = open(path+str(day)+'_min_'+str(args.thr)+'ev.txt', "r")
    next(file)
    for line in file:
        columns = [float(i) for i in line.split()]
        # print(columns)
        col_energy = columns[0]
        energyStored = col_energy
        en_index = storedEn.index( energyStored )
        # filll dictionary from (L, alpha) -> (mean, rms)
        av_Lalpha[en_index].update( {(columns[1],int(columns[2])):(columns[4],columns[6],columns[7])} )
        
    for ev in tree:
        L = ev.L
        energy = ev.energy
        
        # select a day
        if day!=ev.day:
            continue
        # reject SAA
        if ev.field<25000:
            continue
    
        for ia,alpha in enumerate(ev.alpha):

            # match L to L bin
            L_bin = hist1D_L.GetBinLowEdge( hist1D_L.FindBin( L ) )
            # match alpha to pitch bin
            alpha_bin = int(hist1D_alpha.GetBinLowEdge( hist1D_alpha.FindBin( alpha ) )) 
            # match energy to energy bin
            energyStored = energy[ia]
            energy_bin = test_energies.index( energyStored )

            # test if keys exist in dict
            if (L_bin, alpha_bin) in av_Lalpha[energy_bin]:
            
                average = av_Lalpha[energy_bin][(L_bin, alpha_bin)][0]
                rms = av_Lalpha[energy_bin][(L_bin, alpha_bin)][1]
                rmsErr = av_Lalpha[energy_bin][(L_bin, alpha_bin)][2]
                # get daily average in L-alpha cell
                # set RMS99 as threshold for 1% highest fluxes
                xSigma = args.sigma*rms
                if rms==0:
                    # use bin width as threshold 
                    xSigma = args.sigma*(2.*rmsErr)
                # set sigma like threshold
                if args.fitted:
                    xSigma = average + args.sigma*rms 

                flux = getattr(tree,"flux_"+str(L_bin)+"_"+str(alpha_bin))
                
                if len(flux) > ia and flux[ia]>xSigma:
                    if args.debug:
                        print("L-alpha bins : ", L_bin, alpha_bin)
                        print("average :      ", average)
                        print("rms :          ", rms)
                        print('X sigma:       ', xSigma)
                        print('Found flux:    ', flux[ia])
                    
                    hist2D[energy_bin].Fill(L, alpha, flux[ia])
                    hist2D_en[energy_bin].Fill(L, alpha)
                    hist2D_loc[energy_bin].Fill(ev.Long, ev.Lat, flux[ia])
                    hist2D_time[energy_bin].Fill(ev.time,flux[ia])

                    # get minutes from digits in 'time'
                    storedTime = ev.time
                    hour = int(math.floor(storedTime))
                    minute = int(math.floor((storedTime - hour)*60))
                    geoIndex = dataDict[(day, hour, minute)]

                    # fill output tree
                    Ev[0] = count
                    Flux[0] = flux[ia]
                    if rmsErr!=0:
                        # define signal as #counts, rmsErr is half-width of flux distributions
                        Signal[0] = flux[ia]/(2.*rmsErr)
                    else:
                        Signal[0] = 1.
                    Day[0] = day
                    Time[0] = storedTime # daily hours
                    Longitude[0] = ev.Long
                    Latitude[0] = ev.Lat
                    Longitude_geom[0] = ev.geomLong
                    Latitude_geom[0] = ev.geomLat
      
                    # write in L bins
                    binned_L = hist1D_L.GetBinLowEdge(hist1D_L.FindBin(L)) 
                    Lshell[0] = binned_L
                    # write in alpha bins
                    binned_alpha = hist1D_alpha.GetBinLowEdge(hist1D_alpha.FindBin(alpha))
                    Alpha_eq[0] = binned_alpha
                    Energy[0] = energy[ia]
                    GeomInd[0] = geoIndex

                    out_tree.Fill()

                    if Signal[0]>maxCounts:
                        maxCounts=int(Signal[0])
                    count+=1

    for ie in range(en_bins):
        prep2D(hist2D[ie], 'L value', '#alpha_eq [deg]', '#sum#Phi', False)
        prep2D(hist2D_en[ie], 'L value', '#alpha_eq [deg]', '#entries', False)
        prep2D(hist2D_loc[ie], 'Longitude', 'Latitude', '#sum#Phi', False)
        prep2D(hist2D_time[ie], 't [h]', '#Phi', '#entries', False)
    
        inRoot.WriteObject(hist2D[ie],"hist2D_"+str(day)+"_highFlux_energyBin_"+str(round(energies[ie],1)),'kOverwrite')
        inRoot.WriteObject(hist2D_en[ie],"hist2D_"+str(day)+"_highFlux_entries_energyBin_"+str(round(energies[ie],1)),'kOverwrite')
        inRoot.WriteObject(hist2D_loc[ie],"hist2D_"+str(day)+"_highFlux_location_energyBin_"+str(round(energies[ie],1)),'kOverwrite')
        inRoot.WriteObject(hist2D_time[ie],"hist2D_"+str(day)+"_highFlux_time_energyBin_"+str(round(energies[ie],1)),'kOverwrite')
        
    # clean up
    inRoot.cd()
    hist2D.clear()
    hist2D_en.clear()
    hist2D_loc.clear()
    hist2D_time.clear()
    av_Lalpha.clear()

# write out
inRoot.Close()
outRoot.Write()
outRoot.Close()

# add second analysis part:
# 1. fit expenential to 'counts' per day
# 2. use the tau value as measure of the width
# 3. add 'significance' as counts/tau to tree in order to allow for selection
f = r.TFile( outfilename, 'update' )
t3 = f.Get( 'events')
print("Run in total over all: ",t3.GetEntries())
print("with maximum counts of ",maxCounts)

Significance = array('f', [0.])
newBranch = t3.Branch('significance', Significance, 'significance/F')
taus = []

for day in days:
    # loop through high flux tree
    # add histogram per day of counts with maximum counts
    hist1D_counts = r.TH1D('hist1D_counts_'+str(day),'hist1D_counts_'+str(day), int(maxCounts), 0, maxCounts)
    for h in t3:
        if h.day==day:
            hist1D_counts.Fill(h.counts)
    # fit exponential
    f1 = r.TF1("f1_"+str(day),"[0]*exp(-x/[1])",hist1D_counts.GetBinCenter(hist1D_counts.GetMaximumBin()), hist1D_counts.GetBinCenter(hist1D_counts.FindLastBinAbove(10)))
    f1.SetParameters(100,hist1D_counts.GetRMS())
    fresults = hist1D_counts.Fit(f1,"R")
    tau = f1.GetParameter(1)
    f.WriteObject(hist1D_counts,"hist1D_counts_"+str(day),'kOverwrite')
    taus.append(tau)

# use tau to determine significance
for h in t3:
    iday=list(days).index(h.day)
    Significance[0] = h.counts/taus[iday]
    newBranch.Fill()

t3.Write("events", r.TObject.kOverwrite) # save only the new version of the tree
f.Close()
