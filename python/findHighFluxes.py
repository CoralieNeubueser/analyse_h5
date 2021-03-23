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
parser.add_argument('--data', type=str, choices=['hepd','hepp_l','hepp_h'], help='Define input data.', required=True)
parser.add_argument('--thr', type=int, default=100, help='Define minimum statistics used.')
parser.add_argument('--sigma', type=int, default=1, help='Define minimum sigma for flux values > <phi>+sigma*phi_rms.')
parser.add_argument('--fitted', action='store_true', help='Use exponential fit tau value for threshold.')
parser.add_argument('--useVersion', type=str, default='v2', help='Set the version.')
parser.add_argument('--day', type=int, help='Define the day of the fluxes that are suppose to be stored.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

# get L/alpha bins
# go in DeltaL=1. steps to increase stats
l_bins, l_fine_bins = getLbins()
l_bins, l_x_bins = getCorserLbins()
numPbin, Pbins = getPitchBins()
storedEn = []
print(l_x_bins)
print(Pbins)
# retrieve energy bins for either hepd: True, hepp: False
if args.data=='hepd':
    en_bins, energies, en_max = getEnergyBins(args.data, False)
    storedEn = energies
else:
    en_bins, energies, en_max = getEnergyBins(args.data, True)
    storedEn = [round(en,1) for en in energies] 

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
#if energies!=test_energies:
#    energies=test_energies
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
Lshell = array('i', [0])
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
out_tree.Branch( 'L', Lshell, 'L/I' )
out_tree.Branch( 'alpha', Alpha_eq, 'alpha/F' )
out_tree.Branch( 'energy', Energy, 'energy/F' )
out_tree.Branch( 'geomIndex', GeomInd, 'geomIndex/I' )

# get a list of all days, for which data was taken
if not args.day:
    days = getDays(tree)
else:
    days = [args.day]
count = int(0)
# determine max counts for histogram later
maxCounts=0

# read Data file with geom Indices
dataDict = readGeomIndex()

for day in days:
    
    # prepare histograms
    hist2D = [ r.TH2D('hist2D_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_flux_'+str(day)+'_'+str(en)+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float)) for en in energies]
    hist2D_en = [ r.TH2D('hist2D_flux_en_'+str(day)+'_'+str(en)+'MeV', 'hist2D_flux_en_'+str(day)+'_'+str(en)+'MeV', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float)) for en in energies]
    hist2D_loc = [ r.TH2D('hist2D_loc_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_loc_flux_'+str(day)+'_'+str(en)+'MeV', 180, -180,180, 90,-90,90) for en in energies]
    hist2D_time = [ r.TH2D('hist2D_time_flux_'+str(day)+'_'+str(en)+'MeV', 'hist2D_time_flux_'+str(day)+'_'+str(en)+'MeV', (60*24), 0, 24, 100, 0, 0.1) for en in energies]
    av_Lalpha = [ {} for en in energies]

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
        # fill dictionary from (L, alpha) -> (mean, rms)
        av_Lalpha[en_index].update( {(columns[1],int(columns[2])):(columns[4],columns[6],columns[7])} )
        #av_Lalpha[en_index] = {(columns[1],int(columns[2])):(columns[4],columns[6],columns[7])}

    for ev in tree:
        L = ev.L
        energy = ev.energy
        
        # select a day
        if day!=ev.day:
            continue
        # reject SAA
        if ev.field<25000:
            continue
        # L-range
        if math.floor(L)>=10.:
            continue
        # run over max. 1000 events in debug mode
        if args.debug and ev.event>1000:
            print("Attention! In debug mode limited data set of 1000 events are used.")
            break

        # loop over every flux
        for ia,alpha in enumerate(ev.alpha):

            # match L to L bin
            if L<2:
                L_bin = math.floor( ((L-1.)/0.2) )
            else:
                L_bin = l_fine_bins.index(math.floor(L))
            L_binValue = l_fine_bins[L_bin]
            # match alpha to pitch bin
            alpha_bin = int(float(alpha)/20.) 
            alpha_binValue = alpha_bin*20
            # match energy to energy bin
            energyStored = energy[ia]
            energy_bin = test_energies.index( energyStored )

            # test if keys exist in dict of estimated RMS99/taus
            if (L_binValue, alpha_binValue) in av_Lalpha[energy_bin]:
            
                average = av_Lalpha[energy_bin][(L_binValue, alpha_binValue)][0]
                rms = av_Lalpha[energy_bin][(L_binValue, alpha_binValue)][1]
                rmsErr = av_Lalpha[energy_bin][(L_binValue, alpha_binValue)][2]
                # get daily average in L-alpha cell
                # set RMS99 as threshold for 1% highest fluxes
                xSigma = args.sigma*rms
                if rms==0:
                    # use bin width as threshold 
                    xSigma = args.sigma*(2.*rmsErr)
                # set sigma like threshold
                if args.fitted:
                    xSigma = average + args.sigma*rms 
                # get array of fluxes in the L-alpha bin
                flux = getattr(tree,"flux_"+str(l_fine_bins[L_bin])+"_"+str(Pbins[alpha_bin]))
                
                # if flux is above threshold
                if len(flux) > ia and flux[ia]>xSigma:
                    if args.debug:
                        print("L-alpha :      ", ev.L, alpha)
                        print("L-alpha bins : ", L_bin, alpha_bin)
                        print("average :      ", average)
                        print("rms :          ", rms)
                        print('X sigma:       ', xSigma)
                        print('Found flux:    ', flux[ia])
                    
                    hist2D[energy_bin].Fill(L, alpha, flux[ia])
                    hist2D_en[energy_bin].Fill(L, alpha)
                    hist2D_loc[energy_bin].Fill(ev.Long, ev.Lat, flux[ia])
                    hist2D_time[energy_bin].Fill(ev.time, flux[ia])

                    # get minutes from digits in 'time'
                    storedTime = ev.time
                    hour = int(math.floor(storedTime))
                    minute = int(math.floor((storedTime - hour)*60))
                    if (day, hour, minute) in dataDict:
                        geoIndex = dataDict[(day, hour, minute)]
                    else:
                        geoIndex = -1
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
                    binned_L = l_x_bins[l_x_bins.index(math.floor(L))]
                    Lshell[0] = int(binned_L)
                    # write in alpha bins
                    binned_alpha = Pbins[alpha_bin]
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
Tau = array('f', [0.])
Mean = array('f', [0.])
newBranch = t3.Branch('significance', Significance, 'significance/F')
newBranch2 = t3.Branch('tau', Tau, 'tau/F')
newBranch3 = t3.Branch('mean', Mean, 'mean/F')
taus = {}
means = {}
# loop through high flux tree
vecCells = defaultdict(list)
failedFits = 0
fit_summary = [ r.TH2D('hist2D_'+str(day)+'_fit_summary', 'hist2D_'+str(day)+'_fit_summary', l_bins, np.array(l_x_bins), numPbin-1, np.array(Pbins,dtype=float)) for day in days ] 
# fill high fluxes in dict that is used to fill histgram that are fit
for h in t3:
    alpha_index = Pbins.index(h.alpha)
    if h.alpha>80:
        alpha_index = Pbins.index(abs(h.alpha-160))
    vecCells[h.day, l_x_bins.index(math.floor(h.L)), alpha_index].append(h.counts)

for key,value in vecCells.items():
    # key[0]=day, key[1]=L, key[2]=alpha
    # value=[counts] to fit
    # add histogram per day of counts with maximum counts 
    hist1D_counts = r.TH1D('hist1D_counts_'+str(key[0])+'_'+str(l_x_bins[key[1]])+'_'+str(Pbins[key[2]]), 'hist1D_counts_'+str(key[0])+'_'+str(l_x_bins[key[1]])+'_'+str(Pbins[key[2]]), int(maxCounts), 0, maxCounts) 
    # fill histrogram
    for val in value:
        hist1D_counts.Fill(val)
    # if too little stats, continue
    if hist1D_counts.GetEntries()<20:
        continue

    print("START OF FITTING... ")
    print("................f1_"+str(key[0])+"_"+str(key[1])+"_"+str(key[2])+"... ")
    maxbin=hist1D_counts.GetMaximumBin()
    lastbin=hist1D_counts.FindLastBinAbove(1)
    # fit exponential
    f1 = r.TF1("f1_"+str(key[0])+'_'+str(key[1])+'_'+str(key[2]),"[0]*exp(-x/[1])",hist1D_counts.GetBinCenter(maxbin),hist1D_counts.GetBinCenter(lastbin))
    f1.SetParameters(hist1D_counts.GetEntries(),hist1D_counts.GetRMS())
    f1.SetParLimits(1, 1e-5,50)
    # first fitting in maximum range
    fresults = hist1D_counts.Fit(f1,"RSQ")
    mean = hist1D_counts.GetMean()
    tau = 1
    chi2 = 1000
    failed = True

    # if fit worked
    if fresults.Ndf()!=0 and fresults.IsValid():
        chi2 = fresults.Chi2()/fresults.Ndf()
        failed = False
        useFirst = True
        print("First tau: ", fresults.GetParams()[1])

    # if the fit didn't work or the chi2 too high, try again..
    if chi2>10 or failed:
        fitstart=maxbin
        trial = 0
        while trial<10:
            trialbin=trial
            # make sure only fits across 50% of entries is performed
            if hist1D_counts.Integral(trialbin,lastbin)/hist1D_counts.GetEntries()<0.5:
                print("Entries in fit range: ", hist1D_counts.Integral(trialbin,lastbin))
                print("Less than 50%, it will not be continued: ", hist1D_counts.Integral(trialbin,lastbin)/hist1D_counts.GetEntries())
                break
            fres = hist1D_counts.Fit(f1,"SQ","",hist1D_counts.GetBinCenter(trialbin),hist1D_counts.GetBinCenter(lastbin))
            # if fit converged
            if fres.IsValid() and fres.Ndf()!=0:
                # check chi2..
                if fres.Chi2()/fres.Ndf() < chi2:
                    chi2 = fres.Chi2()/fres.Ndf()
                    fitstart = trialbin
                    failed = False
                    useFirst = False
            trial+=1
        #        if args.debug:
        print('Finished after {} trials.'.format(trial))
        if not failed and not useFirst:
            # add last fit with optimised fit range
            fresults = hist1D_counts.Fit(f1,"SQ","",hist1D_counts.GetBinCenter(fitstart),lastbin)

    if not failed:
        tau = fresults.GetParams()[1]
        #        if args.debug:
        print("used tau: ",tau)
        print("fit results with chi2: ",fresults.Chi2()/fresults.Ndf())
    else:
        # if the fit did not converge with good chi2, remove from the histogram
        hist1D_counts.RecursiveRemove( hist1D_counts.FindObject("f1_"+str(key[0])+'_'+str(key[1])+'_'+str(key[2])) )

    f.WriteObject(hist1D_counts,"hist1D_counts_"+str(key[0])+'_'+str(l_x_bins[key[1]])+'_'+str(Pbins[key[2]]),'kOverwrite')
    taus[key]=tau
    means[key]=mean
    fit_summary[list(days).index(key[0])].Fill(l_x_bins[key[1]], Pbins[key[2]], tau)
    if key[2]!=4:
        fit_summary[list(days).index(key[0])].Fill(l_x_bins[key[1]], Pbins[8-key[2]], tau)

for iday,hist in enumerate(fit_summary):
    prep2D(hist, 'L shell', '#alpha_eq [deg]', '#tau', False)
    f.WriteObject(hist,"hist2D_"+str(list(days)[iday])+'_fit_summary','kOverwrite')
vecCells.clear()

# use tau to determine significance
for h in t3:
    iL=l_x_bins.index(math.floor(h.L))
    iP = Pbins.index(h.alpha)
    if h.alpha>80:
        iP = Pbins.index(abs(h.alpha-160))

    if (h.day,iL,iP) in taus:
        Significance[0] = (h.counts-means[h.day,iL,iP])/taus[h.day,iL,iP]
        Tau[0] = taus[h.day,iL,iP]
        Mean[0] = means[h.day,iL,iP]
    else:
        Significance[0] = 0
        Tau[0] = 0
        Mean[0] = 0

    newBranch.Fill()
    newBranch2.Fill()
    newBranch3.Fill()

t3.Write("events", r.TObject.kOverwrite) # save only the new version of the tree
f.Close()
