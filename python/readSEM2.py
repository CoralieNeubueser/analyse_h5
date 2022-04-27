import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
import pathlib
# load defined functions
from utils import *
from drawFunctions import *
from collections import defaultdict

r.gStyle.SetPadRightMargin(0.2)

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--data', type=str, choices=['noaa'], required=True, help='Define patht to data file.')
parser.add_argument('--satellite', type=int, choices=[19], default=19, required=True, help='Define satellite.')
parser.add_argument('--telescope', type=int, choices=[0,90], default=0, required=True, help='Define telescope.')
parser.add_argument('--useVersion', type=str, default='v2', choices=['v1','v2','v2.1'], help='Define wether v1/ (no flux=0) or v2/ (all fluxes), or v2.1/ (all fluxes, summed over energy) is written.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

Rfilename = args.inputFile
integral=-1
if args.integral:
    print("Integrating events over {}s.".format(args.integral))
    integral = args.integral
# define expected energy bins
energy_bins, energyTab, energyMax = getEnergyBins(args.data, False)
# L bins
l_bins, l_x_bins = getLbins()
# pitch bins
p_bins, p_x_bins = getPitchBins()
# earth radius at equator in km
RE = 6378.137 

# read tree
f = r.TFile( Rfilename , 'read' )
t = f.Get("T")
entries = t.GetEntries()
# get measurement rate
times = []
for entr in range(0,10):
    t.GetEntry(entr)
    times.append(t.msec/1e3)
diffTime = int(round(np.amin(np.diff(np.array(times)))))
print("Fluxes are measured every {}s. ".format(diffTime))
xMeasurements=1.
if integral!=-1:
    xMeasurements=int(integral/diffTime)
    if integral%diffTime!=0:
        raise ValueError("Integration window of {}s not possible due to {}s measurement frequency.".format(integral, diffTime))
    print("Fluxes are integrated in {}s by summing over {} measurements.".format(integral, xMeasurements))

print("input files has ",entries, " entries")

# prepare root output
rootName = os.path.split(Rfilename)[1]
det = '{0}_poes{1}_{2}degree'.format(args.data,args.satellite,args.telescope)
outRootDir = sharedOutPath()+"/data/root/"+args.useVersion+"/"+det+"/"
if integral!=-1:
    outRootDir = sharedOutPath()+"/data/root/"+args.useVersion+"/"+det+"/"+str(args.integral)+"s/"
outRootName = outRootDir+rootName

month = rootName[9:15]
day = rootName[9:17]
print(day)

print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )

# pipeline root leafs

L = array( 'f', [ 0. ] )
Lorig = array( 'f', [ 0. ] )
Ldip = array( 'f', [ 0. ] )
N = array('i', [0])
P_vec = r.std.vector(float)()#it was int
A_vec = r.std.vector(float)()
E_vec = r.std.vector(float)()
C = array( 'f', [ 0. ] )
F_vec_en = r.std.vector(float)(energy_bins)
F_vec_pt = r.std.vector(float)(9)
F_vecvec = r.std.vector(float)()#why ?
T = array( 'f', [ 0. ] )
Tday = array( 'i', [0] )
Lo = array( 'i', [ 0 ] ) 
La = array( 'i', [ 0 ] ) 
geomLo = array( 'i', [ 0 ] ) 
geomLa = array( 'i', [ 0 ] ) 
B = array( 'f', [ 0. ] )
B_eq = array( 'f', [ 0. ] )
Ev = array( 'i', [ 0 ] )
Ch_vec = r.std.vector(int)()
Alt = array( 'f', [ 0. ] )

tree.Branch( 'event', Ev, 'event/I' )
tree.Branch( 'channel', Ch_vec)
tree.Branch( 'L', L, 'L/F' )
tree.Branch( 'Lorig', Lorig, 'Lorig/F' )
tree.Branch( 'Ldip', Ldip, 'Ldip/F' )
tree.Branch( 'pitch', P_vec)
tree.Branch( 'alpha', A_vec)
tree.Branch( 'energy', E_vec) 
tree.Branch( 'count', C, 'count/F' )
tree.Branch( 'nflux', N, 'nflux/I' )
tree.Branch( 'flux_en', F_vec_en) 
tree.Branch( 'flux_pt', F_vec_pt) 
tree.Branch( 'flux', F_vecvec) 
tree.Branch( 'time', T, 'time/F' )
tree.Branch( 'day', Tday, 'day/I' )
tree.Branch( 'Long', Lo, 'Long/I' )
tree.Branch( 'Lat', La, 'Lat/I' )
tree.Branch( 'geomLong', geomLo, 'geomLong/I' )
tree.Branch( 'geomLat', geomLa, 'geomLat/I' )
tree.Branch( 'field', B, 'field/F' )
tree.Branch( 'field_eq', B_eq, 'field_eq/F' )
tree.Branch( 'altitude', Alt, 'altitude/F' )

if args.integral:
    Norm = r.std.vector(float)()
    tree.Branch( 'norm', Norm )

vecCells = {} 
for cell_l in range(0,len(l_x_bins)-1):
    for cell_p in range(0,len(p_x_bins)-1):
        vecCells[cell_l,cell_p] = r.std.vector(float)()
        tree.Branch( 'flux_'+str(l_x_bins[cell_l])+'_'+str(p_x_bins[cell_p]), vecCells[cell_l,cell_p]) 

vecCellsEn = {} 
for cell_l in range(0,len(l_x_bins)-1):
    for cell_p in range(0,len(p_x_bins)-1):
         vecCellsEn[cell_l,cell_p] = r.std.vector(float)()
         tree.Branch( 'energy_'+str(l_x_bins[cell_l])+'_'+str(p_x_bins[cell_p]), vecCellsEn[cell_l,cell_p])

print("L-alpha map has {} cells.".format(len(vecCells)))

# write 2d histograms
hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),9,0,180)
hist2D_l_pitch_en_zero=r.TH2D("hist2D_l_pitch_en_zero","hist2D_l_pitch_en_zero",l_bins,np.array(l_x_bins),9,0,180)
hist2D_l_pitch_en_zero_noSAA=r.TH2D("hist2D_l_pitch_en_zero_noSAA","hist2D_l_pitch_en_zero_noSAA",l_bins,np.array(l_x_bins),9,0,180)
hist2D_loc_L=r.TH2D("hist2D_loc_L","hist2D_loc_L",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_Lorig=r.TH2D("hist2D_loc_Lorig","hist2D_loc_Lorig",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

BfieldSum = float(0.)
BeqSum = float(0.)
vec_en = r.std.vector(float)(energy_bins)
vec_nEn = r.std.vector(int)(energy_bins)
vec_pt = r.std.vector(float)(9)
vec_nPt = r.std.vector(int)(9)
vecSum = defaultdict(dict)
vecAlphaL = defaultdict(dict)
vecPitch = defaultdict(dict)
vecChannel = defaultdict(dict)

countFlux = int(0)

#loop on input file entries
countEv = 1
countIntSec = 0
for iev,evt in enumerate(t): 

    if ( countEv%10000 ==0 ) :
        print (countEv," over ",entries)
        #break

    if (math.floor(float(iev)/xMeasurements) < countEv):
        countIntSec += 1

        # calculate daytime in h from msec
        daytime = evt.msec/1e3/60/60
        # lingitude in CSES way..
        longitude = int(evt.lon)
        geomLongitude = int(evt.maglon)
        if evt.lon>180:
            longitude = int(evt.lon) - 360
            geomLongitude = int(evt.maglon) - 360

        Lshell = evt.L_IGRF
        radiusInER = ( RE + evt.alt )/RE
        Ldipole = calculateL(evt.maglat, radiusInER)
        if Lshell==-1:
            Lshell = Ldipole
            #print(Lshell)
        
        # get B field strength at the equator
        Beq = getBeq(Lshell) 
        # in case that Beq>B
        if Beq > abs(evt.btotsat):
            Beq = abs(evt.btotsat)
        # calculate equatorial pitch angle
        alpha_eq = getAlpha_eq( eval('evt.alpha{0}sat'.format(args.telescope)), abs(evt.btotsat), Beq )

        Pvalue = eval('evt.alpha{0}sat'.format(args.telescope))
        # find corresponding L/alpha bin, and use the average Lshell and alpha values 
        # only if L shell values < 10                                                                                                                                                                                       
        if Lshell < 10 and Lshell >= 1:
            Albin=-1
            for ia in range(p_bins):
                if p_x_bins[ia] <= alpha_eq and p_x_bins[ia+1] > alpha_eq:
                    Albin=ia
                    break
            if Albin==-1 :
                print('B:    ', evt.btotsat)
                print('Beq:  ', Beq)
                print('Alpha:', alpha_eq)
                print('L:    ', Lshell)
    

            for ien,en in enumerate(energyTab):
                flux = eval('evt.mep_ele_tel{0}_flux_e{1}'.format(args.telescope,(ien+1)))
                fluxSquared = flux*flux

                if (ien,Albin) in vecSum:
                    vecSum[(ien,Albin)] += flux #Squared
                    vecAlphaL[(ien,Albin)] = [vecAlphaL[(ien,Albin)][0]+alpha_eq, round(vecAlphaL[(ien,Albin)][1]+Lshell,1), vecAlphaL[(ien,Albin)][2]+1.]
                    vecPitch[(ien,Albin)] += Pvalue
                else:
                    vecSum[(ien,Albin)] = flux #Squared
                    vecAlphaL[(ien,Albin)] = [alpha_eq, round(Lshell,1), 1.]
                    vecPitch[(ien,Albin)] = Pvalue

        # add next event
        if countIntSec<xMeasurements:
            continue

    # when finished to loop over integral*s, normalise and write new fluxes
    if args.debug and countIntSec!=1:
        print("Integrated over:", countIntSec)

    # fill flux branches of L-pitch
    for cell,value in vecSum.items():
        # not interested in L shell values > 10
        if vecAlphaL[cell][1]/vecAlphaL[cell][2] >= 10:
            continue
        # find corresponding L/alpha bin, and use the average Lshell and alpha values
        Lbin=-1
        Albin=-1
        for il in range(l_bins):
            if l_x_bins[il] <= vecAlphaL[cell][1]/vecAlphaL[cell][2] and l_x_bins[il+1] > vecAlphaL[cell][1]/vecAlphaL[cell][2]:
                Lbin=il
                break
        for ia in range(p_bins):
            if p_x_bins[ia] <= vecAlphaL[cell][0]/vecAlphaL[cell][2] and p_x_bins[ia+1] > vecAlphaL[cell][0]/vecAlphaL[cell][2]:
                Albin=ia
                break
        if Lbin==-1 or Albin==-1:
            print('Alpha: ',vecAlphaL[cell][0]/vecAlphaL[cell][2])
            print('L:     ',vecAlphaL[cell][1]/vecAlphaL[cell][2])

        # test normalisation by entries in L-alpha not integration window
        vecCells[Lbin,Albin].push_back(value / vecAlphaL[cell][2] ) #countIntSec)
        vecCellsEn[Lbin,Albin].push_back( energyTab[cell[0]] )
        
        # fill histograms                 
        hist2D_l_pitch.Fill(l_x_bins[Lbin], p_x_bins[Albin], value / vecAlphaL[cell][2]) #float(countIntSec))
        hist2D_loc_flux.Fill(longitude, int(evt.lat),  value / vecAlphaL[cell][2]) #/ float(countIntSec))
        
        if cell[0]==0:
            hist2D_l_pitch_en_zero.Fill(l_x_bins[Lbin], p_x_bins[Albin])
        if cell[0]==0 and evt.btotsat>22000:
            hist2D_l_pitch_en_zero_noSAA.Fill(l_x_bins[Lbin], p_x_bins[Albin])

        # fill 2D histograms / event                                                             
        # time of half-orbit
        bint = hist2D_loc_field.GetBin(hist2D_loc_field.GetXaxis().FindBin(longitude),hist2D_loc_field.GetYaxis().FindBin(int(evt.lat)),0)
        if hist2D_loc.GetBinContent(bint)==0.:
            hist2D_loc.SetBinContent(bint, float(daytime))
        # B field of the earth
        if hist2D_loc_field.GetBinContent(bint)==0.:
            hist2D_loc_field.SetBinContent(bint, abs(evt.btotsat))
        # L-shell of the earth
        if hist2D_loc_L.GetBinContent(bint)==0.:
            hist2D_loc_L.SetBinContent(bint, Lshell)
       # L-shell (originally stored) of the earth
        if hist2D_loc_Lorig.GetBinContent(bint)==0. and evt.L_IGRF!=-1:
            hist2D_loc_Lorig.SetBinContent(bint, evt.L_IGRF)

        
    # fill the vector 'flux' and the corresponding 'energy'/'pitch'/'alpha' vectors
    for (key, value) in vecSum.items():
        F_vecvec.push_back(value / float(vecAlphaL[key][2]))#float(countIntSec))
        E_vec.push_back(energyTab[key[0]])
        P_vec.push_back(int(vecPitch[key] / float(vecAlphaL[key][2]))) # normalise if 2 pitch angles ended up in same alpha cell /s
        A_vec.push_back(vecAlphaL[key][0] / vecAlphaL[key][2]) # normalise if 2 pitch angles ended up in same alpha cell /s
        #Ch_vec.push_back(vecChannel[key])
        if value!=0:
            countFlux+=1
        if args.integral:
            Norm.push_back(float(vecAlphaL[key][2]))
    
    Ev[0] = countEv
    L[0] = Lshell
    Lorig[0] = evt.L_IGRF
    Ldip[0] = Ldipole
    T[0] = daytime
    Tday[0] = int(day) 
    C[0] = -1.
    # NOTE: lon,lat and geomag related are float in the native data
    # here converted to int ...
    Lo[0] = longitude
    La[0] = int(evt.lat)
    geomLo[0] = geomLongitude
    geomLa[0] = int(evt.maglat)
    B[0] = abs(evt.btotsat)
    B_eq[0] = abs(Beq)
    N[0] = energy_bins
    Alt[0] = evt.alt
    
    # fill tree
    tree.Fill()
    # clean-up
    vecSum.clear()
    vecPitch.clear()
    vecAlphaL.clear()
    vecChannel.clear()
    E_vec.clear()
    P_vec.clear()
    A_vec.clear()
    Ch_vec.clear()
    F_vec_en.clear()
    vec_en.clear()
    #vec_pt.clear()
    #F_vec_pt.clear()
    F_vecvec.clear()
    #F_vec_pt.resize(9)
    F_vec_en.resize(energy_bins)
    vec_pt.resize(9)
    vec_en.resize(energy_bins)
    for vec in vecCells:
        vecCells[vec].clear()
    for vec in vecCellsEn:
        vecCellsEn[vec].clear()
    if args.integral:
        Norm.clear()

    if iev/xMeasurements<5 and args.debug:
        print("Day:                    ", day)
        print("B field [nT]:           ", evt.btotsat)
        print("L-value:                ", Lshell)
        print("LON/LAT:                ", longitude, int(evt.lat))
        print("GMLON/GMLAT:            ", geomLongitude, int(evt.maglat))
        print("Beq [nT]:               ", round(Beq,2))
        print("Count:                  ", countEv)
        if int(iev/xMeasurements)==4:
            print("ATTENTION!!! in debug mode only 4 events will be processed! use -q")
            break

    # reset counter
    countEv += 1
    countIntSec = 0
    countInt = 0
    countFlux = 0

prep2D(hist2D_l_pitch_en_zero, 'L value', '#alpha_eq [deg]', '#fluxes', False)
prep2D(hist2D_l_pitch_en_zero_noSAA, 'L value', '#alpha_eq [deg]', '#fluxes no SAA', False)
prep2D(hist2D_l_pitch, 'L value', '#alpha_eq [deg]', '#sum#Phi', False)
prep2D(hist2D_loc_L,  'Longitude', 'Latitude', 'L', False)
prep2D(hist2D_loc_Lorig,  'Longitude', 'Latitude', 'L_{orig}', False)
prep2D(hist2D_loc,  'Longitude', 'Latitude', '#entries', False)
prep2D(hist2D_loc_flux,  'Longitude', 'Latitude', '#sum#Phi', False)
prep2D(hist2D_loc_field, 'Longitude', 'Latitude', 'B [nT]', False)

outRoot.Write()
outRoot.Close() 
os.system('chmod -R g+rwx %s'%(outRootName))
