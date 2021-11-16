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
parser.add_argument('--useVersion', type=str, default='v2', choices=['v1','v2','v2.1'], help='Define wether v1/ (no flux=0) or v2/ (all fluxes), or v2.1/ (all fluxes, summed over energy) is written.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

Rfilename = args.inputFile

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
print("input files has ",entries, " entries")

# prepare root output
#outRootDir = os.path.split(filename)[0]
rootName = os.path.split(Rfilename)[1]
outRootName = sharedOutPath()+"/data/root/"+args.useVersion+"/"+args.data+"/"+rootName
month = rootName[9:15]
day = rootName[9:17]
print(day)

print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )

# pipeline root leafs

L = array( 'f', [ 0. ] )
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
tree.Branch( 'Long', Lo, 'Long/I' )# this was int
tree.Branch( 'Lat', La, 'Lat/I' )  # this was int
tree.Branch( 'geomLong', geomLo, 'geomLong/I' )# this was int
tree.Branch( 'geomLat', geomLa, 'geomLat/I' )  # this was int
tree.Branch( 'field', B, 'field/F' )
tree.Branch( 'field_eq', B_eq, 'field_eq/F' )
tree.Branch( 'altitude', Alt, 'altitude/F' )

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
#print("Pitches: [")
#for i in range(0,36):
#    print(str(dset_p[0][i])+', ')

# write 2d histograms
hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),8,0,180)
hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,np.array(l_x_bins),8,0,180)
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
for evt in t: 
    
    if ( countEv%10000 ==0 ) :
     print (countEv," over ",entries)
     #break

    # calculate daytime in h from msec
    daytime = evt.msec/1e3/60/60
    # lingitude in CSES way..
    longitude = int(evt.lon)
    geomLongitude = int(evt.maglon)
    if evt.lon>180:
        longitude = int(evt.lon) - 360
        geomLongitude = int(evt.maglon) - 360


    # get B field strength at the equator
    Beq = getBeq(evt.L_IGRF) 
    # in case that Beq>B
    if Beq > abs(evt.btotsat):
       Beq = abs(evt.btotsat)
       #if args.debug:
       # print(Beq,evt.L_IGRF, evt.L_IGRF)
       # print("Attention! Beq>B")
    # calculate equatorial pitch angle
    alpha_eq = getAlpha_eq( evt.alpha0sat, abs(evt.btotsat), Beq )

    # find corresponding L/alpha bin, and use the average Lshell and alpha values 
    # only if L shell values < 10                                                                                                                                                                                       
    if evt.L_IGRF < 10 and evt.L_IGRF> 1:
        Lbin=-1
        Albin=-1
        for il in range(l_bins):
            if l_x_bins[il] <= evt.L_IGRF  and l_x_bins[il+1] > evt.L_IGRF:
                Lbin=il
                break
        for ia in range(p_bins):
            if p_x_bins[ia] <= alpha_eq and p_x_bins[ia+1] > alpha_eq:
                Albin=ia
                break
        if Lbin==-1 or Albin==-1 :
            print('B:    ', evt.btotsat)
            print('Beq:  ', Beq)
            print('Alpha:', alpha_eq)
            print('L:    ', evt.L_IGRF)
    
        for ien,en in enumerate(energyTab):
            value = eval('evt.mep_ele_tel0_flux_e'+str(ien+1))
            vecCells[Lbin,Albin].push_back(value)
            vecCellsEn[Lbin,Albin].push_back(en)

            # fill histograms                                                                                                                                                                                           
            hist2D_l_pitch.Fill(l_x_bins[Lbin], p_x_bins[Albin], value)
            hist2D_l_pitch_en.Fill(l_x_bins[Lbin], p_x_bins[Albin])
            hist2D_loc_flux.Fill(longitude, int(evt.lat), value)
            # fill 2D histograms / event                                                                                                                                                                                               
            # time of half-orbit
            bint = hist2D_loc_field.GetBin(hist2D_loc_field.GetXaxis().FindBin(longitude),hist2D_loc_field.GetYaxis().FindBin(int(evt.lat)),0)
            if hist2D_loc.GetBinContent(bint)==0.:
                hist2D_loc.SetBinContent(bint, float(daytime))
            # B field of the earth
            if hist2D_loc_field.GetBinContent(bint)==0.:
                hist2D_loc_field.SetBinContent(bint, abs(evt.btotsat))

        
    for ien,en in enumerate(energyTab):
        F_vecvec.push_back(eval('evt.mep_ele_tel0_flux_e'+str(ien+1)))
        E_vec.push_back(en)
        P_vec.push_back(evt.alpha0sat)
        A_vec.push_back(alpha_eq) 
        #Ch_vec.push_back(vecChannel[key])
    
    Ev[0] = countEv
    L[0] = evt.L_IGRF
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

    countEv += 1

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
    # reset counter
    countEv += 1
    countIntSec = 0
    countInt = 0
    countFlux = 0

outRoot.Write()
outRoot.Close() 
