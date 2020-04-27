import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
# load defined functions
from utils import *

r.gStyle.SetPadRightMargin(0.2)

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--data', type=str, choices=['hepd','hepp'], required=True, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

f = h5py.File(filename, 'r')

# read h5 file
if args.debug:
    print(list(f.keys()))
for dset in traverse_datasets(f):
    if args.debug:
        print(dset, f[dset].shape, f[dset].dtype)
    
parameters = dict([('hepd', ['L_parameter', 'HEPD_ele_energy_table', 'HEPD_ele_pitch_table', 'HEPD_ele_energy_pitch', 'UTCTime', 'HEPD_ele_counts','B']),
                   ('hepp', ['L_parameter', 'Energy_Table_Electron', 'PitchAngle', 'A411', 'UTC_TIME', 'Count_Electron', 'CyclotronFrequency_Electron']),
               ])
lonlat = dict([('hepd', ['LonLat', 'LonLat'] ),
               ('hepp', ['GEO_LON', 'GEO_LAT'] )
           ])

dset1 = f[parameters[args.data][0]][()] 
dset_lon = f[lonlat[args.data][0]][()]
dset_lat = f[lonlat[args.data][1]][()]
dset_en = f[parameters[args.data][1]][()]
dset_p = f[parameters[args.data][2]][()]
dset2 = f[parameters[args.data][3]][()]
dset_time = f[parameters[args.data][4]][()]
dset_count = f[parameters[args.data][5]][()]
dset_field = f[parameters[args.data][6]][()]

maxEv = len(dset2)
if args.debug:
    print("Events: ", maxEv)

time_blanc = dset_time[0]
time_blanc_min = dset_time[maxEv-1]
energy_bins = 12
if args.data=='hepp':
    time_blanc = dset_time[0][0]
    time_blanc_min = dset_time[maxEv-1][0]
    energy_bins = 256

time_min = int(str(time_blanc)[-6:-4])*60*60 +  int(str(time_blanc)[-4:-2])*60 +  int(str(time_blanc)[-2:])
time_max = int(str(time_blanc_min)[-6:-4])*60*60 +  int(str(time_blanc_min)[-4:-2])*60 +  int(str(time_blanc_min)[-2:])

badRun=False
# compare time of first and last event and only use here in case that full half-orbit is stored
if (time_max-time_min)/60<30: # in minutes
    badRun=True
    sys.exit(("This file {} is not used, due to incomplete semi-orbit. ").format(filename))

# prepare root output                            
outRootName = os.path.split(filename)[1].replace("h5","root")
outRootName = home()+"/root/"+outRootName
if args.debug:
    print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )
L = array( 'f', [ 0. ] )
P = array( 'f', [ 0. ] )
E = array( 'f', [ 0. ] )
C = array( 'f', [ 0. ] )
F_vec =r.std.vector(float)(energy_bins)
T = array( 'f', [ 0. ] )
Tday = array( 'f', [ 0. ] )
Lo = array( 'i', [ 0 ] )
La = array( 'i', [ 0 ] )
B = array( 'f', [ 0. ] )
            
tree.Branch( 'L', L, 'L/F' )
tree.Branch( 'pitch', P, 'pitch/F' )
tree.Branch( 'energy', E, 'energy/F' )
tree.Branch( 'count', C, 'count/F' )
tree.Branch( 'flux_en', F_vec )
tree.Branch( 'time', T, 'time/F' )
tree.Branch( 'day', Tday, 'day/F' )
tree.Branch( 'Long', Lo, 'Long/I' )
tree.Branch( 'Lat', La, 'Lat/I' )
tree.Branch( 'field', B, 'field/F' )

# write 2d histograms
hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",18,1,10,len(dset_p[0]),np.amin(dset_p[0])-0.5*(dset_p[0][1]-dset_p[0][0]),np.amax(dset_p[0])+0.5*(dset_p[0][8]-dset_p[0][7]))
hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

for iev,ev in enumerate(dset2):
    lonInt = int(0)
    latInt = int(0)
    Bfield = float(0.)
    if args.data=='hepd':
        # fill tree and histograms for HEPD data
        time_calc = 60*60*int(str(dset_time[iev])[-6:-4]) + 60*int(str(dset_time[iev])[-4:-2]) + int(str(dset_time[iev])[-2:])
        time_act = (time_calc-time_min)/60.
        lonInt = int(dset_lon[iev][0])
        latInt = int(dset_lat[iev][1])
        Bfield = dset_field[iev]

    elif args.data=='hepp':
        # fill tree and histos for HEPP data 
        time_calc = 60*60*int(str(dset_time[iev][0])[-6:-4]) + 60*int(str(dset_time[iev][0])[-4:-2]) + int(str(dset_time[iev][0])[-2:])
        time_act = (time_calc-time_min)/60.

        lonInt = int(dset_lon[iev])
        latInt = int(dset_lat[iev])
        # translate cyclotron frequency w=qe*B/(2pi*me) 1/s to B
        qe = 1.602176634e-19 # C = 1.602176634×10−19 As
        # 1Gs = e-4 T = e-4kg/(As2)
        me = 9.109e-31 # kg
        # w seems to have been wrongly calculated in T instead of Gauss, or in 10kHz
        # translate in nT
        Bfield = dset_field[iev]*me/qe*2*np.pi*1e9
        


    # fill 2D histograms / event
    # time of half-orbit
    binx = hist2D_loc.GetXaxis().FindBin(lonInt)
    biny = hist2D_loc.GetYaxis().FindBin(latInt)
    bint = hist2D_loc.GetBin(binx,biny,0)
    if hist2D_loc.GetBinContent(bint)==0.:
        hist2D_loc.SetBinContent(bint, float(time_act))    
    # B field of the earth
    bint = hist2D_loc_field.GetBin(hist2D_loc_field.GetXaxis().FindBin(lonInt),hist2D_loc_field.GetYaxis().FindBin(latInt),0)
    if hist2D_loc_field.GetBinContent(bint)==0.: 
        hist2D_loc_field.SetBinContent(bint, Bfield)

    countInt = int()
    if iev==1 and args.debug:
        print("B field [nT]: ", Bfield)
        print("LON/LAT: {}/{}".format(lonInt,latInt) )
        print("L-value: ", dset1[iev])
        if args.data=='hepd':
            print("Count:   ", dset_count[iev])
    countInt = dset_count[iev]

    # loop through energy bins
    for ie,en in enumerate(ev):
        # loop through pitch
        for ip,flux in enumerate(en):
            # fill tree only for non-zero fluxes
            if flux!=0:
                if iev==1 and args.debug:
                    print("--- Energy:  ", dset_en[0][ie])
                    print("--- Pitch:   ", dset_p[0][ip])
                    print("--- Flux:    ", flux)
                    print("--- Time:    ", dset_time[iev])
                                    
                # HEPP data has stores counts per pitch angle  
                if args.data=='hepp':
                    countInt = dset_count[iev][ip]
                
                # fill histograms
                oldbin = hist2D_l_pitch.FindBin(dset1[iev], dset_p[0][ip]) 
                oldCount = hist2D_l_pitch.GetBinContent(oldbin)
                if oldCount != 0.:
                    hist2D_l_pitch.SetBinContent(oldbin, (oldCount+flux)/2.)
                else:
                    hist2D_l_pitch.SetBinContent(oldbin, flux)

                oldbin1 = hist2D_loc_flux.FindBin(lonInt, latInt)
                oldCount1 = hist2D_loc_flux.GetBinContent(oldbin1)
                if oldCount1 != 0.:
                    hist2D_loc_flux.SetBinContent(oldbin1, (oldCount1+flux)/2.)
                else:
                    hist2D_loc_flux.SetBinContent(oldbin1, flux)

                # fill tree
                L[0] = dset1[iev]
                P[0] = dset_p[0][ip]
                F_vec[ie] = flux
                E[0] = dset_en[0][ie]
                T[0] = time_calc/60/60 # in hours
                Tday[0] = int(str(dset_time[iev][0])[-14:-6])
                C[0] = countInt
                Lo[0] = lonInt
                La[0] = latInt
                B[0] = Bfield
                    
                tree.Fill()


# Print out histrograms                        
outpdf = os.path.split(filename)[1]
outpdf = outpdf.replace("h5","pdf")
if args.data=='hepd':
    outpdf = outpdf.replace("CSES_HEP","map")
elif args.data=='hepp':
    outpdf = outpdf.replace("CSES_01_HEP_1","map")
outpdf = home()+"/plots/"+outpdf

if args.debug:
    print("Writing maps to: ", outpdf)
draw2D(hist2D_l_pitch, "L-value", "pitch [deg]", "#LT electron flux#GT [Hz/(cm^{2}#upoint sr)]", 5e-4, 10, outpdf, True)
outpdf = outpdf.replace("map","loc")
draw2D(hist2D_loc, "longitude [deg]", "latitude [deg]", "#Delta t [min]", 1, 35, outpdf, False)
outpdf = outpdf.replace("loc","loc_flux")
draw2D(hist2D_loc_flux, "longitude [deg]", "latitude [deg]", "#LT electron flux#GT [Hz/(cm^{2}#upoint sr)]", 5e-4, 10, outpdf, True)
outpdf = outpdf.replace("flux","field")
draw2D(hist2D_loc_field, "longitude [deg]", "latitude [deg]", "B [nT]", 17000, 55000, outpdf, False)

outRoot.Write()
outRoot.Close()
