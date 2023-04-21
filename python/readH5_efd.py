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
parser.add_argument('--data', type=str, choices=['efd_vlf'], required=True, help='Define patht to data file.')
parser.add_argument('--useVersion', type=str, default='v2', choices=['v2'], help='Define wether v1/ (no flux=0) or v2/ (all fluxes), or v2.1/ (all fluxes, summed over energy) is written.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile
f = h5py.File(filename, 'r')

# read h5 file
if args.debug:
    print(list(f.keys()))

version=args.useVersion
integral=1

for dset in traverse_datasets(f):
    if args.debug:
        print(dset, f[dset].shape, f[dset].dtype)
    
parameters = dict([
                   # for EFD no L_shell, pitch angle, Bfield.
                   # E_x, E_y, E_z, Freq, UTC, ALTITUDE
                   ('efd_vlf', ['A131_P', 'A132_P', 'A133_P', 'FREQ', 'UTC_TIME','ALTITUDE','WORKMODE'])
               ])
lonlat = dict([
               ('efd_vlf', ['GEO_LON', 'GEO_LAT'] ),
           ])
gmlonlat = dict([
                 ('efd_vlf', ['MAG_LON', 'MAG_LAT'] ),
             ])

dset_lon = f[lonlat[args.data][0]][()]
dset_lat = f[lonlat[args.data][1]][()]
dset_gmlon = f[gmlonlat[args.data][0]][()]
dset_gmlat = f[gmlonlat[args.data][1]][()]
dset_e_x = f[parameters[args.data][0]][()]
dset_e_y = f[parameters[args.data][1]][()]
dset_e_z = f[parameters[args.data][2]][()]
dset_freq = f[parameters[args.data][3]][()]
dset_time = f[parameters[args.data][4]][()]
dset_alt = f[parameters[args.data][5]][()]
dset_mode = f[parameters[args.data][6]][()]

maxEv = len(dset_time)
if args.debug:
    print("Events: ", maxEv)
time_blanc = dset_time[0]
time_blanc_min = dset_time[maxEv-1]

# prepare root output
outRootDir = os.path.split(filename)[0]
print(outRootDir)

rootName = os.path.split(filename)[1].replace("h5","root")
outRootName = sharedOutPath()+"/data/root/"+version+"/"+args.data+"/"+rootName
if integral!=1:
    outRootName = sharedOutPath()+"/data/root/"+version+"/"+args.data+"/"+str(args.integral)+"s/"+rootName

print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )

Night = array( 'i', [ 0 ] )
C = array( 'f', [ 0. ] )
E_x_vec = r.std.vector(float)()
E_y_vec = r.std.vector(float)()
E_z_vec = r.std.vector(float)()
Freq_vec = r.std.vector(float)()
T = array( 'f', [ 0. ] )
Tday = array( 'i', [0] )
Lo = array( 'f', [ 0. ] )
La = array( 'f', [ 0. ] )
geomLo = array( 'f', [ 0. ] )
geomLa = array( 'f', [ 0. ] )
Alt = array( 'f', [ 0. ] )
Ev = array( 'i', [ 0 ] )
Orbit = array( 'i', [ 0 ] )
Mode = array( 'i', [ 0 ] )

tree.Branch( 'event', Ev, 'event/I' )
tree.Branch( 'orbit', Orbit, 'orbit/I' )
tree.Branch( 'ascending', Night, 'ascending/I' )
tree.Branch( 'Ex', E_x_vec) 
tree.Branch( 'Ey', E_y_vec)
tree.Branch( 'Ez', E_z_vec) 
tree.Branch( 'freq', Freq_vec)
tree.Branch( 'time', T, 'time/F' )
tree.Branch( 'day', Tday, 'day/I' )
tree.Branch( 'Long', Lo, 'Long/F' )
tree.Branch( 'Lat', La, 'Lat/F' )
tree.Branch( 'geomLong', geomLo, 'geomLong/F' )
tree.Branch( 'geomLat', geomLa, 'geomLat/F' )
tree.Branch( 'altitude', Alt, 'altitude/F' )
tree.Branch( 'mode', Mode, 'mode/I' )

head, tail = os.path.split(filename)
numbers = re.findall('\d+', tail)
orbit_index = int(numbers[4])

# write 2d histograms
hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)

countFlux = int(0)
countEv = 1

for iev,ev in enumerate(dset_time):
    lon = float(0)
    lat = float(0)
    gmlon = float(0)
    gmlat = float(0)
    day = int(0)
    daytime = float(0.)
    countInt = int(0)
    prevLat = None
                
    if dset_time[iev][0] == -9999:
        print("Something is wrong in input.. date is set to -9999. Take next event. ")
        continue
    # fill tree and histos for HEPP data
    time_file = 60*60*int(str(dset_time[iev][0])[8:10]) + 60*int(str(dset_time[iev][0])[10:12]) + int(str(dset_time[iev][0])[12:14])
    # time from filename is exactly 1 minute off, take the times as stored in file 
    daytime = time_file/60./60. #time_calc/60./60.
    day = int(str(dset_time[iev][0])[:4] + str(dset_time[iev][0])[4:6] + str(dset_time[iev][0])[6:8]) #int(str(time_blanc)[-14:-6])
    
    lon = dset_lon[iev][0] #[0])
    lat = dset_lat[iev][0] #[1])
    if lat>90:
        print("Something is wrong in Lat/Lon of input.. {},{} skip this event.".format(lat, lon))
        continue
    gmlon = dset_gmlon[iev][0]
    gmlat = dset_gmlat[iev][0] #[1])
    
    # determine wether orbit is ascending
    if prevLat:
        if prevLat < lat:
            ascending = 1
        elif prevLat > lat:
            ascending = 0 
    else:
        if iev+1 < len(dset):
            lat_next = dset_lat[iev+1][0]
            if lat < lat_next:
                ascending = 1
            else:
                ascending = 0
        else:
            lat_prev = dset_lat[iev-1][0]
            if lat > lat_prev:
                ascending = 1
            else:
                ascending = 0
        
    if iev<5 and args.debug:
        print("Day:                    ", day)
        print("LON/LAT:                ", lon,lat)
        print("GMLON/GMLAT:            ", gmlon, gmlat)
        print("Count:                  ", dset_count[iev])
        if iev==4:
            print("ATTENTION!!! in debug mode only 4 events will be processed! use -q")
            break
    
    # loop through frequency bins
    for ie,flux in enumerate(dset_e_x[iev]):
        #print(dset_e_y[iev][ie])
        #print(ie, dset_freq[ie], dset_freq[ie][0])
        Freq_vec.push_back(dset_freq[ie][0])
        E_x_vec.push_back(dset_e_x[iev][ie]) 
        E_y_vec.push_back(dset_e_y[iev][ie])
        E_z_vec.push_back(dset_e_z[iev][ie]) 
        
        hist2D_loc_flux.Fill(lon, lat, dset_e_x[iev][ie])                                                                                                                                                                    

    hist2D_loc.Fill(lon, lat, daytime)
    prevLat = lat

    # fill tree with measures / 1s*integral
    # effectively filles the values for the last integral time point
    Ev[0] = countEv
    Night[0] = ascending
    Orbit[0] = orbit_index
    T[0] = daytime # in hours
    Tday[0] = day
    Lo[0] = lon
    La[0] = lat
    geomLo[0] = gmlon
    geomLa[0] = gmlat
    Alt[0] = dset_alt[iev]
    Mode[0] = dset_mode[iev]

    # fill tree
    tree.Fill()
    # clean-up
    E_x_vec.clear()
    E_y_vec.clear()
    E_z_vec.clear()
    Freq_vec.clear()
    # reset counter
    countEv += 1

prep2D(hist2D_loc, 'Longitude', 'Latitude', 'daytime [h]', False)
prep2D(hist2D_loc_flux, 'Longitude', 'Latitude', '#sum Ex', False)

outRoot.Write()
outRoot.Close()
os.system('chmod -R g+rwx %s'%(outRootName))
