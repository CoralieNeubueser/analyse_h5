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

args,_=parser.parse_known_args()

filename = args.inputFile

f = h5py.File(filename, 'r')

# read h5 file
print(list(f.keys()))
for dset in traverse_datasets(f):
    print(dset, f[dset].shape, f[dset].dtype)
    
parameters = dict([('hepd', ['L_parameter', 'HEPD_ele_energy_table', 'HEPD_ele_pitch_table', 'HEPD_ele_energy_pitch', 'UTCTime']),
                   ('hepp', ['L_parameter', 'Energy_Table_Electron', 'PitchAngle', 'A412', 'UTC_TIME']),
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

maxEv = len(dset2)
print("Events: ", maxEv)

time_blanc = dset_time[0]
if args.data=='hepp':
    time_blanc = dset_time[0][0]
time_max = int(str(time_blanc)[-6:-4])*60*60 +  int(str(time_blanc)[-4:-2])*60 +  int(str(time_blanc)[-2:])

# prepare root output                            
outRootName = os.path.split(filename)[1].replace("h5","root")
outRootName = home()+"/root/"+outRootName
print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )
L = array( 'f', [ 0. ] )
P = array( 'f', [ 0. ] )
E = array( 'f', [ 0. ] )
C = array( 'f', [ 0. ] )
T = array( 'f', [ 0. ] )
Clist = []
Lvalue = [1,1.5,2,2.5,3,3.5,4,4.5,5]
for iL in range(1,10):
    for iP in range(0,9):
        Clist.append(array( 'f', [ 0. ] ))
        ind = (iL-1)*9 + iP
        tree.Branch( 'count_'+str(Lvalue[iL-1])+'_'+str(iP), Clist[ind], 'count_'+str(Lvalue[iL-1])+'_'+str(iP)+'/F' )
            
tree.Branch( 'L', L, 'L/F' )
tree.Branch( 'pitch', P, 'pitch/F' )
tree.Branch( 'energy', E, 'energy/F' )
tree.Branch( 'count', C, 'count/F' )
tree.Branch( 'time', T, 'time/F' )

# write 2d histograms
hist2D=r.TH2D("hist2D","hist2D",18,1,10,len(dset_p[0]),np.amin(dset_p[0])-0.5*(dset_p[0][1]-dset_p[0][0]),np.amax(dset_p[0])+0.5*(dset_p[0][8]-dset_p[0][7]))
hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",360,1,360,180,-90,90)

for iev,ev in enumerate(dset2):
    # fill tree and histograms for HEPD data
    if args.data=='hepd':
        time_calc = 60*60*int(str(dset_time[iev])[-6:-4]) + 60*int(str(dset_time[iev])[-4:-2]) + int(str(dset_time[iev])[-2:])
        time_act = (time_calc-time_max)/60.
        hist2D_loc.SetBinContent(int(dset_lon[iev][0]+180), int(dset_lat[iev][1]+90), float(time_act))
    # fill tree and histos for HEPP data
    elif args.data=='hepp':
        time_calc = 60*60*int(str(dset_time[iev][0])[-6:-4]) + 60*int(str(dset_time[iev][0])[-4:-2]) + int(str(dset_time[iev][0])[-2:])
        time_act = (time_calc-time_max)/60.
        hist2D_loc.SetBinContent(int(dset_lon[iev]+180), int(dset_lat[iev]+90), float(time_act))

    # loop through energy bins
    for ie,en in enumerate(ev):
        # loop through pitch
        for ip,count in enumerate(en):
            if iev==1 and ie==0 and ip==0:
                print("L-value: ", dset1[iev])
                print("Pitch:   ", dset_p[0][ip])
                print("Count:   ", count)
                print("Energy: ", dset_en[0][ie])
            if count!=0: # and dset1[iev]<=10:
                oldbin = hist2D.FindBin(dset1[iev], dset_p[0][ip]) 
                oldCount = hist2D.GetBinContent(oldbin)
                hist2D.SetBinContent(oldbin, (oldCount+count)/2.)
                # fill tree
                L[0] = dset1[iev]
                P[0] = dset_p[0][ip]
                C[0] = count
                E[0] = dset_en[0][ie]
                T[0] = time_act
                
                # if L-value <=5
                if dset1[iev]<=5:  
                    # get indices
                    ind_L = int(dset1[iev]*2)
                    ind_tot = (ind_L-1)*9 + ip
                    Clist[ind_tot][0] = count
                    
                tree.Fill()    
                        
outpdf = os.path.split(filename)[1]
outpdf = outpdf.replace("h5","pdf")
if args.data=='hepd':
    outpdf = outpdf.replace("CSES_HEP","map")
elif args.data=='hepp':
    outpdf = outpdf.replace("CSES_01_HEP_1","map")
outpdf = home()+"/plots/"+outpdf
print("Writing maps to: ", outpdf)

draw2D(hist2D, "L-value", "pitch [deg]", "#LT electron rate#GT [Hz/(cm^{2}#upoint sr)]", 5e-4, 10, outpdf, True)

outpdf = outpdf.replace("map","loc")
draw2D(hist2D_loc, "longitude [deg]", "latitude [deg]", "#Delta t [min]", 1, 35, outpdf, False)

outRoot.Write()
outRoot.Close()
