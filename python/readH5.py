import os,sys,argparse,re
import math
import matplotlib as plt
import numpy as np
from array import array
import pathlib
# load defined functions
from utils import *
from drawFunctions import *

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
                   #('hepp', ['L_parameter', 'Energy_Table_Electron', 'PitchAngle', 'A411', 'UTC_TIME', 'Count_Electron', 'CyclotronFrequency_Electron'])
                   ('hepp', ['L_parameter', 'electron_energy_table', 'PitchAngle', 'electron_energyPitchAngleSpectrum', 'UTCTime', 'electron_counts', 'electron_gyrofrequency'])
               ])
lonlat = dict([('hepd', ['LonLat', 'LonLat'] ),
               ('hepp', ['LonLat', 'LonLat'] )
               #('hepp', ['GEO_LON', 'GEO_LAT'] )
           ])
gmlonlat = dict([('hepd', ['GMLonLat', 'GMLonLat'] ),
                 ('hepp', ['GMLonLat', 'GMLonLat'] )
                 #('hepp', ['MAG_LON', 'MAG_LAT'] )
             ])

### data format of HEPP_august_2018: HEPP_L
###['Altitude', 'GMLonLat', 'L_parameter', 'LonLat', 'ProductAttributes', 'UTCTime', 'electron_counts', 'electron_energyPitchAngleSpectrum', 'electron_energy_table', 'electron_gyrofrequency', 'particle_counts', 'proton_counts', 'proton_energyPitchAngleSpectrum', 'proton_energy_table', 'proton_gyrofrequency']

dset1 = f[parameters[args.data][0]][()] 
dset_lon = f[lonlat[args.data][0]][()]
dset_lat = f[lonlat[args.data][1]][()]
dset_gmlon = f[gmlonlat[args.data][0]][()]
dset_gmlat = f[gmlonlat[args.data][1]][()]
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

### prepare
# define expected energy bins, fill for hepd as default
energy_bins, energyTab, energyMax = getEnergyBins(True, False)
# L bins
l_bins, l_x_bins = getLbins()
# pitch bins
p_bins, p_x_bins = getPitchBins()
# earth radius at equator in km
RE = 6378.137 

if args.data=='hepp':
    head, tail = os.path.split(filename)
    times = re.findall('\d+', tail)
    time_blanc = int(str(times[1]+times[2])) #str(dset_time[0])
    time_blanc_min = int(str(times[3]+times[4])) #dset_time[maxEv-1][0]
    energy_bins, energyTab, energyMax = getEnergyBins(False, True)

time_min = int(str(time_blanc)[-6:-4])*60*60 +  int(str(time_blanc)[-4:-2])*60 +  int(str(time_blanc)[-2:])

# prepare root output
outRootDir = os.path.split(filename)[0]
print(outRootDir)

if 'L3_test' in outRootDir:
    useDir = sharedOutPath()+"/data/root/v2/L3_test/"+os.path.split(outRootDir)[1]+'/'
    pathlib.Path(useDir).mkdir(parents=True, exist_ok=True) 
    outRootName = useDir+os.path.split(filename)[1].replace("h5","root")
    
else:
    outRootName = os.path.split(filename)[1].replace("h5","root")
    outRootName = sharedOutPath()+"/data/root/v2/"+outRootName

print("Writing output into root file: ", outRootName)
outRoot = r.TFile( outRootName , 'recreate' )
tree = r.TTree( 'tree', 'tree with histos' )

L = array( 'f', [ 0. ] )
N = array('i', [0])
P_vec = r.std.vector(int)()
A_vec = r.std.vector(float)()
E_vec = r.std.vector(float)()
C = array( 'f', [ 0. ] )
F_vec_en = r.std.vector(float)(energy_bins)
F_vec_pt = r.std.vector(float)(9)
F_vecvec = r.std.vector(float)()
T = array( 'f', [ 0. ] )
Tday = array( 'i', [0] )
Lo = array( 'i', [ 0 ] )
La = array( 'i', [ 0 ] )
B = array( 'f', [ 0. ] )
B_eq = array( 'f', [ 0. ] )
Ev = array( 'i', [ 0 ] )

tree.Branch( 'event', Ev, 'event/I' )
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
tree.Branch( 'Long', Lo, 'Long/I' )
tree.Branch( 'Lat', La, 'Lat/I' )
tree.Branch( 'field', B, 'field/F' )
tree.Branch( 'field_eq', B_eq, 'field_eq/F' )

# define the L-pitch map
vecCells = []
numCells=0
for cell_l in range(0,len(l_x_bins)-1):
    for cell_p in range(0,len(p_x_bins)-1):
        vecCells.append( r.std.vector(float)() )
        tree.Branch( 'flux_'+str(l_x_bins[cell_l])+'_'+str(p_x_bins[cell_p]), vecCells[numCells]) 
        numCells+=1
print("L-alpha map has {} cells.".format(numCells))

# write 2d histograms
hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),len(dset_p[0]),0,180)
hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,np.array(l_x_bins),len(dset_p[0]),0,180)
hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

for iev,ev in enumerate(dset2):
    lonInt = int(0)
    latInt = int(0)
    gmlonInt = int(0)
    gmlatInt = int(0)
    Bfield = float(0.)
    Lshell = float(0.)
    day = int()
    Beq = float()
    # use energy bins as defined in h5 table
    energies = energyTab
    countInt = int(0)
    hepd = False

    if args.data=='hepd':
        hepd = True
        # fill tree and histograms for HEPD data
        time_calc = 60*60*int(str(dset_time[iev])[-6:-4]) + 60*int(str(dset_time[iev])[-4:-2]) + int(str(dset_time[iev])[-2:])
        time_act = (time_calc-time_min)/60.
        daytime = time_calc/60./60.
        day = int(str(dset_time[iev])[-14:-6])
        year = int(str(dset_time[iev])[-14:-10])

        lonInt = int(dset_lon[iev][0])
        latInt = int(dset_lat[iev][1])
        gmlonInt = int(dset_gmlon[iev][0])
        gmlatInt = int(dset_gmlat[iev][1])
        Bfield = dset_field[iev]
        Lshell = dset1[iev]
        countInt = dset_count[iev]
                
    elif args.data=='hepp':
        # fill tree and histos for HEPP data 
        time_calc = time_min + iev #60*60*int(str(dset_time[iev][0])[-6:-4]) + 60*int(str(dset_time[iev][0])[-4:-2]) + int(str(dset_time[iev][0])[-2:])
        time_act = (time_calc-time_min)/60.
        daytime = time_calc/60./60.
        day = int(str(time_blanc)[-14:-6])
        year = int(str(time_blanc)[-14:-10])
        lonInt = int(dset_lon[iev][0])
        latInt = int(dset_lat[iev][1])
        gmlonInt = int(dset_gmlon[iev][0])
        gmlatInt = int(dset_gmlat[iev][1])
        # translate cyclotron frequency w=qe*B/(2pi*me) 1/s to B
        qe = 1.602176634e-19 # C = 1.602176634×10−19 As
        # 1Gs = e-4 T = e-4kg/(As2)
        me = 9.109383701528e-31 # kg
        # w seems to have been wrongly calculated in T instead of Gauss, or in 10kHz
        # translate in nT
        Bfield = dset_field[iev]*me/qe*2*np.pi*1e9
        Lshell = dset1[iev]
        # sum over all channel counts per s
        for channel_count in dset_count[iev]:
            countInt += channel_count

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

    # get B field strength at the equator
    Beq = getBeq(Lshell) 
    # in case that Beq>B
    if Beq > Bfield:
        if args.debug:
            print(Beq,Bfield)
            print(Lshell)
            print("Attention! Beq>B")
        # difference found for HEPP-L, probably due to Bfield caluculations??
        Beq = Bfield
        #continue

    countFlux = int(0)

    if iev==1 and args.debug:
        print("Day:                    ", day)
        print("B field [nT]:           ", Bfield)
        print("L-value:                ", Lshell)
        print("LON/LAT:                ", lonInt,latInt)
        print("GMLON/GMLAT:            ", gmlonInt, gmlatInt)
        print("Beq [nT]:               ", round(Beq,2))
        print("Count:   ", countInt)
    
    # loop through energy bins
    for ie,en in enumerate(ev):
        # loop through pitch
        for ip,flux in enumerate(en):

            # fill tree only for non-zero fluxes
            # fill also 0s, decided 2020/10/26
            # if float(flux)!=0:
            countFlux+=1
            # correct flux by new geometrical factors
            flux = flux*getGeomCorr(hepd, ie)

            # fill pitch/energy/flux vectors
            #print(F_vec_pt[ip])
            F_vec_pt[ip] += flux
            #print("{}\n".format(F_vec_pt[ip]))
            F_vecvec.push_back(flux)
            
            # HEPP data has stores counts per 9 different devices (merge all)  
            if not hepd:
                countInt = dset_count[iev][ip]
                # rebin the energy range from 256 to 16
                maxE = dset_en[0][255]
                minE = dset_en[0][0]
                binE = (maxE-minE)/16.
                newbin = math.floor(ie/16.)
                E_vec.push_back(float(round(energies[newbin],5)))
                F_vec_en[newbin] += flux
                Pvalue = dset_p[0][ip]
                P_vec.push_back(int(Pvalue))
                ie = newbin

            else:
                Pvalue = (dset_p[0][ip]+dset_p[0][ip-1])/2.
                if ip==0:
                    Pvalue = dset_p[0][ip]/2.
                P_vec.push_back(int(Pvalue))
                E_vec.push_back(float(energies[ie]))
                F_vec_en[ie] += float(flux)
                
            # calculate equatorial pitch angle
            alpha_eq = getAlpha_eq( Pvalue, Bfield, Beq )
            A_vec.push_back( alpha_eq )

            if iev==1 and args.debug:
                print("--- Energy bin:      ", ie)
                print("--- Energy:          ", energies[ie])
                print("--- Pitch bin:       ", ip)
                print("--- Pitch:           ", Pvalue)
                print("--- Pitch_eq:        ", alpha_eq)
                print("--- Flux:            ", flux)
                print("--- Day time [h]:    ", time_calc/60/60 )

            # fill flux branches of L-pitch
            icell = 0
            found = False
            for cell_l in range(0,len(l_x_bins)-1):
                for cell_p in range(0,len(p_x_bins)-1):
                    if l_x_bins[cell_l] < Lshell and l_x_bins[cell_l+1] > dset1[iev]:
                        if p_x_bins[cell_p] < alpha_eq and p_x_bins[cell_p+1] > alpha_eq:
                            found = True
                            vecCells[icell].push_back(flux)
                            break
                    icell+=1
                if found==True:
                    break
                    
            # fill histograms
            hist2D_l_pitch.Fill(Lshell, alpha_eq, flux)
            hist2D_l_pitch_en.Fill(Lshell, alpha_eq)
            hist2D_loc_flux.Fill(lonInt, latInt, flux)

    # fill tree with measures / 1s
    Ev[0] = iev
    L[0] = Lshell
    T[0] = daytime # in hours
    Tday[0] = day 
    C[0] = countInt
    Lo[0] = lonInt
    La[0] = latInt
    B[0] = Bfield
    B_eq[0] = Beq
    N[0] = countFlux

    tree.Fill()

    # clean-up
    E_vec.clear()
    P_vec.clear()
    A_vec.clear()
    F_vec_en.clear()
    F_vec_pt.clear()
    F_vecvec.clear()
    F_vec_pt.resize(9)
    F_vec_en.resize(energy_bins)
    for vec in vecCells:
        vec.clear()

prep2D(hist2D_l_pitch, 'L value', '#alpha_{eq} [deg]', '#sum flux', False)
prep2D(hist2D_l_pitch_en, 'L value', '#alpha_{eq} [deg]', '#entries', False)
prep2D(hist2D_loc_flux, 'Longitude', 'latitude', '#sum flux', False)

# Print out histrograms                        
#outpdf = os.path.split(filename)[1]
#outpdf = outpdf.replace("h5","pdf")
#if args.data=='hepd':
#    outpdf = outpdf.replace("CSES_HEP","map")
#elif args.data=='hepp':
#    outpdf = outpdf.replace("CSES_01_HEP_1","map")
#outpdf = home()+"/plots/"+outpdf

#if args.debug:
#    print("Writing maps to: ", outpdf)
#draw2D(hist2D_l_pitch, "L-value", "pitch [deg]", "#LT electron flux#GT [Hz/(cm^{2}#upoint sr)]", 5e-4, 10, outpdf, True)
#outpdf = outpdf.replace("map","loc")
#draw2D(hist2D_loc, "longitude [deg]", "latitude [deg]", "#Delta t [min]", 1, 35, outpdf, False)
#outpdf = outpdf.replace("loc","loc_flux")
#draw2D(hist2D_loc_flux, "longitude [deg]", "latitude [deg]", "#LT electron flux#GT [Hz/(cm^{2}#upoint sr)]", 5e-4, 10, outpdf, True)
#outpdf = outpdf.replace("flux","field")
#draw2D(hist2D_loc_field, "longitude [deg]", "latitude [deg]", "B [nT]", 17000, 55000, outpdf, False)

outRoot.Write()
outRoot.Close()
