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
parser.add_argument('--data', type=str, choices=['hepd','hepp_l','hepp_h'], required=True, help='Define patht to data file.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--useVersion', type=str, default='v2.1', choices=['v1','v2','v2.1'], help='Define wether v1/ (no flux=0) or v2/ (all fluxes), or v2.1/ (all fluxes, summed over energy) is written.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
args,_=parser.parse_known_args()

filename = args.inputFile

f = h5py.File(filename, 'r')

print("Integrating events over {}s.".format(args.integral))

# read h5 file
if args.debug:
    print(list(f.keys()))

version=args.useVersion
integral=1
if args.integral:
    integral=args.integral

for dset in traverse_datasets(f):
    if args.debug:
        print(dset, f[dset].shape, f[dset].dtype)
    
parameters = dict([('hepd', ['L_parameter', 'HEPD_ele_energy_table', 'HEPD_ele_pitch_table', 'HEPD_ele_energy_pitch', 'UTCTime', 'HEPD_ele_counts','B']),
                   ('hepp_l', ['L_parameter', 'Energy_Table_Electron', 'PitchAngle', 'A411', 'UTC_TIME', 'Count_Electron', 'CyclotronFrequency_Electron']),
                   #('hepp', ['L_parameter', 'electron_energy_table', 'PitchAngle', 'electron_energyPitchAngleSpectrum', 'UTCTime', 'electron_counts', 'electron_gyrofrequency'])
                   ('hepp_h', ['L_parameter', 'Energy_Table_Electron', 'PitchAngle', 'A411', 'UTC_TIME', 'Count_Electron', 'CyclotronFrequency_Electron'])

               ])
lonlat = dict([('hepd', ['LonLat', 'LonLat'] ),
               #('hepp', ['LonLat', 'LonLat'] )
               ('hepp_l', ['GEO_LON', 'GEO_LAT'] ),
               ('hepp_h', ['GEO_LON', 'GEO_LAT'] )
           ])
gmlonlat = dict([('hepd', ['GMLonLat', 'GMLonLat'] ),
                 #('hepp', ['GMLonLat', 'GMLonLat'] )
                 ('hepp_l', ['MAG_LON', 'MAG_LAT'] ),
                 ('hepp_h', ['MAG_LON', 'MAG_LAT'] ),
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
energy_bins, energyTab, energyMax = getEnergyBins(args.data, False)
# L bins
l_bins, l_x_bins = getLbins()
# pitch bins
p_bins, p_x_bins = getPitchBins()
# earth radius at equator in km
RE = 6378.137 

if args.data=='hepp_l' or args.data=='hepp_h':
    hepp = True
    head, tail = os.path.split(filename)
    # CSES_01_HEP_1_L02_A4_069070_20190502_144452_20190502_152151_000
    times = re.findall('\d+', tail)
    time_blanc = int(str(times[5]+times[6]))  #int(str(times[1]+times[2])) #str(dset_time[0])
    time_blanc_min = str(times[7]+times[8]) #int(str(times[3]+times[4])) #dset_time[maxEv-1][0]
    energy_bins, energyTab, energyMax = getEnergyBins(args.data, True)

time_min = int(str(time_blanc)[-6:-4])*60*60 +  int(str(time_blanc)[-4:-2])*60 +  int(str(time_blanc)[-2:])

# prepare root output
outRootDir = os.path.split(filename)[0]
print(outRootDir)

if 'L3_test' in outRootDir:
    useDir = sharedOutPath()+"/data/root/"+version+"/L3_test/"+os.path.split(outRootDir)[1]+'/'
    if integral!=1:
        useDir = sharedOutPath()+"/data/root/"+version+"/"+str(args.integral)+"s/L3_test/"+os.path.split(outRootDir)[1]+'/'
    pathlib.Path(useDir).mkdir(parents=True, exist_ok=True) 
    outRootName = useDir+os.path.split(filename)[1].replace("h5","root")
    
else:
    rootName = os.path.split(filename)[1].replace("h5","root")
    outRootName = sharedOutPath()+"/data/root/"+version+"/"+args.data+"/"+rootName
    if integral!=1:
        outRootName = sharedOutPath()+"/data/root/"+version+"/"+args.data+"/"+str(args.integral)+"s/"+rootName

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
geomLo = array( 'i', [ 0 ] )
geomLa = array( 'i', [ 0 ] )
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
tree.Branch( 'geomLong', geomLo, 'geomLong/I' )
tree.Branch( 'geomLat', geomLa, 'geomLat/I' )
tree.Branch( 'field', B, 'field/F' )
tree.Branch( 'field_eq', B_eq, 'field_eq/F' )

# define the L-pitch map
energyBins = getEnergyBins(args.data, False)
if args.data=='hepp_l' or args.data=='hepp_h':
    energyBins = getEnergyBins(args.data, True)

vecCells = {} 
for cell_l in range(0,len(l_x_bins)-1):
    for cell_p in range(0,len(p_x_bins)-1):
        vecCells[cell_l,cell_p] = r.std.vector(float)()
        tree.Branch( 'flux_'+str(l_x_bins[cell_l])+'_'+str(p_x_bins[cell_p]), vecCells[cell_l,cell_p]) 
# not needed when summed over energy
if version!='v2.1':
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
hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),len(dset_p[0]),0,180)
hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,np.array(l_x_bins),len(dset_p[0]),0,180)
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
countInt = int(0)
countFlux = int(0)

countEv = 1
countIntSec = 0
for iev,ev in enumerate(dset2):
    lonInt = int(0)
    latInt = int(0)
    gmlonInt = int(0)
    gmlatInt = int(0)
    Bfield = float(0.)
    Lshell = float(0.)
    day = int(0)
    daytime = float(0.)
    Beq = float(0.)

    # use energy bins as defined in h5 table
    energies = energyTab
    hepd = False

    if (math.floor(float(iev)/float(integral)) < countEv):
        countIntSec += 1

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
            BfieldSum += Bfield
            Lshell = dset1[iev]
            countInt += dset_count[iev]
                
        elif args.data=='hepp_l' or args.data=='hepp_h':
            # fill tree and histos for HEPP data 
            time_calc = time_min + iev #60*60*int(str(dset_time[iev][0])[-6:-4]) + 60*int(str(dset_time[iev][0])[-4:-2]) + int(str(dset_time[iev][0])[-2:])
            time_act = (time_calc-time_min)/60.
            daytime = time_calc/60./60.
            day = int(str(time_blanc)[-14:-6])
            year = int(str(time_blanc)[-14:-10])

            lonInt = int(dset_lon[iev][0]) #[0])
            latInt = int(dset_lat[iev][0]) #[1])
            gmlonInt = int(dset_gmlon[iev][0])
            gmlatInt = int(dset_gmlat[iev][0]) #[1])
            # translate cyclotron frequency w=qe*B/(2pi*me) 1/s to B
            qe = 1.602176634e-19 # C = 1.602176634×10−19 As
            # 1Gs = e-4 T = e-4kg/(As2)
            me = 9.109383701528e-31 # kg
            # w seems to have been wrongly calculated in T instead of Gauss, or in 10kHz
            # translate in nT
            Bfield = dset_field[iev]*me/qe*2*np.pi*1e9
            BfieldSum += Bfield
            Lshell = dset1[iev][0]
            # sum over all channel counts per s
            for channel_count in dset_count[iev]:
                countInt += channel_count

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

        BeqSum += Beq
        
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
                ip_orig = ip
                if version == 'v1' and flux==0:
                    continue

                # rebin the HEPP-L entries from 36 to 9
                if args.data == 'hepp_h':
                    ip = int(ip/4.)

                # fill tree only for non-zero fluxes
                # fill also 0s, decided 2020/10/26
                # correct flux by new geometrical factors
                flux = flux*getGeomCorr(hepd, ie)

                vec_nPt[ip] += 1
            
                # HEPP data has stores counts per 9 different devices (merge all)                 
                # fill energy-flux vector (summ over fluxes over all pitches) 
                if not hepd:
                    if args.data=='hepp_l':
                        countInt = dset_count[iev][ip]
                    else:
                        countInt = dset_count[iev][0]
                    # rebin the energy range from 256 to 16
                    maxE = dset_en[0][255]
                    minE = dset_en[0][0]
                    binE = (maxE-minE)/16.
                    newbin = math.floor(ie/16.)
                    vec_en[newbin] += pow(flux,2)
                    Pvalue = dset_p[0][ip_orig]
                    ie = newbin

                else:
                    Pvalue = (dset_p[0][ip_orig]+dset_p[0][ip_orig-1])/2.
                    if ip==0:
                        Pvalue = dset_p[0][ip_orig]/2.
                    vec_en[ie] += pow(flux,2)
                
                # fill pitch-flux vector (summ over fluxes over all energies, normalise before to MeV, using the energy bin width)
                vec_pt[ip] += pow(float(flux/getEnergyBinWidth(hepd, ie)),2)
                # calculate equatorial pitch angle
                alpha_eq = getAlpha_eq( Pvalue, Bfield, Beq )
                # fill Energy-local pitch matrix
                # store corresponding L/alpha values
                if (ie,ip) in vecSum:
                    vecSum[(ie,ip)] += flux
                    vecAlphaL[(ie,ip)] = [vecAlphaL[(ie,ip)][0]+alpha_eq, vecAlphaL[(ie,ip)][1]+round(Lshell,1), vecAlphaL[(ie,ip)][2]+1.]
                else:
                    vecSum[(ie,ip)] = flux
                    vecAlphaL[(ie,ip)] = [alpha_eq, round(Lshell,1), 1.]

                if iev==1 and args.debug:
                    print("--- Energy bin:      ", ie)
                    print("--- Energy:          ", round(energies[ie],1))
                    print("--- Pitch bin:       ", ip)
                    print("--- Pitch:           ", Pvalue)
                    print("--- Pitch_eq:        ", alpha_eq)
                    print("--- Flux:            ", flux)
                    print("--- Day time [h]:    ", time_calc/60/60 )
        # add next event
        if countIntSec<integral:
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
            print(vecAlphaL[cell][1]/vecAlphaL[cell][2])
        vecCells[Lbin,Albin].push_back(value / countIntSec)
        vecCellsEn[Lbin,Albin].push_back(round(energies[cell[0]],1))
                    
        # fill histograms
        hist2D_l_pitch.Fill(l_x_bins[Lbin], p_x_bins[Albin], float(value)/float(countIntSec))
        hist2D_l_pitch_en.Fill(l_x_bins[Lbin], p_x_bins[Albin])
        hist2D_loc_flux.Fill(lonInt, latInt, float(value)/float(countIntSec))
        # fill 2D histograms / event
        # time of half-orbit
        bint = hist2D_loc_field.GetBin(hist2D_loc_field.GetXaxis().FindBin(lonInt),hist2D_loc_field.GetYaxis().FindBin(latInt),0)
        if hist2D_loc.GetBinContent(bint)==0.:
            hist2D_loc.SetBinContent(bint, float(daytime))
        # B field of the earth
        if hist2D_loc_field.GetBinContent(bint)==0.:
            hist2D_loc_field.SetBinContent(bint, Bfield)

    
    # normalise flux vectors
    numPt = 0
    for ipt,flux_pt in enumerate(vec_pt):
        F_vec_pt[ipt] = math.sqrt(flux_pt) / float(countIntSec)
        numPt += vec_nPt[ipt]

    for ien in range(energy_bins):
        if vec_en[ien]!=0:
            F_vec_en[ien] = math.sqrt(vec_en[ien]) / float(countIntSec)

    # fill the vector 'flux' and the corresponding 'energy'/'pitch'/'alpha' vectors
    for (key, value) in vecSum.items():
        F_vecvec.push_back(value / float(countIntSec))
        E_vec.push_back(round(energies[key[0]],1))
        P_vec.push_back(p_x_bins[key[1]])
        A_vec.push_back(vecAlphaL[key][0]/vecAlphaL[key][2])
        if value!=0:
            countFlux+=1

    # fill tree with measures / 1s*integral
    # effectively filles the values for the last integral time point
    Ev[0] = countEv
    L[0] = Lshell
    T[0] = daytime # in hours
    Tday[0] = day 
    C[0] = countInt
    Lo[0] = lonInt
    La[0] = latInt
    geomLo[0] = gmlonInt
    geomLa[0] = gmlatInt
    B[0] = Bfield
    B_eq[0] = Beq
    N[0] = countFlux
    
    # fill tree
    tree.Fill()
    # clean-up
    vecSum.clear()
    E_vec.clear()
    P_vec.clear()
    A_vec.clear()
    F_vec_en.clear()
    vec_en.clear()
    vec_pt.clear()
    F_vec_pt.clear()
    F_vecvec.clear()
    F_vec_pt.resize(9)
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

prep2D(hist2D_l_pitch, 'L value', '#alpha_{eq} [deg]', '#sum flux', False)
prep2D(hist2D_l_pitch_en, 'L value', '#alpha_{eq} [deg]', '#entries', False)
prep2D(hist2D_loc, 'Longitude', 'Latitude', 'daytime [h]', False)
prep2D(hist2D_loc_flux, 'Longitude', 'Latitude', '#sum flux', False)
prep2D(hist2D_loc_field, 'Longitude', 'Latitude', 'B field [nT]', False)

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
