import os,sys,math
import numpy as np
import ROOT as r
import h5py
from drawFunctions import *

def home():
    return os.getcwd()

# read in h5 file
def traverse_datasets(hdf_file):

    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            keystr=str(prefix+"/"+key)
            path = str(prefix+"/"+key) #f'{prefix}'
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    for path, _ in h5py_dataset_iterator(hdf_file):
        yield path

# define L value bins
# returns the number of bins and a list of bins
def getLbins():
    l_x_bins = []
    for x in range(0,5):
        l_x_bins.append(1.+0.2*x)
    for x in range(0,9):
        l_x_bins.append(2.+float(x))
    l_bins=len(l_x_bins)-1
    return l_bins, l_x_bins

# define pitch bins
# returns the number of bins and a list of bins
def getPitchBins():
    p_x_bins = []
    for p in range(0,10):
        p_x_bins.append(p*20)
    p_bins=len(p_x_bins)
    return p_bins, p_x_bins

# define energy bins
# returns 
# 1. the number of energy bins for either hepd of hepp data
# 2. a list of lower edge energy values
# 3. upper bin edge of last energy bin
def getEnergyBins(hepd, hepp):
    if hepd:
        return 12, [2.0, 6.5, 9.9, 12.2, 14.9, 18.2, 22.3, 27.2, 33.3, 40.7, 49.7, 60.8], 70
    elif hepp:
        return 256, [0.1, 0.11137255, 0.1227451,  0.13411765, 0.1454902,  0.15686275,
                     0.1682353 , 0.17960784, 0.1909804 , 0.20235294, 0.21372551, 0.22509804,
                     0.23647058, 0.24784315, 0.25921568, 0.27058825, 0.28196079, 0.29333335,
                     0.30470589, 0.31607845, 0.32745099, 0.33882353, 0.35019609, 0.36156863,
                     0.37294117, 0.38431373, 0.39568627, 0.40705884, 0.41843137, 0.42980394,
                     0.44117647, 0.45254904, 0.46392158, 0.47529411, 0.48666668, 0.49803922,
                     0.50941181, 0.52078432, 0.53215688, 0.54352945, 0.55490202, 0.56627452,
                     0.57764709, 0.58901966, 0.60039222, 0.61176473, 0.6231373 , 0.63450986,
                     0.64588237, 0.65725493, 0.6686275 , 0.68000007, 0.69137257, 0.70274514,
                     0.71411771, 0.72549027, 0.73686278, 0.74823534, 0.75960791, 0.77098042,
                     0.78235298, 0.79372555, 0.80509812, 0.81647062, 0.82784319, 0.83921576,
                     0.85058826, 0.86196083, 0.87333339, 0.88470596, 0.89607847, 0.90745103,
                     0.9188236 , 0.93019611, 0.94156867, 0.95294124, 0.96431381, 0.97568631,
                     0.98705888, 0.99843144, 1.00980401, 1.02117646, 1.03254902, 1.04392159,
                     1.05529416, 1.06666672, 1.07803929, 1.08941185, 1.10078442, 1.11215687,
                     1.12352943, 1.134902  , 1.14627457, 1.15764713, 1.1690197 , 1.18039227,
                     1.19176471, 1.20313728, 1.21450984, 1.22588241, 1.23725498, 1.24862754,
                     1.26000011, 1.27137268, 1.28274512, 1.29411769, 1.30549026, 1.31686282,
                     1.32823539, 1.33960795, 1.35098052, 1.36235297, 1.37372553, 1.3850981 ,
                     1.39647067, 1.40784323, 1.4192158 , 1.43058836, 1.44196081, 1.45333338,
                     1.46470594, 1.47607851, 1.48745108, 1.49882364, 1.51019621, 1.52156866,
                     1.53294122, 1.54431379, 1.55568635, 1.56705892, 1.57843149, 1.58980405,
                     1.6011765 , 1.61254907, 1.62392163, 1.6352942 , 1.64666677, 1.65803933,
                     1.6694119 , 1.68078434, 1.69215691, 1.70352948, 1.71490204, 1.72627461,
                     1.73764718, 1.74901974, 1.76039219, 1.77176476, 1.78313732, 1.79450989,
                     1.80588245, 1.81725502, 1.82862759, 1.84000003, 1.8513726 , 1.86274517,
                     1.87411773, 1.8854903 , 1.89686286, 1.90823543, 1.919608  , 1.93098044,
                     1.94235301, 1.95372558, 1.96509814, 1.97647071, 1.98784328, 1.99921584,
                     2.01058817, 2.02196074, 2.0333333 , 2.04470587, 2.05607843, 2.067451  ,
                     2.07882357, 2.09019613, 2.1015687 , 2.11294127, 2.12431359, 2.13568616,
                     2.14705873, 2.15843129, 2.16980386, 2.18117642, 2.19254899, 2.20392156,
                     2.21529412, 2.22666669, 2.23803926, 2.24941182, 2.26078439, 2.27215695,
                     2.28352928, 2.29490185, 2.30627441, 2.31764698, 2.32901955, 2.34039211,
                     2.35176468, 2.36313725, 2.37450981, 2.38588238, 2.39725494, 2.40862751,
                     2.42000008, 2.43137264, 2.44274521, 2.45411754, 2.4654901 , 2.47686267,
                     2.48823524, 2.4996078 , 2.51098037, 2.52235293, 2.5337255 , 2.54509807,
                     2.55647063, 2.5678432 , 2.57921576, 2.59058833, 2.6019609 , 2.61333323,
                     2.62470579, 2.63607836, 2.64745092, 2.65882349, 2.67019606, 2.68156862,
                     2.69294119, 2.70431376, 2.71568632, 2.72705889, 2.73843145, 2.74980402,
                     2.76117659, 2.77254891, 2.78392148, 2.79529405, 2.80666661, 2.81803918,
                     2.82941175, 2.84078431, 2.85215688, 2.86352944, 2.87490201, 2.88627458,
                     2.89764714, 2.90901971, 2.92039227, 2.9317646 , 2.94313717, 2.95450974,
                     2.9658823 , 2.97725487, 2.98862743, 3.        ], 3.

# returns the geometrical correction factor for hepd data, per energy bin
def getGeomCorr(hepd, energyBin):
    # geometrical factors
    ele_GF = [ 0.76, 188.26, 326.64, 339.65, 344.99, 331.83, 304.73, 263.56, 217.33, 169.48, 117.31, 71.45 ]
    ele_corr_GF = [ 131.128, 545.639, 560.297, 530.937, 477.827, 413.133, 334.176, 252.3, 204.52, 103.216, 77.5552, 61.1536 ]
    if hepd:
        return ele_GF[energyBin]/ele_corr_GF[energyBin]
    else:
        return 1.

# for HEPD data
# returns the 'new' geometrical factor based on MC for Edep, per energy bin  
def getGeomFactor(energyBin):
    # geometrical factors
    ele_corr_GF = [ 131.128, 545.639, 560.297, 530.937, 477.827, 413.133, 334.176, 252.3, 204.52, 103.216, 77.5552, 61.1536 ]
    return ele_corr_GF[energyBin]
    
# helper to merge not only root tree, but the 2D histograms in a file as well
# specify if neede with bool 'allHists'
def merge(name, listOfFiles, runs, allHists):
    ch=r.TChain("tree")  # creates a chain to process a Tree called "tree"

    # L bins
    l_bins, l_x_bins = getLbins()
    # pitch bins              
    p_bins, p_x_bins = getPitchBins()

    hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),9,0,180)
    hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,np.array(l_x_bins),9,0,180)
    hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

    for iF,inFile in enumerate(listOfFiles):
        # limit number of files used for merge
        if iF>runs:
            break
        ch.Add(inFile)

        if allHists:
            # open file, get histogram, add without sum over entries
            inRoot = r.TFile( inFile , 'read' )
            h1=r.TH2D()
            h1_en=r.TH2D()
            h2=r.TH2D()
            h3=r.TH2D()
            h4=r.TH2D()
            r.gDirectory.GetObject('hist2D_l_pitch', h1)
            r.gDirectory.GetObject('hist2D_l_pitch_en', h1_en)
            r.gDirectory.GetObject('hist2D_loc', h2)
            r.gDirectory.GetObject('hist2D_loc_flux', h3)
            r.gDirectory.GetObject('hist2D_loc_field', h4)

            # fill pitch/L histogram with sum
            for binx in range(0,h1.GetNbinsX()+1):
                for biny in range(0,h1.GetNbinsY()+1):
                    bint = h1.GetBin(binx,biny)
                    cont = h1.GetBinContent(bint)
                    oldcont = hist2D_l_pitch.GetBinContent(bint)
                    if cont!=0.:
                        if oldcont!=0.:
                            hist2D_l_pitch.SetBinContent(bint, oldcont+cont)
                        else:
                            hist2D_l_pitch.SetBinContent(bint, cont)

            # fill pitch/L histogram of entries with sum
            for binx in range(0,h1_en.GetNbinsX()+1):
                for biny in range(0,h1_en.GetNbinsY()+1):
                    bint = h1_en.GetBin(binx,biny)
                    cont = h1_en.GetBinContent(bint)
                    oldcont = hist2D_l_pitch_en.GetBinContent(bint)
                    if cont!=0.:
                        if oldcont!=0.:
                            hist2D_l_pitch_en.SetBinContent(bint, oldcont+cont)
                        else:
                            hist2D_l_pitch_en.SetBinContent(bint, cont)

            # loop through all lon/lat
            for binx in range(0,h3.GetNbinsX()+1):
                for biny in range(0,h3.GetNbinsY()+1):
                    bint = h3.GetBin(binx,biny)
                    # fill flux over location histogram with sum
                    cont = h3.GetBinContent(bint)
                    oldcont = hist2D_loc_flux.GetBinContent(bint)
                    if cont!=0.:
                        if oldcont!=0.:
                            hist2D_loc_flux.SetBinContent(bint, oldcont+cont)
                        else:
                            hist2D_loc_flux.SetBinContent(bint, cont)

                    # fill field/time bins, only if not already has a value 
                    cont2 = h2.GetBinContent(bint)
                    cont4 = h4.GetBinContent(bint)
                    if hist2D_loc.GetBinContent(bint)==0.:
                        hist2D_loc.SetBinContent(bint, cont2)
                    if hist2D_loc_field.GetBinContent(bint)==0.:
                        hist2D_loc_field.SetBinContent(bint, cont4)
            inRoot.Close()
    
    ch.Merge(name)
    outRoot = r.TFile(name, 'update')

    if allHists:
        r.gStyle.SetPadRightMargin(0.5)
        r.gStyle.SetOptStat(0)
        prep2D(hist2D_l_pitch, 'L value', 'pitch [deg]', '#sum{flux}',False)
        prep2D(hist2D_l_pitch_en, 'L value', 'pitch [deg]', '#entries',False)
        prep2D(hist2D_loc, 'longitude', 'latitude', 'time',False)
        prep2D(hist2D_loc_flux, 'longitude', 'latitude', '#LT flux #GT',False)
        prep2D(hist2D_loc_field, 'longitude', 'latitude', 'field [nT]',False)
        
        hist2D_l_pitch.Write()
        hist2D_l_pitch_en.Write()
        hist2D_loc.Write()
        hist2D_loc_flux.Write()
        hist2D_loc_field.Write()
    outRoot.Close()
    return True

# return magnetic field strenght at geomagnetic equator
# ref McIlwain1966 'Magnetic coordinates'
def getBeq(iL):
    M=0.311653 # gauss RE = e-4T
    # translate in nT
    return M/pow(iL,3)*pow(10,5)

# return equatorial pitch angle in degrees
def getAlpha_eq(iAlpha, iB, iBeq):
    
    iAlpha_rad = np.radians(iAlpha)
    if np.sqrt(iBeq/iB) > 1:
        alpha_eq = iAlpha_rad
        print("B={}, Beq={}".format(iB,iBeq))
        print("ATTENTION! iBeq/iB>1. This should not happen!")
        return 999

    else:
        if iAlpha<=90:
            alpha_eq = np.arcsin( np.sin( iAlpha_rad ) * np.sqrt(iBeq/iB) )
        else:
            alpha_eq = np.pi - np.arcsin( np.sin( iAlpha_rad ) * np.sqrt(iBeq/iB) )

        return math.degrees(alpha_eq)
