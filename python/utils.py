import os,sys,math
import subprocess
import time
import numpy as np
import ROOT as r
import h5py
from collections import defaultdict
from drawFunctions import *
import pandas as pd


path_to_INIT = '/opt/exp_software/limadou/set_env_standalone.sh'

def home():
    return os.getcwd()

# define path to shared directory for output file storage
def sharedOutPath():
    return '/storage/gpfs_data/limadou/data/flight_data/analysis/'

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

# define corser L value bins 
# returns the number of bins and a list of bins
def getCorserLbins():
    l_x_bins = []
    for x in range(0,10):
        l_x_bins.append(1.+float(x))
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

# get time bin width for all detectors
# returns the time bin width in seconds
def getTimeBins(det):
    if det=='noaa':
        return 2
    else:
        return 1

# define energy bins
# returns 
# 1. the number of energy bins for either hepd of hepp data
# 2. a list of lower edge energy values
# 3. upper bin edge of last energy bin
def getEnergyBins(det, rebinned):
    if det=='hepd':
        # original upper bin edges: [[ 4. 8.96811 10.9642  13.4045  16.388   20.0355  24.4949  29.9468  36.6122  44.7611  54.7237  66.9037 ]]
        return 12, [2.0, 6.5, 10.0, 12.2, 14.9, 18.2, 22.3, 27.2, 33.3, 40.7, 49.7, 60.8], 66.9
    elif det=='noaa':
        return 4, [0.04, 0.13, 0.287, 0.612], 0.612
    elif det=='hepp_l' and rebinned:
        return 16, [0.1, 0.28125, 0.4625, 0.64375, 0.825, 1.00625, 1.1875, 1.36875, 1.55, 1.73125, 1.9125, 2.09375, 2.275, 2.45625, 2.6375, 2.81875], 3
    elif det=='hepp_l' and not rebinned:
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
                     2.9658823 , 2.97725487, 2.98862743, 3. ], 3.
    elif det=='hepp_h' and rebinned:
        return 16, [ 1.0, 4.388239860534668, 7.776480197906494, 11.164719581604004, 14.552960395812988, 17.941200256347656, 21.329439163208008, 24.717679977416992, 28.105920791625977, 31.494159698486328, 34.88240051269531, 38.2706413269043, 41.658878326416016, 45.047119140625, 48.435359954833984, 51.82360076904297 ], 55.0
    
    elif det=='hepp_h' and not rebinned:
        return 256, [1.0, 
                     1.2117650508880615, 
                     1.4235299825668335, 
                     1.635295033454895, 
                     1.847059965133667, 
                     2.0588250160217285, 
                     2.27059006690979, 
                     2.4823551177978516, 
                     2.694119930267334, 
                     2.9058849811553955, 
                     3.117650032043457, 
                     3.3294150829315186, 
                     3.541179895401001, 
                     3.7529449462890625, 
                     3.964709997177124, 
                     4.1764750480651855, 
                     4.388239860534668, 
                     4.600005149841309, 
                     4.811769962310791, 
                     5.023534774780273, 
                     5.235300064086914, 
                     5.4470648765563965, 
                     5.658830165863037, 
                     5.8705949783325195, 
                     6.082359790802002, 
                     6.294125080108643, 
                     6.505889892578125, 
                     6.717655181884766, 
                     6.929419994354248, 
                     7.1411848068237305, 
                     7.352950096130371, 
                     7.5647149085998535, 
                     7.776480197906494, 
                     7.988245010375977, 
                     8.200010299682617, 
                     8.411774635314941, 
                     8.623539924621582, 
                     8.835305213928223, 
                     9.047069549560547, 
                     9.258834838867188, 
                     9.470600128173828, 
                     9.682365417480469, 
                     9.894129753112793, 
                     10.105895042419434, 
                     10.317660331726074, 
                     10.529424667358398, 
                     10.741189956665039, 
                     10.95295524597168, 
                     11.164719581604004, 
                     11.376484870910645, 
                     11.588250160217285, 
                     11.800015449523926, 
                     12.01177978515625, 
                     12.22354507446289, 
                     12.435310363769531, 
                     12.647074699401855, 
                     12.858839988708496, 
                     13.070605278015137, 
                     13.282369613647461, 
                     13.494134902954102, 
                     13.705900192260742, 
                     13.917664527893066, 
                     14.129429817199707, 
                     14.341195106506348, 
                     14.552960395812988, 
                     14.764724731445312, 
                     14.976490020751953, 
                     15.188255310058594, 
                     15.400019645690918, 
                     15.611784934997559, 
                     15.8235502243042, 
                     16.035314559936523, 
                     16.247079849243164, 
                     16.458845138549805, 
                     16.670610427856445, 
                     16.882375717163086, 
                     17.094139099121094, 
                     17.305904388427734, 
                     17.517669677734375, 
                     17.729434967041016, 
                     17.941200256347656, 
                     18.152965545654297, 
                     18.364730834960938, 
                     18.576494216918945, 
                     18.788259506225586, 
                     19.000024795532227, 
                     19.211790084838867, 
                     19.423555374145508, 
                     19.63532066345215, 
                     19.84708595275879, 
                     20.058849334716797, 
                     20.270614624023438, 
                     20.482379913330078, 
                     20.69414520263672, 
                     20.90591049194336, 
                     21.11767578125, 
                     21.329439163208008, 
                     21.54120445251465, 
                     21.75296974182129, 
                     21.96473503112793, 
                     22.17650032043457, 
                     22.38826560974121, 
                     22.60003089904785, 
                     22.81179428100586, 
                     23.0235595703125, 
                     23.23532485961914, 
                     23.44709014892578, 
                     23.658855438232422, 
                     23.870620727539062, 
                     24.08238410949707, 
                     24.29414939880371, 
                     24.50591468811035, 
                     24.717679977416992, 
                     24.929445266723633, 
                     25.141210556030273, 
                     25.352975845336914, 
                     25.564739227294922, 
                     25.776504516601562, 
                     25.988269805908203, 
                     26.200035095214844, 
                     26.411800384521484, 
                     26.623565673828125, 
                     26.835329055786133, 
                     27.047094345092773, 
                     27.258859634399414, 
                     27.470624923706055, 
                     27.682390213012695, 
                     27.894155502319336, 
                     28.105920791625977, 
                     28.317684173583984, 
                     28.529449462890625, 
                     28.741214752197266, 
                     28.952980041503906, 
                     29.164745330810547, 
                     29.376510620117188, 
                     29.588275909423828, 
                     29.800039291381836, 
                     30.011804580688477, 
                     30.223569869995117, 
                     30.435335159301758, 
                     30.6471004486084, 
                     30.85886573791504, 
                     31.070629119873047, 
                     31.282394409179688, 
                     31.494159698486328, 
                     31.70592498779297, 
                     31.91769027709961, 
                     32.12945556640625, 
                     32.34122085571289, 
                     32.55298614501953, 
                     32.76475143432617, 
                     32.97651672363281, 
                     33.18827819824219, 
                     33.40004348754883, 
                     33.61180877685547, 
                     33.82357406616211, 
                     34.03533935546875, 
                     34.24710464477539, 
                     34.45886993408203, 
                     34.67063522338867, 
                     34.88240051269531, 
                     35.09416580200195, 
                     35.305931091308594, 
                     35.517696380615234, 
                     35.729461669921875, 
                     35.94122314453125, 
                     36.15298843383789, 
                     36.36475372314453, 
                     36.57651901245117, 
                     36.78828430175781, 
                     37.00004959106445, 
                     37.211814880371094, 
                     37.423580169677734, 
                     37.635345458984375, 
                     37.847110748291016, 
                     38.058876037597656, 
                     38.2706413269043, 
                     38.48240661621094, 
                     38.69417190551758, 
                     38.90593338012695, 
                     39.117698669433594, 
                     39.329463958740234, 
                     39.541229248046875, 
                     39.752994537353516, 
                     39.964759826660156, 
                     40.1765251159668, 
                     40.38829040527344, 
                     40.60005569458008, 
                     40.81182098388672, 
                     41.02358627319336, 
                     41.2353515625, 
                     41.44711685180664, 
                     41.658878326416016, 
                     41.870643615722656, 
                     42.0824089050293, 
                     42.29417419433594, 
                     42.50593948364258, 
                     42.71770477294922, 
                     42.92947006225586, 
                     43.1412353515625, 
                     43.35300064086914, 
                     43.56476593017578, 
                     43.77653121948242, 
                     43.98829650878906, 
                     44.2000617980957, 
                     44.41182327270508, 
                     44.62358856201172, 
                     44.83535385131836, 
                     45.047119140625, 
                     45.25888442993164, 
                     45.47064971923828, 
                     45.68241500854492, 
                     45.89418029785156, 
                     46.1059455871582, 
                     46.317710876464844, 
                     46.529476165771484, 
                     46.741241455078125, 
                     46.953006744384766, 
                     47.16476821899414, 
                     47.37653350830078, 
                     47.58829879760742, 
                     47.80006408691406, 
                     48.0118293762207, 
                     48.223594665527344, 
                     48.435359954833984, 
                     48.647125244140625, 
                     48.858890533447266, 
                     49.070655822753906, 
                     49.28242111206055, 
                     49.49418640136719, 
                     49.70595169067383, 
                     49.9177131652832, 
                     50.129478454589844, 
                     50.341243743896484, 
                     50.553009033203125, 
                     50.764774322509766, 
                     50.976539611816406, 
                     51.18830490112305, 
                     51.40007019042969, 
                     51.61183547973633, 
                     51.82360076904297, 
                     52.03536605834961, 
                     52.24713134765625, 
                     52.45889663696289, 
                     52.670658111572266, 
                     52.882423400878906, 
                     53.09418869018555, 
                     53.30595397949219, 
                     53.51771926879883, 
                     53.72948455810547, 
                     53.94124984741211, 
                     54.15301513671875, 
                     54.36478042602539, 
                     54.57654571533203, 
                     54.78831100463867,
                     55.0], 55.0

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
def getGeomFactor(det,energyBin):
    if det=='hepd':
        # geometrical factors
        ele_corr_GF = [ 131.128, 545.639, 560.297, 530.937, 477.827, 413.133, 334.176, 252.3, 204.52, 103.216, 77.5552, 61.1536 ]
        #ele_corr_GF = 12*[0.]
    elif det=='hepp_l_channel_narrow':
        ele_corr_GF = 16*[0.12] # 1. The geometrical factors of HEPP-L is 0.12 cm^2 sr for 5 detectors and 0.73 cm^2 sr for 4 detectors. (Zhenxia priv. comm: 12 May 2021)
    elif det=='hepp_l_channel_wide':
        ele_corr_GF = 16*[0.73] # 1. The geometrical factors of HEPP-L is 0.12 cm^2 sr for 5 detectors and 0.73 cm^2 sr for 4 detectors. (Zhenxia priv. comm: 12 May 2021)
    # the geometrical factor is mixed between the channels.. use the smallest.
    elif det=='hepp_l_channel_all':
        ele_corr_GF = 16*[0.12]
    # unknown geometrical factors.. for HEPP-H and NOAA-POES19
    elif det=='hepp_h':
        ele_corr_GF = 16*[1.] 
    elif det=='noaa':
        # taken from https://www.ngdc.noaa.gov/stp/satellite/poes/docs/NGDC/MEPED%20telescope%20processing%20ATBD_V1.pdf
        # found 17.01.2022
        #ele_corr_GF = [100./1.24, 100./1.44, 100./0.75, 100./0.55]
        ele_corr_GF = [1.24, 1.44, 0.75, 0.55]
    return ele_corr_GF[energyBin]

def getInverseGeomFactor(det,energyBin):
    ele_corr_GF = getGeomFactor(det,energyBin)
    if det=='hepd':
        # reverse correction, and use original factor
        ele_GF = [ 0.76, 188.26, 326.64, 339.65, 344.99, 331.83, 304.73, 263.56, 217.33, 169.48, 117.31, 71.45 ]
        return getGeomCorr(det,energyBin) / ele_GF[energyBin]
    else:
        return 1./ele_corr_GF

# returns the energy bin width, used to normalise fluxes to counts/s/cm2/sr/MeV 
def getEnergyBinWidth(det, energyBin):
    enBorder=np.array([0.5, 4., 8.96811, 10.9642,  13.4045,  16.388,   20.0355,  24.4949,  29.9468,  36.6122,  44.7611,  54.7237,  66.9037 ])
    # HEPP-L fluxes need normalisation to 11keV
    if det=='hepp_l':
        #enBorder=np.array([0.1, 0.28125, 0.4625, 0.64375, 0.825, 1.00625, 1.1875, 1.36875, 1.55, 1.73125, 1.9125, 2.09375, 2.275, 2.45625, 2.6375, 2.81875, 3.0])
        #enBorder=np.array(16*[0.011])
        return 0.01 
    elif det=='hepp_h':
        return 0.21
    else:
        return enBorder[energyBin+1]-enBorder[energyBin]

# return number and maximum of flux bins to be used in determining the averages and RMS99
def getFluxBins(det):
    if det=='hepd':
        return 50, 0.01
    elif det=='hepp_l':
        return 1000, 10
    else:
        return 1000, 0.05

def getCountsBins(det):
    if det=='hepd':
        return 50
    elif det=='hepp_l':
        return int(1e6)
    #elif det=='noaa':
    #    return int(1e6)
    else:
        return -1

# return a list of days
def getDays(tree):
    lst_days = set()
    # loop through ev to fill list of days
    for ev in tree:
        lst_days.add(ev.day)
    return lst_days

# read in daily averages, rms99, and rms99_of_99
# return list of dict
def readAverageFile(fileName,storedEn,method='rms99'):

    av_Lalpha = [ {} for en in storedEn]
    file = open(fileName, "r")
    next(file)
    for line in file:
        columns = [float(i) for i in line.split()]
        col_energy = columns[0]
        energyStored = col_energy
        en_index = storedEn.index( energyStored )

        if method=='gauss':
            # energy L pitch entries mean meanErr rms rmsErr mean_counts meanErr_counts mpv sigma chi2 avGeomIndex
            # fill dictionary from (L, alpha) -> (mean, rms, mpv, mpvErr, sigma, sigaErr) 
            av_Lalpha[en_index].update( {( columns[1],int(columns[2]) ):( columns[4],columns[6],columns[10],columns[11],columns[12],columns[13],columns[14] )} )
        else:
            # fill dictionary from (L, alpha) -> (mean, rms, rms99, rms99Err, rms99_of_99, rmsErr_99_of_99, weight)
            av_Lalpha[en_index].update( {( columns[1],int(columns[2]) ):( columns[4],columns[6],columns[10],columns[11],columns[12],columns[13],columns[14] )} )
    return av_Lalpha

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


# earth radius at equator in km
RE = 6378.137

# calculate L in dipole approximation
def calculateL(geomLat, r):
    return r / pow(np.cos( np.radians(geomLat) ),2) 

# return magnetic field strenght at geomagnetic equator
# ref McIlwain1966 'Magnetic coordinates'
def getBeq(iL):
    M=0.311653 # gauss RE = e-4T
    # translate in nT
    return M/pow(iL,3)*pow(10,5)

# return equatorial pitch angle in degrees
def getAlpha_eq(iAlpha, iB, iBeq):
    if math.isnan(iAlpha):
        print("Equatorial pitch angle not calculatable.. ")
        print("Pitch angle is nan.")
        return 999
    
    iAlpha_rad = np.radians(iAlpha)
    if iBeq/iB > 1.0:
        alpha_eq = iAlpha_rad
        print("B={}, Beq={}".format(iB,iBeq))
        print("ATTENTION! iBeq/iB>1. This should not happen!")
        return 999

    else:
        #print(iBeq, iB)
        if iAlpha<=90:
            alpha_eq = np.arcsin( np.sin( iAlpha_rad ) * np.sqrt(iBeq/iB) )
        else:
            alpha_eq = np.pi - np.arcsin( np.sin( iAlpha_rad ) * np.sqrt(iBeq/iB) )
        return math.degrees(alpha_eq)

#__________________________________________________________
def getCommandOutput(command):
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}

#__________________________________________________________  
def writeExecutionFile(filename, command):
    frun = None
    try:
        frun = open(filename, 'w')
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        time.sleep(10)
        frun = open(filename, 'w')
        print(frun)
    frun.write('#!/bin/bash\n')
    frun.write('unset LD_LIBRARY_PATH\n')
    frun.write('unset PYTHONHOME\n')
    frun.write('unset PYTHONPATH\n')
    frun.write('export JOBDIR=$PWD\n')
    frun.write('source %s\n' % (path_to_INIT))
    frun.write('cd %s\n'%(home()))
    frun.write('source '+home()+'/env.sh\n')
    frun.write(command+'\n')
    os.system('chmod 777 %s'%(filename))
    
    return frun

#__________________________________________________________
def SubmitToCondor(cmd,run,irun):
    runpath, runname = os.path.split(run)
    logdir = home() + '/log'
    frunname = 'job_%s.sh'%(str(runname.replace('.h5','')))
    print(frunname)

    frun = writeExecutionFile(logdir+'/'+frunname, cmd)

    os.system('chmod 777 %s/%s'%(logdir,frunname))

    os.system("mkdir -p %s/out"%logdir)
    os.system("mkdir -p %s/log"%logdir)
    os.system("mkdir -p %s/err"%logdir)
    
    # create also .sub file here                                                                                                                                                                                   
    fsubname = frunname.replace('.sh','.sub')
    fsub = None
    try:
        fsub = open(logdir+'/'+fsubname, 'w')
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        time.sleep(10)
        fsub = open(logdir+'/'+fsubname, 'w')
    
    fsub.write('executable            = %s/%s\n' %(logdir,frunname))
    fsub.write('arguments             = $(ClusterID) $(ProcId)\n')
    fsub.write('output                = %s/out/job.%s.$(ClusterId).$(ProcId).out\n'%(logdir,str(irun)))
    fsub.write('log                   = %s/log/job.%s.$(ClusterId).log\n'%(logdir,str(irun)))
    fsub.write('error                 = %s/err/job.%s.$(ClusterId).$(ProcId).err\n'%(logdir,str(irun)))
    #fsub.write('RequestCpus = 4\n')
    fsub.write('+JobFlavour = "longlunch"\n')
    fsub.write('queue 1\n')
    fsub.close()
    
    cmdBatch="condor_submit -name sn-02 %s/%s \n"%(logdir,fsubname)
    print(cmdBatch)
    p = subprocess.Popen(cmdBatch, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)


#__________________________________________________________
# submit a list                                                                                                             
def SubmitListToCondor(args, irun=1):

    os.system("mkdir -p %s/out"%(home()+'/log/'))
    os.system("mkdir -p %s/log"%(home()+'/log/'))
    os.system("mkdir -p %s/err"%(home()+'/log/'))

    logdir = home()+'/log/'
    # create also .sub file here 
    fsubname = 'condor.sub'
    fsub = None
    try:
        fsub = open(logdir+'/'+fsubname, 'w')
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        time.sleep(10)
        fsub = open(logdir+'/'+fsubname, 'w')

    fsub.write('arguments             = $(ClusterID) $(ProcId)\n')
    fsub.write('output                = %s/out/job.%s.$(ClusterId).$(ProcId).out\n'%(logdir,str(irun)))
    fsub.write('log                   = %s/log/job.%s.$(ClusterId).log\n'%(logdir,str(irun)))
    fsub.write('error                 = %s/err/job.%s.$(ClusterId).$(ProcId).err\n'%(logdir,str(irun)))
    #fsub.write('RequestCpus = 8\n')
    #fsub.write('Request_Memory = 32 Mb')
    fsub.write('+JobFlavour = "longlunch"\n')
    fsub.write("queue executable,seed from %s" % os.path.join(logdir, "arguments.txt"))    
    fsub.close()

    cmdBatch="condor_submit -name sn-02 %s/%s \n"%(logdir,fsubname)
    print(cmdBatch)
    p = subprocess.Popen(cmdBatch, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)


## read in geomagnetic indices
def readGeomIndex():
    ###
    dataPath = sharedOutPath()+"/data/geomIndices/BDR_SymH_AE_Indices_merged2018and2019_FlagV4b_All.txt"
    if open(dataPath, "r"):
        print('Reading IGRF12 magnetic field values from: ', dataPath)
    else:
        print('ERROR.. could not find txt files for magnetic field values..')
    file = open(dataPath, "r")
    next(file)
    lines = file.readlines()[20:]
    newDict = {}
    for line in lines:
        c, d, t, i1, i2, flag = line.strip().split(' ')
        pair = (int(d.replace('-','')), int(t[:2]), int(t[3:5]))
        newDict[pair] = int(flag)

    return newDict

## read txt files with magetic field strenght for 500km altitude, from IGRF12 model
## day in yyyymmdd
def readIGRF(day):
    dataPath = sharedOutPath()+"/data/geomField/"+str(day)+'_field.txt'
    file = open(dataPath,'r')
    next(file)
    lines = file.readlines()
    newDict={}
    for line in lines:
        d, lat, lon, mag = line.strip().split(' ')
        pair = (int(d.replace('-','')), int(lat), int(lon))
        newDict[pair] = float(mag)
    return newDict

def getGeomIndex(dic, day):
    ###
    print(str(day)[0:4])
    if int(str(day)[0:4])>2020:
        return -1
    else:
        allGeomIndices = []
        for line in dic:
            if line[0] == day:
                allGeomIndices.append(dic[line])
    
        return np.mean(np.array(allGeomIndices))

def getAlphaLindex(alpha_v, L_v):
    # find corresponding L/alpha bin, and use the average Lshell and alpha values                                                  
    Lbin=-1
    Albin=-1
    for il in range(getLbins()[0]):
        if getLbins()[1][il] <= L_v and getLbins()[1][il+1] > L_v:
            Lbin=il
            break
    for ia in range(getPitchBins()[0]):
        if getPitchBins()[1][ia] <= alpha_v and getPitchBins()[1][ia+1] > alpha_v:
            Albin=ia
            break
    return ia,il

## read-in of list of days with EQs
def readEQfile(magMin=0):
    df = pd.read_csv('{}/data/earthquakes_2019_2021.csv'.format(sharedOutPath()), usecols=['time','latitude','longitude','mag'], sep=';')
    df_sel = df[df.mag>magMin]
    months = np.unique([int(m[:7].replace('-','')) for m in df_sel.time])
    days = [int(day[:10].replace('-','')) for day in df_sel.time]
    daytimes = [int(t[11:13])+int(t[14:16])/60 for t in df_sel.time]
    return months, days, daytimes 

## get SAA cut dependent on detector
def getSAAcut(det):
    if det=='noaa':
        return 22000
    else:
        return 25000
