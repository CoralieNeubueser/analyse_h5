import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required=True)
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp', action='store_true', help='Analyse HEPP data.')
parser.add_argument('--merge', action='store_true', help='Merge all runs.')
parser.add_argument('--ana', action='store_true', help='Analyse all runs.')
parser.add_argument('--test', action='store_true', help='Analyse test runs.')
parser.add_argument('-q','--quiet', action='store_true', help='Run without printouts.')
args,_=parser.parse_known_args()

runs = args.numRuns
os.system('source /opt/exp_software/limadou/set_env_standalone.sh')

datapaths = []
if args.hepd:
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3h5/*.h5')
    if args.test:
        datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_orig/*.h5')
        datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate/*.h5')
        datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_05_95/*.h5')
        datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate_05_95/*.h5')

# select HEPP data from 22-26.02.2019 (solar quiet period)
elif args.hepp:
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190222*.h5')
    datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190223*.h5')
    datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190224*.h5')
    datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190225*.h5')
    datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190226*.h5')

    
if not args.merge and not args.ana:

    if len(datapaths) < runs:
        print("Only {} files available for reading. ".format(len(datapaths)))
        runs = len(datapaths)

    for irun,run in enumerate(datapaths):
        if irun>(runs-1):
            break
        
        if args.test:
            outRootDir = os.path.split(run)[0]
            outfile = home()+"/root/L3_test/"+os.path.split(outRootDir)[1]+'/'+(os.path.split(run)[1]).replace("h5","root")
        else:
            outfile = home()+"/root/"+(os.path.split(run)[1]).replace("h5","root")
    
        
        #print("Test if output exists: ", outfile)
        if os.path.isfile(outfile):
            print("Output root file already exists... \n read in the next file. ")
            runs=runs+1
        else:
            cmd='python3 python/readH5.py --inputFile '+str(run)
            if not args.quiet:
                cmd+=' --debug'
            if args.hepd:
                cmd+=' --data hepd'
            elif args.hepp:
                cmd+=' --data hepp'
            print(cmd)
            os.system(cmd)

# run analysis on single merged root file
elif args.ana and not args.test:
    mge = 'root/all_hepd_'+str(runs)+'runs.root'
    print("Run analysis on single file: ", str(mge))
    if args.hepp:
        mge = 'root/all_hepp_'+str(runs)+'runs.root'

    # open root file
    runana='python3 python/analyseRoot.py --inputFile '+str(mge)
    os.system(runana)

elif args.ana and args.test:
    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:
        mge = home()+"/root/L3_test/"+t+'/all.root'
        mge = 'root/all_hepd_'+str(runs)+'runs.root'
        print("Run analysis on single file: ", str(mge))

        # open root file 
        runana='python3 python/analyseRoot.py --inputFile '+str(mge)
        os.system(runana)

# merge all root files
elif args.merge and not args.test:
    mge=''
    runList=[]

    if args.hepd:
        # write 
        mge = 'root/all_hepd_'+str(runs)+'runs.root'
        runList = glob.glob('root/CSES_HEP_DDD_*.root')

    elif args.hepp: 
        mge = 'root/all_hepp_'+str(runs)+'runs.root'
        runList = glob.glob('root/CSES_01_HEP_1_*.root')

    if len(runList) < runs:
        print("Only {} files available for merge. ".format(len(runList)))
        runs = len(runList)
        if args.hepd:
            mge = 'root/all_hepd_'+str(runs)+'runs.root'
        elif args.hepp:
            mge = 'root/all_hepp_'+str(runs)+'runs.root'

    print("Merge files in: ", mge)
    merge(mge, runList, runs)

elif args.merge and args.test:

    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:    
        mge = home()+"/root/L3_test/"+t+'/all.root'
        runList = glob.glob('root/L3_test/'+t+'/CSES_*.root')
        runs = len(runList)
        if runs>0:
            print("Merge files in: ", mge)
            merge(mge, runList, runs)

