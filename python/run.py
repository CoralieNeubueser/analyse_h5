import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required='--merge' not in sys.argv)
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp', action='store_true', help='Analyse HEPP data.')
parser.add_argument('--merge', action='store_true', help='Merge all runs.')
parser.add_argument('-q','--quiet', action='store_true', help='Run without printouts.')
args,_=parser.parse_known_args()

runs = args.numRuns
os.system('source /opt/exp_software/limadou/set_env_standalone.sh')

datapath = ''
if args.hepd:
    datapath = '/storage/gpfs_data/limadou/data/flight_data/L3h5/'
elif args.hepp:
    datapath = '/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*'
    
if not args.merge:
    for irun,run in enumerate(glob.glob(datapath+'*.h5')):
        if irun>(runs-1):
            break
        outfile = home()+"/root/"+(os.path.split(run)[1]).replace("h5","root")

        print("Test if output exists: ", outfile)
        if os.path.isfile(outfile):
            print("Output root file already exists... \n run analysis on the next run. ")
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

else:
    mge=''
    runList=[]

    if args.hepd:
        # write 
        mge = 'root/all_hepd.root'
        runList = glob.glob('root/CSES_HEP_DDD_*.root')
    elif args.hepp: 
        mge = 'root/all_hepp.root'
        runList = glob.glob('root/CSES_01_HEP_1_*.root')

    print("Merge all files in: ", mge)
    merge(mge, runList)
    
