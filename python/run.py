import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required=True)
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp', action='store_true', help='Analyse HEPP data.')
parser.add_argument('--merge', action='store_true', help='Merge all runs.')
parser.add_argument('--day', type=int, required=False, help='Merge orbits of a specific day [yyyymmdd].')
parser.add_argument('--month', type=int, required=False, help='Merge orbits of a specific month [yyyymm].')
parser.add_argument('--allHists', action='store_true', help='Merge all runs, including the histograms.')
parser.add_argument('--ana', action='store_true', help='Analyse all runs.')
parser.add_argument('--test', action='store_true', help='Analyse test runs.')
parser.add_argument('--submit', action='store_true', help='Submit to HTCondor batch farm.')
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

elif args.hepp:

    # get HEPP data of quiet period 1.-5.08.2018
    datapaths = glob.glob('/home/LIMADOU/cneubueser/public/HEPP_august_2018/*.h5')
    
    # select HEPP data from 22-26.02.2019 (solar quiet period) 
    #datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190222*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190223*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190224*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190225*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190226*.h5')

# run on single semi-orbits
if not args.merge and not args.ana:

    if len(datapaths) < runs:
        print("Only {} files available for reading. ".format(len(datapaths)))
        runs = len(datapaths)

    for irun,run in enumerate(datapaths):
        if irun>(runs-1):
            break

        # test if half-oribit is complete
        # find int in string
        OrbitDateTime = re.findall('\d+', run)
        # calculate time between start and stop of orbit
        duration = abs(int(OrbitDateTime[7]) - int(OrbitDateTime[5]))
        if args.hepp:
            duration = abs(int(OrbitDateTime[5]) - int(OrbitDateTime[3]))

        if duration < 3000: ## half-orbit not completed
            print('Not full semi-orbit recorded, but only: '+str(duration)+'[mmss]')
            print('Try next run..')
            continue

        if args.test:
            outRootDir = os.path.split(run)[0]
            outfile = sharedOutPath()+"data/root/L3_test/"+os.path.split(outRootDir)[1]+'/'+(os.path.split(run)[1]).replace("h5","root")
        else:
            outfile = sharedOutPath()+"data/root/"+(os.path.split(run)[1]).replace("h5","root")
            print(outfile)

        # Test if output exists
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

            if args.submit:
                SubmitToCondor(cmd, run, irun)
            else:
                os.system(cmd)

# run analysis on single merged root file
elif args.ana and not args.test:
    mge = sharedOutPath()+'data/root/all_hepd_'+str(runs)+'runs.root'
    print("Run analysis on single file: ", str(mge))
    if args.hepp:
        mge = 'data/root/all_hepp_'+str(runs)+'runs.root'

    # open root file
    runana='python3 python/analyseRoot.py --inputFile '+str(mge)
    os.system(runana)

elif args.ana and args.test:
    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:
        mge = sharedOutPath()+"data/root/L3_test/"+t+'/all.root'
        mge = 'data/root/all_hepd_'+str(runs)+'runs.root'
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
        mge = sharedOutPath()+'data/root/hepd/all_hepd.root'
        findOld = None
        runList = []
        if args.day:
            runList = glob.glob(sharedOutPath()+'data/root/CSES_HEP_DDD_*'+str(args.day)+'*.root')
            mge = sharedOutPath()+'data/root/hepd/all_hepd_'+str(args.day)+'_'+str(len(runList))+'_runs.root'
            findOld = glob.glob(sharedOutPath()+'data/root/hepd/all_hepd_'+str(args.day)+'*.root')
            runs = len(runList)
        elif args.month:
            runList = glob.glob(sharedOutPath()+'data/root/CSES_HEP_DDD_*'+str(args.month)+'*.root')
            mge = sharedOutPath()+'data/root/hepd/all_hepd_'+str(args.month)+'_'+str(len(runList))+'_runs.root'
            findOld = glob.glob(sharedOutPath()+'data/root/hepd/all_hepd_'+str(args.month)+'*')
            runs = len(runList)
        else:
            runList = glob.glob(sharedOutPath()+'data/root/CSES_HEP_DDD_*.root')
            runs = len(runList)
    elif args.hepp: 
        mge = sharedOutPath()+'data/root/hepp/all_hepp.root'
        runList = glob.glob(sharedOutPath()+'data/root/CSES_01_HEP_1_*.root')
    
    oldruns=0
    print(findOld)
    if len(findOld)==1:
        head, tail = os.path.split(findOld[0])
        oldruns=list(map(int, re.findall(r'\d+', tail)))[1]
    elif len(findOld)>1:
        for old in findOld:
            os.system('rm {}'.format(old))

    # merge files only if not already exists and existing file has less inputs  
    if runs>0 and oldruns<runs:
        print("Merge files in: ", mge)
        merge(mge, runList, runs, args.allHists)

        if args.hepd:
            cmd = 'python3 python/writeDayAverages.py --inputFile '+mge+' --data hepd --fit '
            if args.day:
                cmd += '--day '+str(args.day)
            if args.submit:
                SubmitToCondor(cmd, mge, 1)
            else:
                os.system(cmd)
    
elif args.merge and args.test:

    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:    
        mge = "data/root/L3_test/"+t+'/all.root'
        runList = glob.glob(sharedOutPath()+'data/root/L3_test/'+t+'/CSES_*.root')
        runs = len(runList)
        if runs>0:
            print("Merge files in: ", mge)
            merge(mge, runList, runs, args.allHists)

