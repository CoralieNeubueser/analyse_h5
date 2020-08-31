import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required=True)
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp', action='store_true', help='Analyse HEPP data.')
parser.add_argument('--merge', action='store_true', help='Merge all runs.')
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
                runpath, runname = os.path.split(run)
                logdir = home() + '/log'
                frunname = 'job_%s.sh'%(str(runname.replace('.h5','')))
                print(frunname)
                frun = None
                try:
                    frun = open(logdir+'/'+frunname, 'w')
                except IOError as e:
                    print("I/O error({0}): {1}".format(e.errno, e.strerror))
                    time.sleep(10)
                    frun = open(logdir+'/'+frunname, 'w')
                print(frun)

                os.system('chmod 777 %s/%s'%(logdir,frunname))
                frun.write('#!/bin/bash\n')
                frun.write('unset LD_LIBRARY_PATH\n')
                frun.write('unset PYTHONHOME\n')
                frun.write('unset PYTHONPATH\n')
                frun.write('export JOBDIR=$PWD\n')
                frun.write('source %s\n' % (path_to_INIT))
                frun.write('cd %s\n'%(home()))
                frun.write(cmd+'\n')
            
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
                fsub.write('RequestCpus = 4\n')        
                fsub.write('+JobFlavour = "espresso"\n')
                #fsub.write('+AccountingGroup = "group_u_FCC.local_gen"\n')
                fsub.write('queue 1\n')
                fsub.close()

                cmdBatch="condor_submit -name sn-01.cr.cnaf.infn.it %s/%s \n"%(logdir,fsubname)

                print(cmdBatch)
                p = subprocess.Popen(cmdBatch, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

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
        mge = 'data/root/all_hepd_'+str(runs)+'runs.root'
        runList = glob.glob('data/root/CSES_HEP_DDD_*.root')

    elif args.hepp: 
        mge = 'data/root/all_hepp_'+str(runs)+'runs.root'
        runList = glob.glob('data/root/CSES_01_HEP_1_*.root')

    if len(runList) < runs:
        print("Only {} files available for merge. ".format(len(runList)))
        runs = len(runList)
        if args.hepd:
            mge = 'data/root/all_hepd_'+str(runs)+'runs.root'
        elif args.hepp:
            mge = 'data/root/all_hepp_'+str(runs)+'runs.root'

    print("Merge files in: ", mge)
    merge(mge, runList, runs, args.allHists)

elif args.merge and args.test:

    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:    
        mge = "data/root/L3_test/"+t+'/all.root'
        runList = glob.glob(sharedOutPath()+'data/root/L3_test/'+t+'/CSES_*.root')
        runs = len(runList)
        if runs>0:
            print("Merge files in: ", mge)
            merge(mge, runList, runs, args.allHists)

