import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required=True)
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp_l', action='store_true', help='Analyse HEPP-L data.')
parser.add_argument('--hepp_h', action='store_true', help='Analyse HEPP-H data.')
parser.add_argument('--merge', action='store_true', help='Merge all runs.')
parser.add_argument('--select', action='store_true', help='Merge all runs.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--day', type=int, required=False, help='Merge orbits of a specific day [yyyymmdd].')
parser.add_argument('--month', type=int, required=False, help='Merge orbits of a specific month [yyyymm].')
parser.add_argument('--allHists', action='store_true', help='Merge all runs, including the histograms.')
parser.add_argument('--ana', action='store_true', help='Analyse all runs.')
parser.add_argument('--test', action='store_true', help='Analyse test runs.')
parser.add_argument('--useVersion', type=str, default='v2.1', help='Define wether flux=0 is stored.')
parser.add_argument('--submit', action='store_true', help='Submit to HTCondor batch farm.')
parser.add_argument('-q','--quiet', action='store_true', help='Run without printouts.')
args,_=parser.parse_known_args()

runs = args.numRuns
os.system('source /opt/exp_software/limadou/set_env_standalone.sh')

det = 'hepd'
datapaths = []
if args.hepd:
    #datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3h5/*.h5')
    # run all orbits in August 2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3h5/*201905*.h5')
    if args.test:
        #datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_orig/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_05_95/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate_05_95/*.h5')
        # new test-sample corrected reconstructed energy 
        datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3_repro/*.h5')
elif args.hepp_l:
    det = 'hepp_l'
    # get HEPP data of quiet period 1.-5.08.2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1_L02*201903*.h5') #('/storage/gpfs_data/limadou/data/flight_data/analysis/data/h5/HEPP_august_2018/*.h5')
    # ('/home/LIMADOU/cneubueser/public/HEPP_august_2018/*.h5')
    
    # select HEPP data from 22-26.02.2019 (solar quiet period) 
    #datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190222*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190223*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190224*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190225*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190226*.h5')
elif args.hepp_h:
    det = 'hepp_h'
    # get HEPP data of quiet period 1.-5.08.2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_2_L02*201903*.h5') #('/storage/gpfs_data/limadou/data/flight_data/analysis/data/h5/HEPP_august_2018/*.h5')                                        
# sort files by time
datapaths.sort(key=os.path.getmtime)

# define the submission file for condor
args_file = open(home()+'/log/arguments.txt', "w")
# run on single semi-orbits
if not args.merge and not args.ana and not args.select:
    
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
        if args.hepp_l or args.hepp_h:
            #            print(OrbitDateTime)
            #            print(OrbitDateTime[5])
            #            print(OrbitDateTime[3])
            duration = abs(int(OrbitDateTime[6]) - int(OrbitDateTime[4]))

        if duration < 3000: ## half-orbit not completed
            print('Not full semi-orbit recorded, but only: '+str(duration)+'[mmss]')
            print('Try next run..')
            continue

        buildPath = sharedOutPath()+"data/root/"+args.useVersion+"/"+det+'/'
        if args.integral:
            buildPath += str(args.integral)+"s/"    
        
        if args.test:
            buildPath = sharedOutPath()+"data/root/"+args.useVersion+"/"
            outRootDir = os.path.split(run)[0]
            outfile = buildPath+"/L3_test/"+os.path.split(outRootDir)[1]+'/'+(os.path.split(run)[1]).replace("h5","root")
        else:
            outfile = buildPath+(os.path.split(run)[1]).replace("h5","root")
        print(outfile)

        # Test if output exists
        if os.path.isfile(outfile):
            print("Output root file already exists... \n read in the next file. ")
            runs=runs+1
        else:
            if not os.path.isdir(buildPath):
                os.makedirs(buildPath)
            cmd='python3 python/readH5.py --inputFile '+str(run)
            if args.useVersion is 'v2.1':
                cmd='python3 python/readH5_v2.1.py --inputFile '+str(run)
            if args.integral:
                cmd+=' --integral '+str(args.integral)
            if not args.quiet:
                cmd+=' --debug'
            if args.hepd:
                cmd+=' --data hepd'
            elif args.hepp_l:
                cmd+=' --data hepp_l'
            elif args.hepp_h:
                cmd+=' --data hepp_h'
            if args.useVersion is not 'v2':
                cmd+=' --useVersion '+args.useVersion
            print(cmd)

            if args.submit:
                exefilename = 'job_%s.sh'%(str(os.path.split(run)[1].replace('.h5','')))
                exefile = writeExecutionFile(home()+'/log/'+exefilename, cmd)
                print("Write execution file to:", exefilename)
                args_file.write("%s\n"%(home()+'/log/'+exefilename))
            else:
                os.system(cmd)
                
    args_file.close()
    if args.submit:
        SubmitListToCondor(args_file)

        
# run analysis on single merged root file
elif args.ana and not args.test:
    mge = sharedOutPath()+'data/root/'+args.useVersion+'/all_hepd_'+str(runs)+'runs.root'
    print("Run analysis on single file: ", str(mge))

    # open root file
    runana='python3 python/analyseRoot.py --inputFile '+str(mge)
    os.system(runana)

elif args.ana and args.test:
    tests=['L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:
        mge = sharedOutPath()+"data/root/"+args.useVersion+"/L3_test/"+t+'/all.root'
        mge = 'data/root/'+args.useVersion+'/all_hepd_'+str(runs)+'runs.root'
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
        mge = sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd.root'
        findOld = None
        runList = []
        if args.day:
            runList = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/CSES_HEP_DDD_*'+str(args.day)+'*.root')
            mge = sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(args.day)+'_'+str(len(runList))+'_runs.root'
            findOld = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(args.day)+'*.root')
            runs = len(runList)
        elif args.month:
            runList = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/CSES_HEP_DDD_*'+str(args.month)+'*.root')
            mge = sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(args.month)+'_'+str(len(runList))+'_runs.root'
            findOld = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(args.month)+'*')
            runs = len(runList)
        else:
            runList = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/CSES_HEP_DDD_*.root')
            runs = len(runList)

    elif args.hepp_l or args.hepp_h: 
        index=1
        if args.hepp_h:
            index=2
        mge = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/all_'+det+'.root'
        runList = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/CSES_01_HEP_'+str(index)+'_L02*.root') #CSES_01_HEP_1_*.root')
        if args.day:
            runList = sorted( glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/CSES_01_HEP_'+str(index)+'_L02_*'+str(args.day)+'*.root'), key=lambda x:float(x[-46:-41]) )
            mge = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/all_'+det+'_'+str(args.day)+'_'+str(len(runList))+'_runs.root'
            findOld = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/all_'+det+'_'+str(args.day)+'*.root')
            runs = len(runList)

    #print(runs)
    oldruns=0
    #print(findOld)
    if len(findOld)>0:
        for old in findOld:
            head, tail = os.path.split(old)
            oldruns_tmp=list(map(int, re.findall(r'\d+', tail)))[1]
            if oldruns_tmp<runs:
                os.system('rm {}'.format(old))
            if oldruns_tmp>oldruns:
                oldruns=oldruns_tmp

    # merge files only if not already exists and existing file has less inputs  
    #if runs>0 and oldruns<runs:
    print("Merge files in: ", mge)
    merge(mge, runList, runs, args.allHists)

    #cmd  = 'hadd -f -k '+mge+' '
    #for ifile in runList:
    #    cmd += ifile+' '
    #cmd += '\n'
    cmd = 'python3 python/writeDayAverages.py --drawHistos --inputFile '+mge
    if args.hepd:
        cmd += ' --data hepd '
    elif args.hepp_l:
        cmd += ' --data hepp_l '
    elif args.hepp_h:
        cmd += ' --data hepp_h '

    if args.day:
        cmd += '--day '+str(args.day)
        if args.submit:
            SubmitToCondor(cmd, mge, 1)
        else:
            os.system(cmd)
    
        

elif args.merge and args.test:

    tests=['L3_repro'] #'L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:    
        mge = sharedOutPath()+"data/root/"+args.useVersion+"/L3_test/"+t+'/all.root'
        if args.day:
            mge = sharedOutPath()+"data/root/"+args.useVersion+"/L3_test/"+t+'/all_'+str(args.day)+'.root'
            runList = sorted( glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/L3_test/'+t+'/CSES_*'+str(args.day)+'*.root'), key=lambda x:float(x[-10:-4]))
        #runList = sorted( glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/L3_test/'+t+'/CSES_*.root'), key=lambda x:float(x[-10:-4]))
        runs = len(runList)
        #print(runList)
        if runs>0:
            print("Merge files in: ", mge)
            merge(mge, runList, runs, args.allHists)
            #cmd  = 'hadd -f -k '+mge+' '
            #for ifile in runList:
            #    cmd += ifile+' '
            #    cmd += '\n'
 
        cmd = 'python3 python/writeDayAverages.py --drawHistos --inputFile '+mge
        if args.hepd:
            cmd += ' --data hepd '
        if args.day:
            cmd += '--day '+str(args.day)
        if args.submit:
            SubmitToCondor(cmd, mge, 1)
        else:
            os.system(cmd)

elif args.select:

    findFile = []
    detPath = det
    fileSnip = 'all_'+det+'_'
    if args.test:
        detPath = 'L3_test/L3_repro'
        fileSnip = 'all_'
    if args.month:
        for iday in range(1,32):
            strday = str(iday)
            if iday < 10:
                strday = '0'+str(iday)
            print( sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+strday+'*.root' )
            found = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+strday+'*.root')
            #print(found)
            for every in found:
                findFile.append( every )
    elif args.day:
        findFile = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.day)+'*.root')

    #print(findFile)
    for ind,foundFile in enumerate(findFile):
        ind+=1
        strDay = str(ind)
        if ind<10:
            strDay = '0'+str(ind)
        dayint = ind
        cmd = 'python3 python/findHighFluxes.py --inputFile '+foundFile+' --data '+det
        if args.day:
            cmd += ' --day '+str(args.day)
            dayint = args.day
        elif args.month:
            cmd += ' --day '+str(args.month)+strDay
            dayint = int(str(args.month)+strDay)

        if args.submit:
            exefilename = 'job_%s.sh'%(str(fileSnip+str(dayint)+'.root'))
            exefile = writeExecutionFile(home()+'/log/'+exefilename, cmd)
            print("Write execution file to:", exefilename)
            args_file.write("%s\n"%(home()+'/log/'+exefilename))
        else:
            os.system(cmd)

    args_file.close()
    if args.submit:
        SubmitListToCondor(args_file)

#    if args.submit:
#        exefilename = 'job_%s.sh'%(str('all_'+str(args.day)+'.root'))
#        exefile = writeExecutionFile(home()+'/log/'+exefilename, cmd)
#        print("Write execution file to:", exefilename)
#        args_file.write("%s\n"%(home()+'/log/'+exefilename))
#        args_file.close()
#        SubmitListToCondor(args_file)
#    else:
#        os.system(cmd)
