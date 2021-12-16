import os, sys, argparse, re, glob
from utils import *

parser = argparse.ArgumentParser()
### DEFINE WHICH DETECTOR TO USE
parser.add_argument('--hepd', action='store_true', help='Analyse HEPD data.')
parser.add_argument('--hepp_l', action='store_true', help='Analyse HEPP-L data.')
parser.add_argument('--hepp_h', action='store_true', help='Analyse HEPP-H data.')
parser.add_argument('--noaa', action='store_true', help='Analyse NOAA data.')
parser.add_argument('--channel', type=str, required= '--hepp_l' in sys.argv,  choices=['narrow','wide'], help='Choose narrow or wide channels to be read.')
### ACTIONS
### defaults action reads h5 and writes out root files
parser.add_argument('--numRuns', type=int, help='Define number of runs to be analysed.', required=True)
parser.add_argument('--EQlist', action='store_true', help='Run analysis on specific days, defined by list of EQs (2019-2021).')
parser.add_argument('--clean', action='store_true', help='Clean up HEPP files of overlapping orbits.')
parser.add_argument('--merge', action='store_true', help='Merge all runs, per month and writes out RMS99 threshold.')
parser.add_argument('--select', action='store_true', help='Select 1% highest fluxes, writes out new root file.')
parser.add_argument('--ana', action='store_true', help='Analyse all runs.')
parser.add_argument('--cluster', action='store_true', help='Run clustering on 1% highest fluxes.')
parser.add_argument('--cut', type=str, default='99perc', choices=['99perc','weights','z_score_more2','z_score_more3','dummy_cut'], help='Define the cut type.')
parser.add_argument('--window',type=int, default=10, help='Define window length in s.')
parser.add_argument('--seeds', type=int, default=4, help='Define numer of seeds necessary to build cluster.')
parser.add_argument('--onlyIn', nargs='+', choices=['L','Pitch','Energy'], help='Define the parameter space in which to cluster. KEEP THE ORDER')
parser.add_argument('--doubleSeed', action='store_true', help='Require two subsequent seeds for cluster building.')

### OPTIONS
parser.add_argument('--originalE', action='store_true', help='Use fine energy binning.')
parser.add_argument('--fit', action='store_true', help='Fit flux distributions with exponential.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--draw', action='store_true', help='Allows to write out root file with distributions that were used to extract averages, and RMS99.')
parser.add_argument('--allHists', action='store_true', help='Merge all runs, including the histograms.')
parser.add_argument('--test', action='store_true', help='Analyse test runs.')

### define specific period to run over
parser.add_argument('--day', type=int, nargs='+', required=False, help='Merge orbits of a specific day [yyyymmdd].')
parser.add_argument('--month', type=int, required=False, help='Merge orbits of a specific month [yyyymm].')

parser.add_argument('--useVersion', type=str, default='v2', choices=['v1','v2','v2.1','v3'], help='Define wether flux=0 is stored.')
parser.add_argument('--submit', action='store_true', help='Submit to HTCondor batch farm.')
parser.add_argument('-q','--quiet', action='store_true', help='Run without printouts.')
args,_=parser.parse_known_args()

runs = args.numRuns
os.system('source /opt/exp_software/limadou/set_env_standalone.sh')

if args.EQlist:
    readEQlist()

det = 'hepd'
data = 'hepd'
datapaths = []
if args.hepd:
    # run all orbits in August 2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3h5/*'+str(args.month)+'*.h5')
    if args.day:
        for d in args.day:
            datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3h5/*'+str(d)+'*.h5')
    if args.test:
        #datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_orig/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_05_95/*.h5')
        #datapaths += glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3h5_rate_05_95/*.h5')
        # new test-sample corrected reconstructed energy 
        datapaths = glob.glob('/storage/gpfs_data/limadou/data/flight_data/L3_test/L3_repro/*.h5')
elif args.hepp_l:
    data = 'hepp_l'
    det = 'hepp_l_channel_'+args.channel
    # get HEPP data of quiet period 1.-5.08.2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1_L02*_'+str(args.month)+'*.h5') #('/storage/gpfs_data/limadou/data/flight_data/analysis/data/h5/HEPP_august_2018/*.h5')
    if args.day:
        for d in args.day:
            datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1_L02*_'+str(d)+'*.h5')
    # select HEPP data from 22-26.02.2019 (solar quiet period) 
    #datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190222*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190223*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190224*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190225*.h5')
    #datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_1*20190226*.h5')
elif args.hepp_h:
    data = 'hepp_h'
    det = 'hepp_h'
    # get HEPP data of quiet period 1.-5.08.2018
    datapaths = glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_2_L02*'+str(args.month)+'*.h5')
    if args.day:
        for d in args.day:
            datapaths += glob.glob('/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/*HEP_2_L02*'+str(d)+'*.h5')
elif args.noaa:
    data = 'noaa'
    det = 'noaa'
    datapaths = glob.glob('/storage/gpfs_data/limadou/vitalelimadou/run/data/poes_n19_'+str(args.month)+'*_proc.nc.root')
    if args.day:
        for d in args.day:
            datapaths += glob.glob('/storage/gpfs_data/limadou/vitalelimadou/run/data/poes_n19_'+str(d)+'_proc.nc.root')

# sort files by time
datapaths.sort(key=os.path.getmtime)

# define the submission file for condor
args_file = open(home()+'/log/arguments.txt', "w")

# run on single semi-orbits
if not args.merge and not args.ana and not args.select and not args.clean and not args.cluster:
    
    if len(datapaths) < runs:
        print("Only {} files available for reading. ".format(len(datapaths)))
        runs = len(datapaths)
    
    # run 10 recos per job 
    totCmd = [""]*math.ceil(len(datapaths)/10)
    toProcessFiles=0
    for irun,run in enumerate(datapaths):
        #if irun>(runs-1):
        #    break
        
        # test if half-oribit is complete
        # find int in string
        OrbitDateTime = re.findall('\d+', run)
        # calculate time between start and stop of orbit
        duration = 3000
        if args.hepd:
            duration = abs(int(OrbitDateTime[7]) - int(OrbitDateTime[5]))
        elif args.hepp_l or args.hepp_h:
            duration = abs(int(OrbitDateTime[6]) - int(OrbitDateTime[4]))

        if duration < 3000: ## half-orbit not completed
            print('Not full semi-orbit recorded, but only: '+str(duration)+'[mmss]')
            print('Try next run..')
            continue

        buildPath = sharedOutPath()+"data/root/"+args.useVersion+"/"+det+'/'
        if args.originalE:
            buildPath += "originalEnergyBins/"
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
            if args.useVersion == 'v2.1':
                cmd='python3 python/readH5_v2.1.py --inputFile '+str(run)
            elif args.useVersion == 'v3':
                cmd='python3 python/readH5_v3.py --inputFile '+str(run)
            if det=='noaa':
                cmd='python3 python/readSEM2.py --data noaa --inputFile '+str(run)
                if args.useVersion == 'v3':
                    cmd='python3 python/readSEM2_v3.py --data noaa --inputFile '+str(run)
            if args.integral:
                cmd+=' --integral '+str(args.integral)
            if not args.quiet:
                cmd+=' --debug'
            if args.hepd:
                cmd+=' --data hepd'
            elif args.hepp_l:
                cmd+=' --data hepp_l --channel '+args.channel
            elif args.hepp_h:
                cmd+=' --data hepp_h'
            if args.useVersion is not 'v2':
                cmd+=' --useVersion '+args.useVersion
            if args.originalE:
                cmd+=' --useOriginalEnergyBinning '
            print(cmd)
            totCmd[math.floor(toProcessFiles/10)] += cmd+'\n'
            toProcessFiles+=1

    for ijob,job in enumerate(totCmd):
        if job=="":
            continue
        if args.submit:
            exefilename = 'job_%s_%s_%s_%s.sh'%(args.useVersion,det,str(os.path.split(run)[1].replace('.h5','')),ijob)
            if args.integral:
                exefilename = exefilename.replace('_000','_int_'+str(args.integral)+'s')
            exefile = writeExecutionFile(home()+'/log/'+exefilename, job)
            print("Write execution file to:", exefilename)
            args_file.write("%s\n"%(home()+'/log/'+exefilename))
        else:
            os.system(job)
                
    args_file.close()
    if args.submit and toProcessFiles:
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
    mge=['']
    runList=[]
    findOld = []

    if args.hepd:
        # write 
        mge = sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd.root'
        runList = []
        if args.day:
            for d in args.day:
                runList = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/CSES_HEP_DDD_*'+str(d)+'*.root')
                mge = sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(d)+'_'+str(len(runList))+'_runs.root'
                findOld = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/hepd/all_hepd_'+str(d)+'*.root')
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
        pathToFind = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/'
        if args.originalE:
            pathToFind += 'originalEnergyBins/'
        elif args.integral:
            pathToFind = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/'+str(args.integral)+'s/'

        mge = pathToFind+'all_'+det+'.root'
        runList = glob.glob(pathToFind+'CSES_01_HEP_'+str(index)+'_L02*.root') #CSES_01_HEP_1_*.root')

        if args.day:
#            pathToFind = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/'
#            if args.originalE:
#                pathToFind += 'originalEnergyBins/'
#            elif args.integral:
#                pathToFind = sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/'+str(args.integral)+'s/'
            for d in args.day:
                runList = sorted( glob.glob(pathToFind+'CSES_01_HEP_'+str(index)+'_L02_*'+str(d)+'*.root'), key=lambda x:float(x[-46:-41]) )
                mge = pathToFind+'all_'+det+'_'+str(d)+'_'+str(len(runList))+'_runs.root'
                findOld = glob.glob(pathToFind+'all_'+det+'_'+str(d)+'*.root')
                runs = len(runList)

    elif args.noaa:
        runs=0
        pathToFind=sharedOutPath()+'data/root/'+args.useVersion+'/'+det+'/'
        if args.day:
            for d in args.day:
                mge = pathToFind+'poes_n19_'+str(d)+'_proc.nc.root' 
        
    oldruns=0
    if len(findOld)>0:
        for old in findOld:
            head, tail = os.path.split(old)
            oldruns_tmp=list(map(int, re.findall(r'\d+', tail)))[1]
            if oldruns_tmp<runs:
                os.system('rm {}'.format(old))
            if oldruns_tmp>oldruns:
                oldruns=oldruns_tmp

    # merge files only if not already exists and existing file has less inputs  
    if runs>0 and oldruns<runs:
        print("Merge files in: ", mge)
        merge(mge, runList, runs, args.allHists)

    #cmd  = 'hadd -f -k '+mge+' '
    #for ifile in runList:
    #    cmd += ifile+' '
    #cmd += '\n'
    cmd = 'python3 python/writeDayAverages.py --useVersion '+args.useVersion+' --inputFile '+mge+' --data '+det+' '
    if args.originalE:
        cmd += '--originalEnergyBins '
    if args.draw:
        cmd += '--drawHistos '
    if args.integral:
        cmd += '--integral {} '.format(args.integral)
    if args.day:
        for index,d in enumerate(args.day):
            cmd2 = cmd + '--day '+str(d)
        if args.submit:
            SubmitToCondor(cmd2, mge, 1)
        else:
            os.system(cmd2)
    

elif args.merge and args.test:

    tests=['L3_repro'] #'L3h5_orig', 'L3h5_rate', 'L3h5_05_95', 'L3h5_rate_05_95']
    for t in tests:    
        mge = sharedOutPath()+"data/root/"+args.useVersion+"/L3_test/"+t+'/all.root'
        if args.day:
            mge = sharedOutPath()+"data/root/"+args.useVersion+"/L3_test/"+t+'/all_'+str(args.day)+'.root'
            runList = sorted( glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/L3_test/'+t+'/CSES_*'+str(args.day)+'*.root'), key=lambda x:float(x[-10:-4]))
        runs = len(runList)
        if runs>0:
            print("Merge files in: ", mge)
            merge(mge, runList, runs, args.allHists)
 
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
    if det=='noaa':
        fileSnip = 'poes_n19_'
    if args.test:
        detPath = 'L3_test/L3_repro'
        fileSnip = 'all_'
    if not args.day:
        for iday in range(1,32):
            strday = str(iday)
            if iday < 10:
                strday = '0'+str(iday)
            print( sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+strday+'*.root' )
            found = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+strday+'*.root')
            for every in found:
                findFile.append( every )
    elif args.day:
        for d in args.day:
            findFile = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(d)+'*.root')

    for ind,foundFile in enumerate(findFile):
        ind+=1
        strDay = str(ind)
        if ind<10:
            strDay = '0'+str(ind)
        dayint = ind
        cmd = 'python3 python/findHighFluxes.py --inputFile '+foundFile+' --data '+det
        if args.useVersion == 'v3':
            cmd = 'python3 python/findHighFluxes_v3.py --inputFile '+foundFile+' --data '+det
        if args.day:
            if d in args.day:
                cmd += ' --day '+str(d)
                dayint = d
        elif args.month:
            cmd += ' --day '+str(args.month)+strDay
            dayint = int(str(args.month)+strDay)

        if args.submit:
            exefilename = 'job_%s_%s.sh'%(args.useVersion,str(fileSnip+str(dayint)+'.root'))
            if args.integral:
                exefilename = 'job_%s_%s_int_%s.sh'%(args.useVersion,str(fileSnip+str(dayint)+'.root'),args.integral)
            exefile = writeExecutionFile(home()+'/log/'+exefilename, cmd)
            print("Write execution file to:", exefilename)
            args_file.write("%s\n"%(home()+'/log/'+exefilename))
        else:
            os.system(cmd)

    args_file.close()
    if args.submit:
        SubmitListToCondor(args_file)

elif args.cluster:

    findFile = []
    detPath = det
    fileSnip = 'all_highFluxes_'+det+'_'
    thresholdDir = sharedOutPath()+'data/thresholds/'+args.useVersion+'/'+detPath+'/'
    thresholdFile = thresholdDir+'thresholds_'+str(args.month)+'.pkl'
    clusterInput = ''
    cmdTot=''

    if args.month:
        for iday in range(1,32):
            strday = str(iday)
            if iday < 10:
                strday = '0'+str(iday)
            found = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+strday+'*.root')
            for every in found:
                findFile.append( every )
    else:
        findFile = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.day)+'*.root')
        thresholdFile = thresholdDir+'thresholds_'+str(args.day)+'.pkl'
        clusterInput = findFile[0]

    # 1. write thresholds 
    # check if directory exists
    if not os.path.exists( thresholdDir ):
        print("Directory is created: ", thresholdDir)
        os.makedirs( thresholdDir )
    # check if threshold file already exists, if so, skip this step
    if os.path.isfile( thresholdFile ):
        print("Files with thresholds already exists. ")
    # do nothing if 99.99% threshold is used for seeds in clusters 
    elif args.cut=='99perc' or args.cut=='weights' or args.cut=='dummy_cut':
        print("Use in tree stored rms99of99, of weights as thresholds.")
    # determine thresholds..
    else:
        cmd = 'python3 python/threshold_estimation.py --INFILES '
        for ifile in findFile:
            cmd += ifile+' '
        cmd += ' --OUTDIR '+thresholdDir+' --OUTFILE '+os.path.basename(thresholdFile)+' --PLOTON 1\n'
        cmdTot=cmd

    # 2. merge into month
    if args.month:
        mge = sharedOutPath()+"data/root/"+args.useVersion+'/'+detPath+'/'+fileSnip+str(args.month)+'.root'
        clusterInput = mge
        if not os.path.isfile( mge ):
            os.system('hadd -f -k '+mge+' '+str(findFile).replace(']','').replace('[','').replace(',',' '))

    # 3. run clustering
    alreadyDone=False
    clusteredIn = 'clustered_inTime'
    if args.onlyIn:
        if isinstance(args.onlyIn, list):
            for par in args.onlyIn:
                clusteredIn += par
        else:
            clusteredIn += args.onlyIn
    else:
        clusteredIn += 'Only'
    clusterOutdir = sharedOutPath()+'data/root/{0}/{1}/{2}/{3}/{4}s_window/{5}_seeds/'.format(args.useVersion,detPath,clusteredIn,args.cut,args.window,args.seeds)
    if args.doubleSeed:
        clusterOutdir = sharedOutPath()+'data/root/{0}/{1}/{2}/{3}/{4}s_window/{5}_subs_seeds/'.format(args.useVersion,detPath,clusteredIn,args.cut,args.window,args.seeds)

    if not os.path.exists( clusterOutdir ):
        print("Directory is created: ", clusterOutdir)
        os.makedirs( clusterOutdir )
    elif os.path.isfile( clusterOutdir+os.path.basename(clusterInput) ):
        overwrite = input('Clustering already done, do you want to over-write? [y/n] ')
        if overwrite=="n":
            alreadyDone=True

    if not alreadyDone:
        os.system('cp '+clusterInput+' '+clusterOutdir)
        cmd2 = 'python3 python/cluster_finding.py --IN {0}{1} --OUT ./ --DET {2} --CUT {3} --CUTfile {4} --WINDOW {5} --MINNSEED {6}'.format(clusterOutdir,os.path.basename(clusterInput),data,args.cut,thresholdFile,args.window,args.seeds)
        if args.doubleSeed:
            cmd2 += ' --DOUBLESEED '
        if args.onlyIn:
            if isinstance(args.onlyIn, list):
                cmd2 +=' --ONLYIN '
                for par in args.onlyIn:
                    cmd2 += par+' '
            else:
                cmd2 = cmd2+' --ONLYIN '+args.onlyIn

        cmdTot += cmd2
        print(cmdTot)

        # submit to Condor 
        if args.submit:
            exefilename = 'job_cluster_%s_%s_%s_%s.sh'%(os.path.basename(clusterInput),args.cut,str(args.window),str(args.seeds))
            exefile = writeExecutionFile(home()+'/log/'+exefilename, cmdTot)
            print("Write execution file to:", exefilename)
            args_file.write("%s\n"%(home()+'/log/'+exefilename))
            args_file.close()
        
            SubmitListToCondor(args_file)
        # or run locally
        else:
            os.system(cmdTot)


if args.clean:

    detPath = det
    if args.originalE:
        detPath+='/originalEnergyBins'
    elif args.integral:
        detPath+='/'+str(args.integral)+'s'

    datapaths = glob.glob(sharedOutPath()+'data/root/'+args.useVersion+'/'+detPath+'/*000.root')
    if not (args.hepp_l or args.hepp_h):
        print("Clean-up only necessary for HEPP data. ")
        datapaths = []

    # cleanup of double orbit index
    lastOrbit={}
    orbit_file={}

    removeOrbits=0
    for orbit in datapaths:
        h_t = os.path.split(orbit)
        size = os.path.getsize(orbit)
        print( 'Current orbit: ',h_t[1], size )
        numbers=re.findall('\d+', h_t[1])
        orbit_index = int(numbers[4])
        print(orbit_index)

        if orbit_index in lastOrbit:
            print('Orbit index found.')

            if size<lastOrbit[orbit_index]:
                print('found orbit with larger size. Remove this one.')
                os.remove(orbit)
                removeOrbits+=1
            else:
                print('found orbit, but with smaller size (current, previous):')
                print(size, lastOrbit[orbit_index])
                print('remove the other file: ', orbit_file[orbit_index])
                os.remove(orbit_file[orbit_index])
                removeOrbits+=1

        lastOrbit[orbit_index]=size
        orbit_file[orbit_index]=orbit

    print("Number of orbits that get removed: ",removeOrbits )
