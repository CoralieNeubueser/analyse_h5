import os,sys,argparse,re
import math
from multiprocessing import Pool, Manager
import multiprocessing as mp
from functools import partial
import matplotlib as plt
import numpy as np
from array import array
# load defined functions  
from utils import *

r.gStyle.SetOptStat(0)
r.gStyle.SetPadRightMargin(0.05);
r.gStyle.SetPadTopMargin(0.05);

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', type=str, help='Define patht to data file.')
parser.add_argument('--data', type=str, choices=['hepd','hepp_l_channel_narrow','hepp_l_channel_wide','hepp_h','noaa'], required=True, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
parser.add_argument('--threshold', type=int, default=100, help='Pick a number as minimum statistic in histograms.')
parser.add_argument('--drawHistos', action='store_true', help='Tell if histograms should be drawn.')
#parser.add_argument('--fit', action='store_true', help='Use an exponential function to fit the distributions.')
#parser.add_argument('--fitFunction', type=str, default='Exp', choices=['Exp','Poisson'], help='Define wich function to use for fitting of distributions.')
parser.add_argument('--day', type=int, default=None, help='Specify a day.')
parser.add_argument('--useVersion', type=str, default='v2', help='Specify a data input version.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
#parser.add_argument('--integrateEn', action='store_true', help='Merge all fluxes over all energy bins.')
parser.add_argument('--originalEnergyBins', action='store_true', help='Use original energy binning.')
args,_=parser.parse_known_args()

# retrieve binning on L/pitch
numLbin, Lbins = getLbins()
numPbin, Pbins = getPitchBins()
# retrieve threshold
threshold = args.threshold
debug = args.debug
det = args.data
hepd = (args.data == 'hepd')
hepp_l = False
if 'hepp_l' in args.data:
    hepp_l = True
    det = 'hepp_l'
hepp_h = (args.data == 'hepp_h')
rebin = False
if (hepp_l or hepp_h) and not args.originalEnergyBins:
    rebin=True

print(len(Lbins[0:numLbin-1]), Lbins[0:numLbin-1])
print(len(Pbins[0:numPbin-1]), Pbins[0:numPbin-1])

colors = [920, 843, 416, 600, 616, 432, 900, 800, 1,
          920-9, 843-9, 416-9, 600-9, 616-9, 432-9, 900-9, 800-9, 1]

# this function draws the single histograms per L-alpha cell, and determines the mean/rms etc.
# it is called in parallel 
def getParallelMeans(strHist):
      rootfilename = strHist[0]
      inRoot = r.TFile( rootfilename , 'read' )
      tree = inRoot.tree
      plot = strHist[1]
      cut = strHist[2]
      opt = strHist[3]
      geo = strHist[9] 
      plot_divGeomF = strHist[10]
      branch = strHist[11]
      energyIndex = strHist[12]
      Lindex = strHist[13]

      # draw the histogram
      tree.Draw(plot, cut, opt)
      hist = r.gDirectory.Get(strHist[4])
      tree.Draw(plot_divGeomF, cut, opt)
      hist_counts_name = strHist[4].replace('hist_','hist_counts_')
      hist_counts = r.gDirectory.Get(hist_counts_name)

      if not hist:
            print('Drawing didnt work...')
            return
      else:
            entr = hist.GetEntries()
            # if the histogram more entries than the threshold  
            if ( entr > threshold ):
                countEntr_noZeros = hist_counts.Integral(2, hist_counts.FindLastBinAbove(0))
                
                if debug:
                    print("Draw options: ", plot,plot_divGeomF,cut,opt)
                    print("Flux histogram has {} entries.".format(entr))
                    print("Counts histogram has {} entries, excluding 0 bin.".format(countEntr_noZeros))

                mean = hist.GetMean()
                meanErr = hist.GetMeanError()
                rms = hist.GetRMS()
                rmsErr = hist.GetRMSError()

                mean_counts = hist_counts.GetMean()
                meanErr_counts = hist_counts.GetMeanError()
                
                # Set range for rms calculation, ecluding 0 counts
                maximumBin = hist_counts.GetBinCenter(hist_counts.FindLastBinAbove(0))
                hist_counts_trunc = hist_counts.Clone('{}_maxX_{}'.format(hist_counts_name,str(maximumBin))) 
                hist_counts_trunc.GetXaxis().SetRange(2, hist_counts.FindLastBinAbove(0))
                #hist_counts.GetXaxis().SetRange(0, hist_counts.FindLastBinAbove(0))
                rms_counts = hist_counts_trunc.GetRMS()
                rmsErr_counts = hist_counts_trunc.GetRMSError()

                # test a RMS99 implementation, adding a 99.99% threshold estimate
                content=0
                mbin=0
                mbin_99_of_99=0
                if countEntr_noZeros>0:
                    # skip 1 bin
                    for ibin in range(2, hist_counts.GetNbinsX()):
                        content+=hist_counts.GetBinContent(ibin)
                        if content/countEntr_noZeros<0.99:
                            mbin=ibin
                        elif content/countEntr_noZeros<0.9999:
                            mbin_99_of_99=ibin
                        else:
                            break

                rms99 = hist_counts.GetBinCenter(mbin) + hist_counts.GetBinWidth(mbin)/2.
                rmsErr99 = hist_counts.GetBinWidth(mbin)/2.
                rms99_of_99 = hist_counts.GetBinCenter(mbin_99_of_99) + hist_counts.GetBinWidth(mbin_99_of_99)/2.
                rmsErr99_of_99 = hist_counts.GetBinWidth(mbin_99_of_99)/2.

                if rms99==0:
                    rms99 = rmsErr99
                if rms99_of_99==0:
                    rms99_of_99 = rmsErr99_of_99

                # detemine weight=rms/rms_99
                # set range of histogram
                hist_counts_rms99 = hist_counts.Clone('{}_maxX_{}_noZeros'.format(hist_counts_name,str(maximumBin)))
                hist_counts_rms99.GetXaxis().SetRange(2, mbin)
                rms_counts_99 = hist_counts_rms99.GetRMS()
                weight = 1
                if rms_counts_99!=0:
                    weight = rms_counts/rms_counts_99

                # write mean etc. into txt file
                with open(strHist[6], 'a') as txtFile:
                    line = strHist[5]+'{} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.1f} {:.1f} {:.1f} {:.1f} {:.2f} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, mean_counts, meanErr_counts, rms99, rmsErr99, rms99_of_99, rmsErr99_of_99, weight, geo)
                    txtFile.writelines(line)
            
                plotName = hist_counts_name 

                # prepare the histogram for drawing and add to list
                hist.SetNameTitle('{0}_maxX_{1} full'.format(plotName,maximumBin), '{0}_maxX_{1} full'.format(plotName,maximumBin)) 
                hist_counts.SetLineColor(colors[strHist[8]])
                hist_counts.SetLineWidth(1)
                hist_counts.SetLineStyle(7)
                hist_counts_trunc.SetNameTitle('{0}_maxX_{1} trunc'.format(plotName,maximumBin), '{0}_maxX_{1} trunc'.format(plotName,maximumBin))
                hist_counts_trunc.SetLineColor(colors[strHist[8]])
                hist_counts_trunc.SetLineWidth(2)
                hist_counts_trunc.SetLineStyle(2) 
                hist_counts_rms99.SetNameTitle('{0}_maxX_{1} rms99'.format(plotName,maximumBin), '{0}_maxX_{1} rms99'.format(plotName,maximumBin))
                hist_counts_rms99.SetLineColor(colors[strHist[8]])
                hist_counts_rms99.SetLineStyle(1)
                hist_counts_rms99.SetLineWidth(3)
                
                if args.drawHistos:
                    strHist[7].append(hist_counts)
                    strHist[7].append(hist_counts_trunc)
                    strHist[7].append(hist_counts_rms99)
                
filename = args.inputFile
# read tree                
inRoot = r.TFile( filename , 'read' )
lst = []
if args.day:
      lst = [int(args.day)]
else:
      lst = getDays(inRoot.tree)
print("To test days: ", lst)

en_bins, energies, en_max = 1, [0.], 0.
en_bins, energies, en_max = getEnergyBins(det, rebin)
print(energies)
# test if pre-defined energy values are the same as in the root tree
# important for the energy selection!
test_energies = set()
for ev in inRoot.tree:
    if len(test_energies) < en_bins:
        for e in ev.energy:
            test_energies.add(e)
    else:
        break
test_energies = sorted(test_energies, key=float)
if energies!=test_energies:
    energies=test_energies
print("For energies: ", energies)

flux_bins, flux_binWidth = getFluxBins(det)
count_bins = getCountsBins(det)

outFilePath = '{0}/data/averages/{1}/{2}/'.format(sharedOutPath(),args.useVersion,args.data)
if args.integral:
    outFilePath = '{0}/data/averages/{1}/{2}/{3}s/'.format(sharedOutPath(), args.useVersion, args.data, args.integral)
#if args.integrateEn:
#    outFilePath += 'integratedEnergies/' 
#if args.fit:
#      outFilePath += 'fitted'+args.fitFunction+'/'
if args.originalEnergyBins:
    outFilePath += 'originalEnergyBins/'
if not os.path.exists(outFilePath):
    os.makedirs(outFilePath)
print('Files will be stored in: ', outFilePath)

# read Data file with geom Indices
data = readGeomIndex()

for d in lst:
      print("Day: ", d)
      # list of commands to be run in parallel
      commands = []
      # list of histograms filled in parallel
      th1ds = Manager().list()
      # stacks and legends to be used for the --drawHistos option
      thstacks = [[r.THStack()] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      tlegends = [[r.TLegend()] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      maxxs = [[0] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      # write averages for all L-p cells / day 
      outFileName = '{0}{1}_min_{2}ev.txt'.format(outFilePath,d,threshold)
      # get average geomagnetic index of this day
      meanGeomIndex = getGeomIndex(data, d)
      print("Write averages in: ", outFileName)
      outFile = open(outFileName, 'w')
      outFile.write('energy L pitch entries mean meanErr rms rmsErr mean_counts meanErr_counts rms99 rmsErr99 rms99_of_99 rmsErr99_of_99 weight avGeomIndex\n')
      outFile.close()
      count = 0
      
      for ien,en in enumerate(energies):
            # get geometrcal factor for meaningful histogram binning, energy dependent flux histo binning for HEPD
            geomFactor=1.
            flux_binWidth = getInverseGeomFactor(args.data,ien)
            geomFactor = getGeomFactor(args.data,ien)

            for iL,L in enumerate(Lbins[0:numLbin-1]):
                  tlegends[ien][iL] = r.TLegend(0.6,0.5,0.9,.9, 'threshold = {} entries'.format(threshold))
                  thstacks[ien][iL] = r.THStack('stack_{0}_{1}'.format(en,L), 'stack_{0}_{1}'.format(en,L)) 
                  
                  for iP,P in enumerate(Pbins[0:numPbin-1]):
                        writeOut = str('{} {} {} '.format(round(en,1), L, P))
                        histName = 'hist_day_{0}_energy_{1}_L_{2}_p_{3}'.format(d,en,L,P) 
                        histNameCounts = 'hist_counts_day_{0}_energy_{1}_L_{2}_p_{3}'.format(d,en,L,P)
                        histCmdCounts = 'flux_{0}_{1}*{2}>>{3}({4},0,{4})'.format(L,P,geomFactor,histNameCounts,count_bins)
                        if count_bins==-1:
                            histCmdCounts ='flux_{0}_{1}*{2}>>{3}'.format(L,P,geomFactor,histNameCounts)
                        lst_comm = [ filename, 'flux_{0}_{1}>>{2}'.format(L,P,histName), 'field>{4} && day=={0} && energy_{1}_{2}=={3}'.format(d,L,P,en,getSAAcut(det)), 'goff', histName, writeOut, outFileName, th1ds, iP, meanGeomIndex, histCmdCounts, '{0}_{1}_{2}_{3}'.format(d,L,P,en), ien, iL ]

                        commands.append(lst_comm)
                        count += 1

                        if debug and count>5:
                            break
            
                  if debug and count>5:
                      break
            if debug and count>5:
                break

      print('Run through {} energy X pitch x L bins.'.format(count))
      pool = Pool(processes=8) # Generally, set to 2*num_cores you have
      # run in parallel
      pool.map( getParallelMeans, commands) 
      pool.close()
      pool.join()
      print("Got {} histograms.. efficiency with thr={} of {:.2f}".format(len(th1ds)/3, threshold, len(th1ds)/3/count*100))
      
      if args.drawHistos:
            # write root file with histograms
            rootOutFileName = outFileName.replace('txt','root')
            outrootfile = r.TFile.Open(rootOutFileName, "RECREATE")
            # loop through list of histograms
            for ih,h in enumerate(th1ds):
                  name = h.GetName()
                  print(name)
                  energyValue = float(re.findall(r"[+-]?\d+\.\d+", name)[0])
                  lValue = float(re.findall(r"[+-]?\d+\.\d+", name)[1])
                  pValue = int(re.findall(r'\d+', name)[5])
                  maxValue = 500
                  if len((re.findall(r'[+-]?\d+\.\d+', name)))>2:
                      maxValue = float(re.findall(r'[+-]?\d+\.\d+', name)[2])
                  label = name[-5:]
                  # find index of energy and L
                  energyIndex = energies.index(energyValue)
                  lIndex = Lbins.index(lValue)
                  # fill hist in corresponding stacks
                  tlegends[energyIndex][lIndex].AddEntry(h, 'L='+str(lValue)+', #alpha_{eq}='+str(pValue)+', '+label, 'l')
                  thstacks[energyIndex][lIndex].Add(h)
                  if maxxs[energyIndex][lIndex] < maxValue:
                      maxxs[energyIndex][lIndex] = maxValue
                  h.Write()

            print('legend E entries: ',len(tlegends))
            print('legend L entries: ',len(tlegends[0]))

            for iest in range(len(energies)):
                  for ifinal in range(len(Lbins[0:numLbin-1])):
                        st = thstacks[iest][ifinal]
                        # draw stack, only if contains histograms (entries>threshold)
                        if not st.GetNhists()>0:
                              continue
                        else:
                              can = r.TCanvas('Energy={0}, L={1}-{2}'.format(energies[iest],Lbins[ifinal],Lbins[ifinal+1]) )
                              can.SetLogy()
                              st.Draw("nostack, histe")
                              st.Draw("nostack, f same")
                              st.GetXaxis().SetTitle('counts [1/s]')
                              st.GetYaxis().SetTitle('# entries')
                              st.SetMinimum(1)
                              st.GetXaxis().SetLimits(0, maxxs[iest][ifinal])
                              tlegends[iest][ifinal].Draw()
                              can.Modified()
                              head, tail = os.path.split( filename.replace('root','pdf').replace('all', 'day_{0}_energy_{1}_L_{2}'.format(d, energies[iest], Lbins[ifinal])) )
                              tail = 'day_{0}_energy_{1}_L_{2}.pdf'.format(d, energies[iest], Lbins[ifinal])
                              if threshold!=100:
                                    head = '{0}/{1}ev/'.format(head,threshold)
                              if not os.path.exists(head):
                                    os.makedirs(head)

                              outCurves = '{}/{}'.format(head,tail)
                              #can.Print(outCurves)
                              can.Write()
            outrootfile.Close()    
