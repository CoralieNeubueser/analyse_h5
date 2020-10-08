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
parser.add_argument('--data', type=str, choices=['hepd','hepp'], required=True, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
parser.add_argument('--threshold', type=int, default=100, help='Pick a number as minimum statistic in histograms.')
parser.add_argument('--drawHistos', action='store_true', help='Tell if histograms should be drawn.')
parser.add_argument('--fit', action='store_true', help='Use an exponential function to fit the distributions.')
parser.add_argument('--day', type=int, help='Specify a day.')
args,_=parser.parse_known_args()

# retrieve binning on L/pitch
numLbin, Lbins = getLbins()
numPbin, Pbins = getPitchBins()
# retrieve threshold
threshold = args.threshold
debug = args.debug
hepd = args.data == 'hepd'
hepp = args.data == 'hepp'

print(len(Lbins[0:numLbin]), Lbins[0:numLbin])
print(len(Pbins[0:numPbin-1]), Pbins[0:numPbin-1])

colors = [920, 843, 416, 600, 616, 432, 900, 800, 1,
          920-9, 843-9, 416-9, 600-9, 616-9, 432-9, 900-9, 800-9, 1]

# this function draws the single histograms per L-alpha cell, and determines the mean/rms etc.
# it is called in parallel 
def getParallelMeans(strHist):
      rootfilename = strHist[0]
      inRoot = r.TFile( rootfilename , 'read' )
      plot = strHist[1]
      cut = strHist[2]
      opt = strHist[3]
      exp = strHist[9]
      if args.debug:
            print("Draw options: ", plot,cut,opt)
      # draw the histogram
      inRoot.tree.Draw(plot, cut, opt)
      hist = r.gDirectory.Get(strHist[4])
      if not hist:
            print('Drawing didnt work...')
            return
      else:
            entr = hist.GetEntries()
            # if the histogram more entries than the threshold  
            if ( entr > threshold ):
                  mean = hist.GetMean()
                  rms = hist.GetRMS()
                  meanErr = hist.GetMeanError()
                  rmsErr = hist.GetRMSError()
                  chi2 = 0 
                  converged = False

                  # fit with exponential
                  if exp:
                        maximum = hist.GetBinCenter(hist.GetMaximumBin())
                        fit_range = (maximum, (maximum+3*rms)) 
                        f1 = r.TF1("f1_"+str(strHist[4]),"[0]*exp(-x/[1])",maximum,(maximum+10*rms))
                        f1.SetParameters(1,rms)
                        f1.SetLineColor(colors[strHist[8]])
                        chi2 = 1000
                        trialStart = 0
                        tau = rms
                        tauErr = rmsErr
                        trialStart = 0
                        while trialStart<3:
                              # make sure that the tail is not the major determinator of the fit
                              if hist.GetBinContent(hist.GetMaximumBin()+trialStart) < entr/5.:
                                    #print("Maximum bin has too low stats.. continue")
                                    trialStart+=1
                                    continue
                              peakPos = hist.GetBinCenter(hist.GetMaximumBin()+trialStart)
                              trialEnd = 0
                              while trialEnd<5:
                                    maximumFlux = peakPos + ((2 + trialEnd)*rms) 
                                    maxEntr = hist.GetBinCenter(hist.FindLastBinAbove(1))
                                    fresults = hist.Fit(f1, "QNS", "goff", peakPos, maximumFlux)
                                    if fresults.IsValid() and fresults.Ndf()>0:
                                          tau = f1.GetParameter(1)
                                          tauErr = f1.GetParError(1)
                                          # only converging fits, with tau>0 and tauErr/tau<10% considered
                                          if (tau>0.) and (tau<maxEntr) and (tauErr!=0) and (tauErr/tau<0.2) and (fresults.Chi2()/fresults.Ndf()<chi2):
                                                chi2 = fresults.Chi2()/fresults.Ndf()
                                                fit_range = (peakPos, maximumFlux)
                                                converged = True
                                    trialEnd += 1
                              trialStart += 1

                        # set tau as rms
                        rms = tau
                        rmsErr = tauErr
                        # if fit converged, add function to histogram for draw option
                        if converged and args.drawHistos:
                              hist.Fit(f1, "Q", "goff", fit_range[0], fit_range[1])

                  # write mean etc. into txt file
                  with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {} {} {} {} {} {}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, chi2, int(converged==True))
                        txtFile.writelines(line)
                  # prepare the histogram for drawing and add to list
                  hist.SetName(strHist[4]) 
                  hist.SetLineColor(colors[strHist[8]])
                  hist.SetLineWidth(2)
                  strHist[7].append(hist)

filename = args.inputFile
# read tree                
inRoot = r.TFile( filename , 'read' )
lst = []
if args.day:
      lst = [int(args.day)]
else:
      getDays(inRoot.tree)
en_bins, energies, en_max = getEnergyBins(hepd, hepp)
print("To test days: ", lst)

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

outFilePath = home()+'/data/averages/'+str(args.data)+'/'
if args.fit:
      outFilePath += 'fittedExp/'
if not os.path.exists(outFilePath):
    os.makedirs(outFilePath)
print(outFilePath)

for d in lst:
      print("Day: ", d)
      # list of commands to be run in parallel
      commands = []
      # list of histograms filled in parallel
      th1ds = Manager().list() 
      # stacks and legends to be used for the --drawHistos option
      thstacks = [[r.THStack()] * len(Lbins[0:numLbin]) for x in range(len(energies))]
      tlegends = [[r.TLegend()] * len(Lbins[0:numLbin]) for x in range(len(energies))]
      # write averagesfor all L-p cells / day 
      outFileName = outFilePath+str(d)+'_min_'+str(threshold)+'ev.txt'
      print("Write averages in: ", outFileName)
      outFile = open(outFileName, 'w')
      outFile.write('energy L pitch entries mean meanErr rms rmsErr chi2 fit\n')
      outFile.close()
      count = 0
      
      for ien,en in enumerate(energies):
            # get geometrcal factor for meaningful histogram binning 
            if hepd:
                  binWidth = 1./getGeomFactor(ien)
            else:
                  # hepp data experience much higher flux values in lowest energy bin
                  # but this will not be taken into account in the number of bins of the histogram
                  binWidth = 10
            for iL,L in enumerate(Lbins[0:numLbin]):
                  tlegends[ien][iL] = r.TLegend(0.6,0.5,0.9,.9, 'threshold = '+str(threshold)+' entries')
                  thstacks[ien][iL] = r.THStack('stack_'+str(en)+'_'+str(L), 'stack_'+str(en)+'_'+str(L)) 

                  for iP,P in enumerate(Pbins[0:numPbin-1]):
                        writeOut = str('{} {} {} '.format(round(energies[ien],1), L, P))
                        histName = 'hist_day_'+str(d)+'_energy_'+str(en)+'_L_'+str(L)+'_p_'+str(P)
                        lst_comm = [ filename, str('flux_'+str(L)+'_'+str(P)+'>>'+str(histName)+'(50,0,'+str(50*binWidth)+')'), 'field>25000 && energy=='+str(en)+' && day=='+str(d), 'goff', histName, writeOut, outFileName, th1ds, iP, args.fit ]
                        commands.append(lst_comm)
                        count += 1

                        if debug and count > 0:
                            break

      print('Run through {} energy X pitch x L bins.'.format(count))
      pool = Pool(processes=10) # Generally, set to 2*num_cores you have
      # run in parallel
      pool.map( getParallelMeans, commands) 
      pool.close()
      pool.join()
      print("Got {} histograms.. efficiency with thr={} of {:.2f}".format(len(th1ds), args.threshold, len(th1ds)/count*100))

      if args.drawHistos:
            if args.fit:
                  r.gStyle.SetOptStat(1111)
            # loop through list of histograms
            for h in th1ds:
                  name = h.GetName()
                  energyValue = float(re.findall(r"[+-]?\d+\.\d+", name)[0])
                  lValue = float(re.findall(r"[+-]?\d+\.\d+", name)[1])
                  pValue = int(re.findall(r'\d+', name)[5])
                  # find index of energy and L
                  energyIndex = energies.index(energyValue)
                  lIndex = Lbins.index(lValue)
                  # fill hist in corresponding stacks
                  tlegends[energyIndex][lIndex].AddEntry(h, 'L='+str(lValue)+', #alpha_{eq}='+str(pValue), 'l')
                  thstacks[energyIndex][lIndex].Add(h)
                  
            print('legend E entries: ',len(tlegends))
            print('legend L entries: ',len(tlegends[1]))

            for iest in range(len(energies)):
                  for ifinal in range(len(Lbins[0:numLbin])):
                        st = thstacks[iest][ifinal]
                        # draw stack, only if contains histograms (entries>threshold)
                        if not st.GetNhists()>0:
                              continue
                        else:
                              can = r.TCanvas()
                              can.SetLogy()
                              st.Draw("nostack, histe")
                              st.Draw("nostack, f same")
                              st.GetXaxis().SetTitle('flux [counts/(s#upoint cm^{2}#upoint sr)]')
                              st.GetYaxis().SetTitle('# entries')
                              st.SetMinimum(1)
                              tlegends[iest][ifinal].Draw()
                              can.Modified()
                              outCurves = filename.replace('root','pdf').replace('all', 'day_'+str(d)+'_energy_'+str(energies[iest])+'_L_'+str(Lbins[ifinal]))
                              can.Print(outCurves)
