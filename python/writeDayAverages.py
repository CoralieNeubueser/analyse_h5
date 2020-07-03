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
args,_=parser.parse_known_args()

# retrieve binning on L/pitch
numLbin, Lbins = getLbins()
numPbin, Pbins = getPitchBins()
# retrieve threshold
threshold = args.threshold
debug = args.debug

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
      # draw the histogram
      inRoot.tree.Draw(plot, cut, opt)
      hist = r.gDirectory.Get(strHist[4])
      if not hist:
            print('Drawing didnt work...')
            return
      else:
            # if the histogram more entries than the threshold  
            if ( hist.GetEntries() > threshold ):
                  # write mean etc. into txt file
                  with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {} {} {} {}\n'.format(hist.GetEntries(), hist.GetMean(), hist.GetMeanError(), hist.GetRMS(), hist.GetRMSError())
                        txtFile.writelines(line)
                  # prepare the histogram for drawing and add to list
                  hist.SetName(strHist[4]) 
                  hist.SetLineColor(colors[strHist[8]])
                  hist.SetLineWidth(2)
                  strHist[7].append(hist)

filename = args.inputFile
# read tree                
inRoot = r.TFile( filename , 'read' )
lst = getDays(inRoot.tree)
en_bins, energies, en_max = getEnergyBins(True, False)
print("To test days: ", lst)
print("For energies: ", energies)

outFilePath = home()+'/averages/'+str(args.data)+'/'
if not os.path.exists(outFilePath):
    os.makedirs(outFilePath)

for d in lst:
      print("Day: ", d)
      # list of commands to be run in parallel
      commands = []
      # list of histograms filled in parallel
      th1ds = Manager().list() 
      # stacks and legends to be used for the --drawHistos option
      thstacks = [[r.THStack()] * len(Lbins[0:numLbin]) for x in range(len(energies))]
      tlegends = [[r.TLegend()] * len(Lbins[0:numLbin]) for x in range(len(energies))]
      tlines = [[[]] * len(Lbins[0:numLbin]) for x in range(len(energies))]
      # write averagesfor all L-p cells / day 
      outFileName = outFilePath+str(d)+'.txt'
      outFile = open(outFileName, 'w')
      outFile.write('energy L pitch entries mean meanErr rms rmsErr \n')
      outFile.close()
      count = 0
      
      for ien,en in enumerate(energies):
            # get geometrcal factor for meaningful histogram binning 
            binWidth = 1./getGeomFactor(ien)

            for iL,L in enumerate(Lbins[0:numLbin]):
                  tlegends[ien][iL] = r.TLegend(0.6,0.5,0.95,.9, 'threshold = '+str(threshold)+' entries')
                  thstacks[ien][iL] = r.THStack('stack_'+str(en)+'_'+str(L), 'stack_'+str(en)+'_'+str(L)) 

                  for iP,P in enumerate(Pbins[0:numPbin-1]):
                        writeOut = str('{} {} {} '.format(en, L, P))
                        histName = 'hist_day_'+str(d)+'_energy_'+str(en)+'_L_'+str(L)+'_p_'+str(P)
                        lst_comm = [ filename, str('flux_'+str(L)+'_'+str(P)+'>>'+str(histName)+'(50,0,'+str(50*binWidth)+')'), 'field>25000 && energy=='+str(en)+' && day=='+str(d), 'goff', histName, writeOut, outFileName, th1ds, iP ]
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
            # loop through list of histograms
            for h in th1ds:
                  name = h.GetName()
                  energyValue = float(re.findall(r"[+-]?\d+\.\d+", name)[0])
                  lValue = float(re.findall(r"[+-]?\d+\.\d+", name)[1])
                  pValue = int(re.findall(r'\d+', name)[5])
                  # find index of energy and L
                  energyIndex = energies.index(energyValue)
                  lIndex = Lbins.index(lValue)
                  # add line for 5 sigma cut
                  fiveSig = h.GetMean() + 5*h.GetRMS()
                  line = r.TLine(fiveSig,h.GetMinimum(),fiveSig,h.GetMaximum());
                  line.SetLineColor(colors[Pbins.index(pValue)])
                  # fill hist in corresponding stacks
                  tlegends[energyIndex][lIndex].AddEntry(h, 'L='+str(lValue)+', p='+str(pValue), 'l')
                  thstacks[energyIndex][lIndex].Add(h)
                  tlines[energyIndex][lIndex].append(line)

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
                              st.GetXaxis().SetTitle('flux [counts/(s#upoint cm^{2}#upoint sr#upoint MeV)]')
                              st.GetYaxis().SetTitle('# entries')
                              st.SetMinimum(1)
                              tlegends[iest][ifinal].Draw()
                              for li in tlines[iest][ifinal]:
                                    li.Draw()
                              can.Modified()
                              outCurves = filename.replace('root','pdf').replace('all', 'day_'+str(d)+'_energy_'+str(energies[iest])+'_L_'+str(Lbins[ifinal]))
                              can.Print(outCurves)
