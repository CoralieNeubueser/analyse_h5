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
parser.add_argument('--fitFunction', type=str, default='Exp', choices=['Exp','Poisson'], help='Define wich function to use for fitting of distributions.')
parser.add_argument('--day', type=int, default=None, help='Specify a day.')
parser.add_argument('--useVersion', type=str, default='v2.1', help='Specify a data input version.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--integrateEn', action='store_true', help='Merge all fluxes over all energy bins.')
args,_=parser.parse_known_args()

# retrieve binning on L/pitch
numLbin, Lbins = getLbins()
numPbin, Pbins = getPitchBins()
# retrieve threshold
threshold = args.threshold
debug = args.debug
hepd = (args.data == 'hepd')
hepp = (args.data == 'hepp')

print(len(Lbins[0:numLbin]), Lbins[0:numLbin])
print(len(Pbins[0:numPbin-1]), Pbins[0:numPbin-1])

colors = [920, 843, 416, 600, 616, 432, 900, 800, 1,
          920-9, 843-9, 416-9, 600-9, 616-9, 432-9, 900-9, 800-9, 1]

def poisson(x, p):
    """Define poissonian function to fit distribution"""
    newx = round(float(x[0])*float(p[2]),0)
    if newx>0:
        #return math.exp(newx*math.log(p[1]*p[2]) - math.lgamma(newx+1) - p[1]*p[2]) * p[0]
        return p[0] * (math.pow(p[1]*p[2],int(newx))/math.factorial(int(newx))) * math.exp(-(p[1]*p[2]))
    else:
        #if  newx = 0 and mu = 0,  1 is returned
        if (p[1] >= 0):
            return p[0]*math.exp(-p[1])
        else:
            #return a nan for mu < 0 since it does not make sense
            return 0

# this function draws the single histograms per L-alpha cell, and determines the mean/rms etc.
# it is called in parallel 
def getParallelMeans(strHist):
      rootfilename = strHist[0]
      inRoot = r.TFile( rootfilename , 'read' )
      plot = strHist[1]
      cut = strHist[2]
      opt = strHist[3]
      tryFit = strHist[9]
      geo = strHist[10] 
      geomF = strHist[11]
      # draw the histogram
      inRoot.tree.Draw(plot, cut, opt)
      hist = r.gDirectory.Get(strHist[4])
      if not hist:
            print('Drawing didnt work...')
            return
      else:
            entr = hist.GetEntries()
            if args.debug:
                  print("Draw options: ", plot,cut,opt)
                  print("Histogram has {} entries.".format(entr))
            # if the histogram more entries than the threshold  
            if ( entr > threshold ):
                  mean = hist.GetMean()
                  rms = hist.GetRMS()
                  meanErr = hist.GetMeanError()
                  rmsErr = hist.GetRMSError()
                  chi2 = 0 
                  converged = False

                  # fit the distribution
                  if tryFit:
                        clone = hist.Clone("clone")
                        maximumBin = hist.GetMaximumBin()
                        maximum = hist.GetBinCenter(maximumBin)
                        # use a larger RMS, that is not effected by the large number of 0s, in order to increase the fit range
                        clone.GetXaxis().SetRange((maximumBin+1), clone.GetNbinsX())
                        rms_without_zeros = clone.GetRMS()
                        entr_without_zeros = clone.GetEntries()
                        fit_range = (maximum, (maximum+3*rms_without_zeros))
                        # use an exponential decay
                        if args.fitFunction=='Exp':
                            f1 = r.TF1("f1_"+str(strHist[4]),"[0]*exp(-x/[1])",maximum,(maximum+10*rms_without_zeros))
                            f1.SetParameters(1,rms_without_zeros)
                            # use a Poissonian
                        else:
                            f1 = r.TF1("f1_"+str(strHist[4]), poisson, maximum, (maximum+10*rms_without_zeros),3);
                            f1.SetParameters(10, rms_without_zeros*float(geomF), float(geomF))
                            f1.SetParLimits(1, 0, 100)
                            f1.FixParameter(2, float(geomF))
                        f1.SetLineColor(colors[strHist[8]])
                        chi2 = 1000
                        trialStart = 0
                        tau = rms
                        tauErr = rmsErr
                        trialStart = 0
                        trialStartMax = 4
                        if args.fitFunction=='Poisson':
                              trialStartMax = 1

                        while trialStart<trialStartMax:
                              if trialStart>0:
                                    # Set new range on the histgram to extract 2nd/3rd maximum
                                    clone.GetXaxis().SetRange((maximumBin+1), clone.GetNbinsX())
                                    maximumBin = clone.GetMaximumBin()
                              #peakCont = clone.GetBinContent( maximumBin )
                              ## make sure that the tail is not the major determinator of the fit
                              #if peakCont < entr_without_zeros/5.:
                                    ## print("Maximum bin has too low stats.. continue")
                                    #trialStart+=1
                                    #continue
                              # set range 
                              peakPos = hist.GetBinCenter(maximumBin)
                              if args.fitFunction=='Poisson':
                                  peakPos = 0
                              trialEnd = 0
                              while trialEnd<5:
                                    maximumFlux = peakPos + ((2 + trialEnd)*rms_without_zeros) 
                                    fresults = hist.Fit(f1, "QNS", "goff", peakPos, maximumFlux)
                                    if fresults.IsValid() and fresults.Ndf()>0:
                                          tau = f1.GetParameter(1)
                                          tauErr = f1.GetParError(1)
                                          # only converging fits, with tau>0 and tauErr/tau<10% considered
                                          if (tau>0.) and (tau<10*rms_without_zeros) and (tauErr!=0) and (tauErr/tau<0.2) and (fresults.Chi2()/fresults.Ndf()<chi2):
                                                chi2 = fresults.Chi2()/fresults.Ndf()
                                                fit_range = (peakPos, maximumFlux)
                                                converged = True
                                    trialEnd += 1
                              trialStart += 1

                        # set tau as rms
                        if converged:
                              rms = tau
                              rmsErr = tauErr
                              # if fit converged, add function to histogram for draw option
                              if args.drawHistos:
                                    hist.Fit(f1, "Q", "goff", fit_range[0], fit_range[1])
                  # test a RMS99 implementation
                  content=0
                  mbin=0
                  for ibin in range(hist.GetNbinsX()):
                        content+=hist.GetBinContent(ibin)
                        if content/entr<0.99:
                              mbin=ibin
                        else:
                              break
                  rms = hist.GetBinCenter(mbin)+hist.GetBinWidth(mbin)/2.
                  rmsErr = hist.GetBinWidth(mbin)/2.
                  if rms==0:
                        rms = hist.GetBinWidth(mbin)

                  # write mean etc. into txt file
                  with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {:.5f} {:.6f} {:.5f} {:.6f} {:.1f} {} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, chi2, int(converged==True), geo)
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
      lst = getDays(inRoot.tree)
print("To test days: ", lst)

outFilePath = sharedOutPath()+'/data/averages/'+args.useVersion+'/'+str(args.data)+'/'
if args.integral:
    outFilePath = sharedOutPath()+'/data/averages/'+args.useVersion+'/'+str(args.data)+'/'+str(args.integral)+'s/'
if args.fit:
      outFilePath += 'fitted'+args.fitFunction+'/'
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
    thstacks = [r.THStack()] * len(Lbins[0:numLbin])
    tlegends = [r.TLegend()] * len(Lbins[0:numLbin])
    # write averages for all L-p cells / day 
    outFileName = outFilePath+str(d)+'_min_'+str(threshold)+'ev.txt'
    # get average geomagnetic index of this day
    meanGeomIndex = getGeomIndex(data, d)
    print("Write averages in: ", outFileName)
    outFile = open(outFileName, 'w')
    outFile.write('energy L pitch entries mean meanErr rms rmsErr chi2 fit avGeomIndex\n')
    outFile.close()
    count = 0
      
    # get geometrcal factor for meaningful histogram binning 
    geomFactor=1.
    if hepd:
        binWidth = 1./getGeomFactor(0)
        geomFactor = getGeomFactor(0)
    else:
        # hepp data experience much higher flux values in lowest energy bin
        # but this will not be taken into account in the number of bins of the histogram 
        binWidth = 10
        
    for iL,L in enumerate(Lbins[0:numLbin]):
        tlegends[iL] = r.TLegend(0.6,0.5,0.9,.9, 'threshold = '+str(threshold)+' entries')
        thstacks[iL] = r.THStack('stack_'+str(L), 'stack_'+str(L)) 
        
        for iP,P in enumerate(Pbins[0:numPbin-1]):
            writeOut = str('{} {} '.format(L, P))
            histName = 'hist_day_L_'+str(L)+'_p_'+str(P)
            lst_comm = [ filename, str('flux_'+str(L)+'_'+str(P)+'>>'+str(histName)+'(50,0,'+str(50*binWidth)+')'), 'field>25000 && day=='+str(d), 'goff', histName, writeOut, outFileName, th1ds, iP, args.fit, meanGeomIndex, geomFactor ]
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
            lValue = float(re.findall(r"[+-]?\d+\.\d+", name)[1])
            pValue = int(re.findall(r'\d+', name)[5])
            # find index of energy and L
            lIndex = Lbins.index(lValue)
            # fill hist in corresponding stacks
            tlegends[lIndex].AddEntry(h, 'L='+str(lValue)+', #alpha_{eq}='+str(pValue), 'l')
            thstacks[lIndex].Add(h)
                  
            print('legend L entries: ',len(tlegends[0]))

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
                    st.GetXaxis().SetTitle('flux [counts/(s#upoint cm^{2}#upoint sr#upoint MeV)]')
                    st.GetYaxis().SetTitle('# entries')
                    st.SetMinimum(1)
                    tlegends[iest][ifinal].Draw()
                    can.Modified()
                    head, tail = os.path.split( filename.replace('root','pdf').replace('all', 'day_'+str(d)+'_L_'+str(Lbins[ifinal])) )
                    tail = 'day_'+str(d)+'_L_'+str(Lbins[ifinal])+'.pdf'
                    if args.fit:
                        head = head+'/fitted'+args.fitFunction+'/'
                    if args.threshold!=100:
                        head = head+'/'+str(args.threshold)+'ev/'
                    if not os.path.exists(head):
                        os.makedirs(head)

                    outCurves = head+'/'+tail
                    can.Print(outCurves)
