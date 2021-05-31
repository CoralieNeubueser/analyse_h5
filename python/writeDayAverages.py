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
parser.add_argument('--data', type=str, choices=['hepd','hepp_l_channel_narrow','hepp_l_channel_wide','hepp_h'], required=True, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
parser.add_argument('--threshold', type=int, default=100, help='Pick a number as minimum statistic in histograms.')
parser.add_argument('--drawHistos', action='store_true', help='Tell if histograms should be drawn.')
parser.add_argument('--fit', action='store_true', help='Use an exponential function to fit the distributions.')
parser.add_argument('--fitFunction', type=str, default='Exp', choices=['Exp','Poisson'], help='Define wich function to use for fitting of distributions.')
parser.add_argument('--day', type=int, default=None, help='Specify a day.')
parser.add_argument('--useVersion', type=str, default='v2', help='Specify a data input version.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
parser.add_argument('--integrateEn', action='store_true', help='Merge all fluxes over all energy bins.')
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
            if debug:
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
                    lastbin=hist.FindLastBinAbove(1)
                    fit_range = (maximum, hist.GetBinCenter(lastbin))
                    # use an exponential decay
                    if args.fitFunction=='Exp':
                        f1 = r.TF1("f1_"+str(strHist[4]),"[0]*exp(-x/[1])",0,hist.GetBinCenter(lastbin))
                        f1.SetParameters(entr,rms_without_zeros)
                        f1.SetParLimits(1, 1e-8,hist.GetBinCenter(lastbin))
                        # use a Poissonian
                    else:
                        f1 = r.TF1("f1_"+str(strHist[4]), poisson, 0, hist.GetBinCenter(lastbin),3);
                        f1.SetParameters(entr, rms_without_zeros*float(geomF), float(geomF))
                        f1.SetParLimits(1, 0, 100)
                        f1.FixParameter(2, float(geomF))
                        
                    f1.SetLineColor(colors[strHist[8]])

                    # first fitting in maximum range
                    fresults = hist.Fit(f1,"RSQ")
                    mean = hist.GetMean()
                    tau = rms
                    tauErr = rmsErr
                    chi2 = 1000
                    failed = True

                    # if fit worked
                    if fresults.Ndf()!=0 and fresults.IsValid():
                        chi2 = fresults.Chi2()/fresults.Ndf()
                        failed = False
                        useFirst = True
                        print("First tau: ", fresults.GetParams()[1])

                    # if the fit didn't work or the chi2 too high, try again..
                    if chi2>10 or failed:
                        fitstart=hist.GetMaximumBin()
                        trial = 0
                        while trial<10:
                            trialbin=trial
                            # make sure only fits across 50% of entries is performed
                            if hist.Integral(trial,lastbin)/hist.GetEntries()<0.5:
                                if debug:
                                    print("Entries in fit range: ", hist.Integral(trial,lastbin))
                                    print("Less than 50%, it will not be continued: ", hist.Integral(trial,lastbin)/hist.GetEntries())
                                break

                            fres = hist.Fit(f1,"SQ","",hist.GetBinCenter(trialbin),hist.GetBinCenter(lastbin))
                            # if fit converged
                            if fres.IsValid() and fres.Ndf()!=0:
                                # check chi2.. 
                                if fres.Chi2()/fres.Ndf() < chi2:
                                    chi2 = fres.Chi2()/fres.Ndf()
                                    fitstart = trialbin
                                    failed = False
                                    useFirst = False
                            trial+=1
                        print('Finished after {} trials.'.format(trial))
                    
                    if not failed and not useFirst:
                        # add last fit with optimised fit range    
                        fresults = hist.Fit(f1,"SQ","",hist.GetBinCenter(fitstart),hist.GetBinCenter(lastbin))

                    if not failed:
                        converged = True
                        tau = fresults.GetParams()[1]
                        print("used tau: ",tau)
                        print("fit results with chi2: ",fresults.Chi2()/fresults.Ndf())
                    else:
                        # if the fit did not converge with good chi2, remove from the histogram
                        hist.RecursiveRemove( hist.FindObject("f1_"+str(strHist[4])) )

                    rms = tau
                    rmsErr = tauErr

                else:  
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
                    line = strHist[5]+'{} {:.8f} {:.8f} {:.8f} {:.8f} {:.1f} {} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, chi2, int(converged==True), geo)
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

en_bins, energies, en_max = 1, [0.], 0.
if args.integrateEn==False:
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

outFilePath = sharedOutPath()+'/data/averages/'+args.useVersion+'/'+str(args.data)+'/'
if args.integral:
    outFilePath = sharedOutPath()+'/data/averages/'+args.useVersion+'/'+str(args.data)+'/'+str(args.integral)+'s/'
if args.integrateEn:
    outFilePath += 'integratedEnergies/' 
if args.fit:
      outFilePath += 'fitted'+args.fitFunction+'/'
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
      # write averages for all L-p cells / day 
      outFileName = outFilePath+str(d)+'_min_'+str(threshold)+'ev.txt'
      # get average geomagnetic index of this day
      meanGeomIndex = getGeomIndex(data, d)
      print("Write averages in: ", outFileName)
      outFile = open(outFileName, 'w')
      outFile.write('energy L pitch entries mean meanErr rms rmsErr chi2 fit avGeomIndex\n')
      outFile.close()
      count = 0
      
      for ien,en in enumerate(energies):
            # get geometrcal factor for meaningful histogram binning, energy dependent flux histo binning for HEPD
            geomFactor=1.
            flux_binWidth = getInverseGeomFactor(args.data,ien)
            geomFactor = getGeomFactor(args.data,ien)

            for iL,L in enumerate(Lbins[0:numLbin-1]):
                  tlegends[ien][iL] = r.TLegend(0.6,0.5,0.9,.9, 'threshold = '+str(threshold)+' entries')
                  thstacks[ien][iL] = r.THStack('stack_'+str(en)+'_'+str(L), 'stack_'+str(en)+'_'+str(L)) 

                  for iP,P in enumerate(Pbins[0:numPbin-1]):
                        writeOut = str('{} {} {} '.format(round(en,1), L, P))
                        histName = 'hist_day_'+str(d)+'_energy_'+str(en)+'_L_'+str(L)+'_p_'+str(P)
                        if args.integrateEn:
                            lst_comm = [ filename, str('flux_'+str(L)+'_'+str(P)+'>>'+str(histName)+'('+str(flux_bins)+',0,'+str(flux_bins*flux_binWidth)+')'), 'field>25000 && day=='+str(d), 'goff', histName, writeOut, outFileName, th1ds, iP, args.fit, meanGeomIndex, geomFactor ]
                        else:
                            lst_comm = [ filename, str('flux_'+str(L)+'_'+str(P)+'>>'+str(histName)+'('+str(flux_bins)+',0,'+str(flux_bins*flux_binWidth)+')'), 'field>25000 && energy_'+str(L)+'_'+str(P)+'=='+str(en)+' && day=='+str(d), 'goff', histName, writeOut, outFileName, th1ds, iP, args.fit, meanGeomIndex, geomFactor ]
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
            print('legend L entries: ',len(tlegends[0]))

            for iest in range(len(energies)):
                  for ifinal in range(len(Lbins[0:numLbin-1])):
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
                              head, tail = os.path.split( filename.replace('root','pdf').replace('all', 'day_'+str(d)+'_energy_'+str(energies[iest])+'_L_'+str(Lbins[ifinal])) )
                              tail = 'day_'+str(d)+'_energy_'+str(energies[iest])+'_L_'+str(Lbins[ifinal])+'.pdf'
                              if args.fit:
                                    head = head+'/fitted'+args.fitFunction+'/'
                              if args.threshold!=100:
                                    head = head+'/'+str(args.threshold)+'ev/'
                              if not os.path.exists(head):
                                    os.makedirs(head)

                              outCurves = head+'/'+tail
                              can.Print(outCurves)
