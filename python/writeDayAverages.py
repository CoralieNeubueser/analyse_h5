
import os,sys,argparse,re
import math
from multiprocessing import Pool, Manager, Array
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
parser.add_argument('--data', type=str, choices=['hepd','hepp_l_channel_narrow','hepp_l_channel_wide','hepp_l_channel_all','hepp_h','noaa_poes19_0degree','noaa_poes19_90degree'], required=True, help='Define patht to data file.')
parser.add_argument('--debug', action='store_true', help='Run in debug mode.')
parser.add_argument('--threshold', type=int, default=100, help='Pick a number as minimum statistic in histograms.')
parser.add_argument('--drawHistos', action='store_true', help='Tell if histograms should be drawn.')
parser.add_argument('--sigma', type=str, default='rms', choices=['rms','rms50','rms99','gauss'], help='Define which method to use for fitting of distributions.')
parser.add_argument('--day', type=int, default=None, help='Specify a day.')
parser.add_argument('--useVersion', type=str, default='v2', choices=['v1','v2','v2.1','v2.2','v3'], help='Specify a data input version.')
parser.add_argument('--integral', type=int, help='Define the time window for integration in seconds.')
#parser.add_argument('--integrateEn', action='store_true', help='Merge all fluxes over all energy bins.')
parser.add_argument('--originalEnergyBins', action='store_true', help='Use original energy binning.')
args,_=parser.parse_known_args()

version=args.useVersion
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
      rfilename = strHist[0]
      inRoot = r.TFile( rfilename , 'read' )
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
      hist_name = strHist[4]
      hist = r.gDirectory.Get(hist_name)
      tree.Draw(plot_divGeomF, cut, opt)
      hist_counts_name = hist_name.replace('hist_','hist_counts_')
      hist_counts = r.gDirectory.Get(hist_counts_name)
      if len(strHist)>14:
          tree.Draw(strHist[14], strHist[16], opt)
          hist2 = r.gDirectory.Get(strHist[15])
          hist.Add( hist2 )
          tree.Draw(strHist[17], strHist[16], opt)
          hist_counts_name2 = strHist[15].replace('hist_','hist_counts_')
          hist_counts2 = r.gDirectory.Get(hist_counts_name2)
          hist_counts.Add( hist_counts2 )
      
      if not hist:
            print('Drawing didnt work...')
            return
      else:
            entr = hist.GetEntries()
            # if the histogram more entries than the threshold  
            if ( entr > threshold ):
                countEntr_noZeros = hist_counts.Integral(2, hist_counts.FindLastBinAbove(0))
                
                mean = hist.GetMean()
                meanErr = hist.GetMeanError()
                rms = hist.GetRMS()
                rmsErr = hist.GetRMSError()

                mean_counts = hist_counts.GetMean()
                meanErr_counts = hist_counts.GetMeanError()
                
                maximumBin = hist.GetBinCenter(hist.FindLastBinAbove(0))
                maximumBin_counts = hist_counts.GetBinCenter(hist.FindLastBinAbove(0))

                if debug:
                    print("Draw options: ", plot,plot_divGeomF,cut,opt)
                    print("Flux histogram has {} entries.".format(entr))
                    print("Counts histogram has {} entries, excluding 0 bin.".format(countEntr_noZeros))
                    #print("Mean/MeanErr/RMS/RMSErr/Mean_counts/MeanErr_counts: {0}/{1}/{2}/{3}/{4}/{5}".format(mean,meanErr,rms,rmsErr,mean_counts,meanErr_counts))
                    if len(strHist)>14:
                        print("2nd draw options: ", strHist[14], strHist[16])
                

                if args.sigma=='rms' or args.sigma=='rms99':
                    # Set range for rms calculation, ecluding 0 counts
                    hist_counts_trunc = hist_counts.Clone('{}_maxX_{}'.format(hist_counts_name,str(maximumBin_counts))) 
                    hist_counts_trunc.GetXaxis().SetRange(2, hist_counts.FindLastBinAbove(0))
                    #hist_counts.GetXaxis().SetRange(0, hist_counts.FindLastBinAbove(0))
                    rms_counts = hist_counts_trunc.GetRMS()
                    rmsErr_counts = hist_counts_trunc.GetRMSError()

                    # test a RMS99 implementation, adding a 99.99% threshold estimate
                    content=0
                    mbin=2 # skip 1 bin
                    mbin_99_of_99=2
                    if countEntr_noZeros>0:
                        while content/countEntr_noZeros<0.99:
                            content = content + hist_counts.GetBinContent(mbin)
                            mbin += 1
                        content=0
                        while content/countEntr_noZeros<0.9999:
                            content = content + hist_counts.GetBinContent(mbin_99_of_99)
                            mbin_99_of_99 += 1

                    rms99 = hist_counts.GetBinCenter(mbin) + hist_counts.GetBinWidth(mbin)/2.
                    rmsErr99 = hist_counts.GetBinWidth(mbin)/2.
                    rms99_of_99 = hist_counts.GetBinCenter(mbin_99_of_99) + hist_counts.GetBinWidth(mbin_99_of_99)/2.
                    rmsErr99_of_99 = hist_counts.GetBinWidth(mbin_99_of_99)/2.

                    if rms99==0:
                        rms99 = rmsErr99
                    if rms99_of_99==0:
                        rms99_of_99 = rmsErr99_of_99

                    # determine weight=rms/rms_99
                    # set range of histogram
                    hist_counts_rms99 = hist_counts.Clone('{}_maxX_{}_noZeros'.format(hist_counts_name,str(maximumBin_counts)))
                    hist_counts_rms99.GetXaxis().SetRange(2, mbin)
                    rms_counts_99 = hist_counts_rms99.GetRMS()
                    weight = 1
                    if rms_counts_99!=0:
                        weight = rms_counts/rms_counts_99

                    # write mean etc. into txt file
                    with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.1f} {:.1f} {:.1f} {:.1f} {:.2f} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, mean_counts, meanErr_counts, rms99, rmsErr99, rms99_of_99, rmsErr99_of_99, weight, geo)
                        txtFile.writelines(line)

                    # prepare the histogram for drawing and add to list
                    hist.SetNameTitle('{0}_maxX_{1} full'.format(hist_name,maximumBin), '{0}_maxX_{1} full'.format(hist_name,maximumBin))
                    
                    hist_counts.SetLineColor(colors[strHist[8]])
                    hist_counts.SetLineWidth(1)
                    hist_counts.SetLineStyle(7)
                    hist_counts_trunc.SetNameTitle('{0}_maxX_{1} trunc'.format(hist_counts_name,maximumBin), '{0}_maxX_{1} trunc'.format(hist_counts_name,maximumBin))
                    hist_counts_trunc.SetLineColor(colors[strHist[8]])
                    hist_counts_trunc.SetLineWidth(2)
                    hist_counts_trunc.SetLineStyle(2)
                    hist_counts_rms99.SetNameTitle('{0}_maxX_{1} rms99'.format(hist_counts_name,maximumBin), '{0}_maxX_{1} rms99'.format(hist_counts_name,maximumBin))
                    hist_counts_rms99.SetLineColor(colors[strHist[8]])
                    hist_counts_rms99.SetLineStyle(1)
                    hist_counts_rms99.SetLineWidth(3)

                    hist.SetLineColor(colors[strHist[8]])
                    hist.SetLineWidth(2)
                    #hist.SetLineStyle(7)

                    if args.sigma=='rms':
                        strHist[7].append(hist)
                    else:
                        strHist[7].append(hist_counts)
                        strHist[7].append(hist_counts_trunc)
                        strHist[7].append(hist_counts_rms99)
                
                elif args.sigma=='rms50':

                    # test a RMS50 implementation on fluxes
                    content=0
                    mbin=0
                    totContent = hist.GetEntries()
                    while content/totContent<0.5:
                        content = content + hist.GetBinContent(mbin)
                        mbin += 1

                    rms50 = hist.GetBinCenter(mbin) + hist.GetBinWidth(mbin)/2.
                    rmsErr50 = hist.GetBinWidth(mbin)/2.
                    rangeFrac = rms50/maximumBin

                    # Set range for rms calculation at x rms50                                                                                                                                                                                
                    hist_trunc5 = hist.Clone('{}_maxX_{}_5rms50_{}'.format(hist_name,maximumBin,str(mbin*5)))
                    hist_trunc5.GetXaxis().SetRange(0, mbin*5)
                    hist_trunc4 = hist.Clone('{}_maxX_{}_4rms50_{}'.format(hist_name,maximumBin,str(mbin*4)))
                    hist_trunc4.GetXaxis().SetRange(0, mbin*4)
                    hist_trunc3 = hist.Clone('{}_maxX_{}_3rms50_{}'.format(hist_name,maximumBin,str(mbin*3)))
                    hist_trunc3.GetXaxis().SetRange(0, mbin*3)
                    hist.SetNameTitle('{0}_maxX_{1} full'.format(hist_name,maximumBin), '{0}_maxX_{1} full'.format(hist_name,maximumBin))
                    
                    mean_5rms50 = hist_trunc5.GetMean()
                    mean_4rms50 = hist_trunc4.GetMean()
                    mean_3rms50 = hist_trunc3.GetMean()

                    rms_5rms50 = hist_trunc5.GetRMS()
                    rms_4rms50 = hist_trunc4.GetRMS()
                    rms_3rms50 = hist_trunc3.GetRMS()

                    # write mean etc. into txt file
                    with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.3f} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, mean_5rms50, rms_5rms50, mean_4rms50, rms_4rms50, mean_3rms50, rms_3rms50, rangeFrac, geo)
                        txtFile.writelines(line)

                    hist.SetLineColor(colors[strHist[8]])
                    hist.SetLineWidth(2)
                    hist_trunc5.SetLineColor(colors[strHist[8]])
                    hist_trunc5.SetFillColor(colors[strHist[8]])
                    hist_trunc5.SetLineWidth(2)
                    hist_trunc5.SetFillStyle(3365)                                                                                                                                                                                                                                
                    hist_trunc4.SetLineColor(colors[strHist[8]])
                    hist_trunc4.SetFillColor(colors[strHist[8]])
                    hist_trunc4.SetLineWidth(2) 
                    hist_trunc4.SetFillStyle(3356)                                                                                                                                                                                                                                
                    hist_trunc3.SetLineColor(colors[strHist[8]])
                    hist_trunc3.SetFillColor(colors[strHist[8]])
                    hist_trunc3.SetLineWidth(2)                   
                    hist_trunc3.SetFillStyle(3144)

                    strHist[7].append(hist)
                    strHist[7].append(hist_trunc5)
                    strHist[7].append(hist_trunc4)
                    strHist[7].append(hist_trunc3)


                elif args.sigma=='gauss':

                    fit = r.TF1('gauss','gaus',hist.GetBinCenter(hist.FindFirstBinAbove(1)), maximumBin)
                    fit.SetParameters(hist.GetEntries(), hist.GetBinCenter(hist.GetMaximumBin()), rms)
                    fit.SetParLimits(1,0,1e6)
                    fit.SetParLimits(2,0,1e3)
                    fit.SetLineColor(colors[strHist[8]])
                    hist.Fit(fit,"RQ")
                    
                    shrinkRange=0
                    chi2Results=[]
                    chi2ndf = -1
                    sigmaMPV = 1
                    if fit.GetNDF()!=0:
                        chi2Results.append(fit.GetChisquare()/fit.GetNDF())
                        chi2ndf = fit.GetChisquare()/fit.GetNDF()
                        if fit.GetParameter(1)!=0:
                            sigmaMPV = fit.GetParameter(2)/fit.GetParameter(1)
                        else:
                            sigmaMPV = fit.GetParameter(2)

                    while (chi2ndf==-1 or chi2ndf>1 or sigmaMPV>1e3) and shrinkRange<5:
                        shrinkRange+=1
                        hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), hist.GetBinCenter(hist.FindLastBinAbove(shrinkRange)))
                        hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), fit.GetParameter(1)+3*fit.GetParameter(2))
                        if fit.GetParError(2)/fit.GetParameter(2)>0.2:
                            hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), fit.GetParameter(1)+2*fit.GetParameter(2))

                        if fit.GetNDF()!=0:
                            chi2Results.append(fit.GetChisquare()/fit.GetNDF())
                            chi2ndf = fit.GetChisquare()/fit.GetNDF()
                            if fit.GetParameter(1)!=0:
                                sigmaMPV =  fit.GetParameter(2)/fit.GetParameter(1)
                            else:
                                sigmaMPV =  fit.GetParameter(2)
                        else:
                            chi2Results.append(999)
                            chi2ndf = -1

                    bestResult = np.amin(np.array(chi2Results))
                    bestResultKey = chi2Results.index(bestResult)
                    print(hist_name)
                    if debug:
                        print('Best key found with:     ', bestResultKey)
                        print('Best results found with: ', bestResult)

                    # Repeat fit for best range
                    hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), hist.GetBinCenter(hist.FindLastBinAbove(bestResultKey)))
                    hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), fit.GetParameter(1)+3*fit.GetParameter(2))
                    if bestResult>1:
                        hist.Fit(fit,"RQ","", hist.GetBinCenter(hist.FindFirstBinAbove(1)), fit.GetParameter(1)+2*fit.GetParameter(2))
                    
                    if fit.GetParameter(1)!=0:
                        sigmaMPV = fit.GetParameter(2)/fit.GetParameter(1)
                    else:
                        sigmaMPV = fit.GetParameter(2)
                    hist.SetLineColor(colors[strHist[8]])
                    #hist.SetLineStyle(2)

                    g2=r.TF1("g2","gaus",hist.GetBinCenter(hist.FindFirstBinAbove(0)), maximumBin)
                    g2.FixParameter(0,fit.GetParameter(0))
                    g2.FixParameter(1,fit.GetParameter(1))
                    g2.FixParameter(2,fit.GetParameter(2))
                    g2.SetLineColor(colors[strHist[8]])
                    g2.SetLineStyle(2)
                    
                    hist.GetListOfFunctions().Add(g2)
                    hist.SetNameTitle('{0}_maxX_{1} full'.format(hist_name,maximumBin), '{0}_maxX_{1} full'.format(hist_name,maximumBin)) 
                    mpv = mean
                    mpvErr = meanErr
                    sig = rms
                    sigErr = rmsErr
                    chi2=-1
                    
                    #                   if sigmaMPV<1e3 and fit.GetNDF()!=0:
                    
                    if fit.GetNDF()!=0:                                                                                                                                                                                                             

                        mpv = g2.GetParameter(1)
                        mpvErr = fit.GetParError(1)
                        sig = g2.GetParameter(2)
                        sigErr = fit.GetParError(2)
                        chi2 = fit.GetChisquare()/fit.GetNDF()

                    # write mean etc. into txt file
                    with open(strHist[6], 'a') as txtFile:
                        line = strHist[5]+'{} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.2f} {:.3f}\n'.format(hist.GetEntries(), mean, meanErr, rms, rmsErr, mean_counts, meanErr_counts, mpv, mpvErr, sig, sigErr, chi2, geo)
                        txtFile.writelines(line)

                    strHist[7].append(hist)

#
#                elif args.sigma=='poisson':
#                    poiss=r.TF1("poiss","TMath::Poisson(x,[1])",hist.GetBinCenter(hist.FindFirstBinAbove(1)), maximumBin)
#                    poiss.SetParameter(1, hist.GetBinCenter(hist.GetMaximumBin()))
#                    #poiss.SetParLimits(1,0,1e6)
#                    #poiss.SetParLimits(2,0,1e3)
#                    poiss.SetLineColor(colors[strHist[8]])
#                    hist.Fit(fit,"RQ")
#
#                    strHist[7].append(hist)

# read tree
filename = args.inputFile
inRoot = r.TFile( filename , 'read' )
tree = inRoot.tree

lst = []
if args.day:
      lst = [int(args.day)]
else:
      lst = getDays(inRoot.tree)
print("To test days:                ", lst)

en_bins, energies, en_max = 1, [0.], 0.
en_bins, energies, en_max = getEnergyBins(det, rebin)

#####
# test if pre-defined energy values are the same as in the root tree
# important for the energy selection!
test_energies = set()
for ev in tree:
    if len(test_energies) < en_bins:
        for e in ev.energy:
            test_energies.add(e)
    else:
        break
test_energies = sorted(test_energies, key=float)
if energies!=test_energies:
    energies=test_energies
print("For energies:                ", energies)
#####

#####
# get maximum flux values per energy/L-alpha cell
allMaxs = {}
for entry in tree:
    for iL,L in enumerate(Lbins[0:numLbin-1]):
        for iP,P in enumerate(Pbins[0:numPbin-1]):
            flux = getattr(tree, 'flux_{0}_{1}'.format(L,P))
            energy_fl =  getattr(tree, 'energy_{0}_{1}'.format(L,P))
            for ifl,fl in enumerate(flux):
                ien = energies.index(energy_fl[ifl])
                if not (ien,iL,iP) in allMaxs:
                    allMaxs[(ien,iL,iP)] = fl
                elif fl > allMaxs[(ien,iL,iP)]:
                    allMaxs[(ien,iL,iP)] = fl
print("{} cells are populated. ".format(len(allMaxs.keys())))
#####

#flux_bins, flux_binWidth = getFluxBins(det)
#count_bins = getCountsBins(det)
#count_bins_corr = count_bins/10
#if args.integral:
#    count_bins_corr = count_bins
#print("Bins for flux distributions: ",count_bins)

writeOutPar = 'mean_counts meanErr_counts rms99 rmsErr99 rms99_of_99 rmsErr99_of_99 weight'

outFilePath = '{0}/data/averages/{1}/{2}/'.format(sharedOutPath(),args.useVersion,args.data)
if args.integral:
    outFilePath = '{0}/data/averages/{1}/{2}/{3}s/'.format(sharedOutPath(), args.useVersion, args.data, args.integral)
if args.sigma=='gauss':
    outFilePath += args.sigma+'/'
    writeOutPar = 'mean_counts meanErr_counts mpv mpvErr sigma sigmaErr chi2'
elif args.sigma=='rms50':
    outFilePath += args.sigma+'/'
    writeOutPar = 'mean_5rms50 rms_5rms50 mean_4rms50 rms_4rms50 mean_3rms50 rms_3rms50 rms50_frac'
if args.originalEnergyBins:
    outFilePath += 'originalEnergyBins/'
if not os.path.exists(outFilePath):
    os.makedirs(outFilePath)
print('Files will be stored in:     ', outFilePath)

# read Data file with geom Indices
data = readGeomIndex()

for d in lst:
      print("Day:                   ", d)
      # list of commands to be run in parallel
      commands = []
      # list of histograms filled in parallel
      hists = Manager().list()
      # stacks and legends to be used for the --drawHistos option
      thstacks = [[r.THStack()] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      tlegends = [[r.TLegend()] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      maxxs = [[0] * len(Lbins[0:numLbin-1]) for x in range(len(energies))]
      # write averages for all L-p cells / day 
      outFileName = '{0}{1}_min_{2}ev.txt'.format(outFilePath,d,threshold)
      # get average geomagnetic index of this day
      meanGeomIndex = getGeomIndex(data, d)
      print("Write averages in:     ", outFileName)
      outFile = open(outFileName, 'w')
      outFile.write('energy L pitch entries mean meanErr rms rmsErr {0} avGeomIndex\n'.format(writeOutPar))
      outFile.close()
      count = 0
      
      for ien,en in enumerate(energies):
            # get geometrcal factor for meaningful histogram binning, energy dependent flux histo binning for HEPD
            geomFactor=1.
            flux_binWidth = getInverseGeomFactor(args.data,ien)
            if args.integral:
                flux_binWidth = flux_binWidth/(args.integral/getTimeBins(args.data))
            geomFactor = getGeomFactor(args.data,ien)

            for iL,L in enumerate(Lbins[0:numLbin-1]):
                  tlegends[ien][iL] = r.TLegend(0.6,0.5,0.9,.9, 'threshold = {} entries'.format(threshold))
                  thstacks[ien][iL] = r.THStack('stack_{0}_{1}'.format(en,L), 'stack_{0}_{1}'.format(en,L)) 
                  
                  for iP,P in enumerate(Pbins[0:numPbin-1]):
                        if not (ien,iL,iP) in allMaxs:
                            continue
                        count_bins = math.ceil(allMaxs[(ien,iL,iP)]/flux_binWidth)
                        #print("maximum flux value: {} number of bins: {}".format(allMaxs[(ien,iL,iP)], count_bins))
                        writeOut = str('{} {} {} '.format(round(en,1), L, P))
                        histName = 'hist_day_{0}_energy_{1}_L_{2}_p_{3}_1'.format(d,en,L,P)
                        histNameCounts = 'hist_counts_day_{0}_energy_{1}_L_{2}_p_{3}_1'.format(d,en,L,P)
                        histCmd = 'flux_{0}_{1}>>{2}({3},0,{4})'.format(L,P,histName,count_bins,count_bins*flux_binWidth)
                        histCmdCounts = 'flux_{0}_{1}*{2}>>{3}({4},0,{5})'.format(L,P,1/flux_binWidth,histNameCounts,count_bins,count_bins)
                        drawOptions = 'field>{4} && day=={0} && energy_{1}_{2}=={3}'.format(d,L,P,en,getSAAcut(det))
                        if args.data=='hepp_l_channel_all':
                            histName = ['hist_day_{0}_energy_{1}_L_{2}_p_{3}_1'.format(d,en,L,P),'hist_day_{0}_energy_{1}_L_{2}_p_{3}_2'.format(d,en,L,P)]
                            histNameCounts = ['hist_counts_day_{0}_energy_{1}_L_{2}_p_{3}_1'.format(d,en,L,P), 'hist_counts_day_{0}_energy_{1}_L_{2}_p_{3}_2'.format(d,en,L,P)]
                            histCmd = ['flux_{0}_{1}>>{2}({3},0,{3})'.format(L,P,histName[0],count_bins), 'flux_{0}_{1}>>{2}({3},0,{3})'.format(L,P,histName[1],count_bins)]
                            histCmdCounts = ['flux_{0}_{1}*{2}>>{3}({4},0,{5})'.format(L,P,getGeomFactor('hepp_l_channel_narrow',ien),histNameCounts[0],count_bins,count_bins), 'flux_{0}_{1}*{2}>>{3}({4},0,{5})'.format(L,P,getGeomFactor('hepp_l_channel_narrow',ien),histNameCounts[1],count_bins_corr,count_bins)]
                            drawOptions = ['field>{4} && day=={0} && energy_{1}_{2}=={3} && channel%2==0'.format(d,L,P,en,getSAAcut(det)), 'field>{4} && day=={0} && energy_{1}_{2}=={3} && channel%2!=0'.format(d,L,P,en,getSAAcut(det))]

                        if count_bins==-1:
                            histCmdCounts = 'flux_{0}_{1}*{2}>>{3}'.format(L,P,geomFactor,histNameCounts)
                            histCmd = 'flux_{0}_{1}>>{2}'.format(L,P,histName)
                        lst_comm = [ filename, histCmd, drawOptions, 'goff', histName, writeOut, outFileName, hists, iP, meanGeomIndex, histCmdCounts, '{0}_{1}_{2}_{3}'.format(d,L,P,en), ien, iL ]
                        if args.data=='hepp_l_channel_all':
                            lst_comm = [ filename, histCmd[0], drawOptions[0], 'goff', histName[0], writeOut, outFileName, hists, iP, meanGeomIndex, histCmdCounts[0], '{0}_{1}_{2}_{3}'.format(d,L,P,en), ien, iL, histCmd[1], histName[1], drawOptions[1], histCmdCounts[1] ]

                        commands.append(lst_comm)
                        count += 1

                        if debug and count>5:
                            break
                  if debug and count>5:
                      break
            if debug and count>5:
                break

      print('Run through {} energy X pitch x L bins.'.format(count))
      pool = Pool(processes=16) # Generally, set to 2*num_cores you have
      # run in parallel
      pool.map( getParallelMeans, commands) 
      pool.close()
      pool.join()
      th1ds=hists
      print("Got {} histograms..s".format(count))
      os.system('chmod -R g+rwx %s'%(outFileName))
      
      if args.drawHistos:

          # write root file with histograms
          rootOutFileName = outFileName.replace('txt','root')
          outrootfile = r.TFile.Open(rootOutFileName, "RECREATE")
          outrootfile.mkdir( "plots" )
          outrootfile.mkdir( "histos" )
          outrootfile.cd( "histos" )
          # loop through list of histograms
          for ih,h in enumerate(th1ds):
              name = h.GetName()
              print(name)
              foundValues = re.findall(r"[+-]?\d+\.\d+", name)
              energyValue = float(foundValues[0])
              lValue = float(foundValues[1])
              pValue = int(re.findall(r'\d+', name)[5])
              maxValue = 500
              if len(foundValues)>2:
                  maxValue = float(foundValues[2])
              label = name[-8:]
              # find index of energy and L
              energyIndex = energies.index(energyValue)
              lIndex = Lbins.index(lValue)
              # fill hist in corresponding stacks
              tlegends[energyIndex][lIndex].AddEntry(h, 'L='+str(lValue)+', #alpha_{eq}='+str(pValue)+', '+label, 'l')
              thstacks[energyIndex][lIndex].Add(h)
              if maxxs[energyIndex][lIndex] < maxValue:
                  maxxs[energyIndex][lIndex] = maxValue
              #print(maxValue)
              #h.Write()
              #getlength = len("hist_day_20210801_energy_0.03999999910593033_L_1.2_p_140")
              #h.SetName(name[:getlength])
              h.Write()
                  
          print('legend E entries: ',len(tlegends))
          print('legend L entries: ',len(tlegends[0]))

          for iest in range(len(energies)):
              for ifinal in range(len(Lbins[0:numLbin-1])):
                  st = thstacks[iest][ifinal]
                  # draw stack, only if contains histograms (entries>threshold)
                  if st.GetNhists()>0:
                      #print(iest,ifinal)
                      #print(energies[iest],Lbins[ifinal])
                      can = r.TCanvas('Energy={0}, L={1}-{2}'.format(energies[iest],Lbins[ifinal],Lbins[ifinal+1]) )
                      can.SetLogy()
                      st.Draw("nostack histe")
                      st.Draw("nostack f same")
                      if args.sigma=='rms99':
                          st.GetXaxis().SetTitle('counts/s')
                      else:
                          st.GetXaxis().SetTitle('flux [counts/s/cm^{2}/sr]')
                          if hepp_l:
                              st.GetXaxis().SetTitle('flux [counts/s/cm^{2}/sr/MeV]')
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
                      outrootfile.cd( "plots" )
                      can.Write()
                  else:
                      continue

          outrootfile.Close()    
          os.system('chmod -R g+rwx %s'%(rootOutFileName))
