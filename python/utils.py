import os,sys
import numpy as np
import ROOT as r
import h5py

def traverse_datasets(hdf_file):

    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            keystr=str(prefix+"/"+key)
            path = str(prefix+"/"+key) #f'{prefix}'
            if isinstance(item, h5py.Dataset): # test for dataset                                                                                                                                                         
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)                                                                                                                                                 
                yield from h5py_dataset_iterator(item, path)

    for path, _ in h5py_dataset_iterator(hdf_file):
        yield path


def draw2D(hist2d, xtitle, ytitle, ztitle, zmin, zmax, out, log):
    r.gROOT.SetBatch(True)
    can=r.TCanvas('can_'+str(hist2d.GetName()))
    if log:
        can.SetLogz()

    hist2d.GetXaxis().SetTitle(xtitle)
    hist2d.GetYaxis().SetTitle(ytitle)
    hist2d.GetZaxis().SetTitle(ztitle)
    hist2d.SetMaximum(zmax)
    hist2d.SetMinimum(zmin)
    hist2d.Draw("colz")
    
    can.Print(out)

def prep2D(hist2d, xtitle, ytitle, ztitle, logz):

    hist2d.GetXaxis().SetTitle(xtitle),
    hist2d.GetYaxis().SetTitle(ytitle)
    hist2d.GetZaxis().SetTitle(ztitle)

    if logz:
        can.SetLogz()

    return hist2d

def home():
    return os.getcwd()

def merge(name, listOfFiles, runs):
    ch=r.TChain("tree")  # creates a chain to process a Tree called "tree"

    # L bins
    l_x_bins = []
    for x in range(0,5):
        l_x_bins.append(1.+0.2*x)
    for x in range(0,9):
        l_x_bins.append(2.+float(x))
    l_bins=len(l_x_bins)-1
    # pitch bins              
    p_x_bins = []
    for p in range(0,10):
        p_x_bins.append(p*20)
    p_bins=len(p_x_bins)

    hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",l_bins,np.array(l_x_bins),9,0,180)
    hist2D_l_pitch_en=r.TH2D("hist2D_l_pitch_en","hist2D_l_pitch_en",l_bins,np.array(l_x_bins),9,0,180)
    hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

    for iF,inFile in enumerate(listOfFiles):
        # limit number of files used for merge
        if iF>runs:
            break
        ch.Add(inFile)
        # open file, get histogram, add without sum over entries
        inRoot = r.TFile( inFile , 'read' )
        h1=r.TH2D()
        h1_en=r.TH2D()
        h2=r.TH2D()
        h3=r.TH2D()
        h4=r.TH2D()
        r.gDirectory.GetObject('hist2D_l_pitch', h1)
        r.gDirectory.GetObject('hist2D_l_pitch_en', h1_en)
        r.gDirectory.GetObject('hist2D_loc', h2)
        r.gDirectory.GetObject('hist2D_loc_flux', h3)
        r.gDirectory.GetObject('hist2D_loc_field', h4)

        # fill pitch/L histogram with sum
        for binx in range(0,h1.GetNbinsX()+1):
            for biny in range(0,h1.GetNbinsY()+1):
                bint = h1.GetBin(binx,biny)
                cont = h1.GetBinContent(bint)
                oldcont = hist2D_l_pitch.GetBinContent(bint)
                if cont!=0.:
                    if oldcont!=0.:
                        hist2D_l_pitch.SetBinContent(bint, oldcont+cont)
                    else:
                        hist2D_l_pitch.SetBinContent(bint, cont)

        # fill pitch/L histogram of entries with sum
        for binx in range(0,h1_en.GetNbinsX()+1):
            for biny in range(0,h1_en.GetNbinsY()+1):
                bint = h1_en.GetBin(binx,biny)
                cont = h1_en.GetBinContent(bint)
                oldcont = hist2D_l_pitch_en.GetBinContent(bint)
                if cont!=0.:
                    if oldcont!=0.:
                        hist2D_l_pitch_en.SetBinContent(bint, oldcont+cont)
                    else:
                        hist2D_l_pitch_en.SetBinContent(bint, cont)

        # loop through all lon/lat
        for binx in range(0,h3.GetNbinsX()+1):
            for biny in range(0,h3.GetNbinsY()+1):
                bint = h3.GetBin(binx,biny)
                # fill flux over location histogram with sum
                cont = h3.GetBinContent(bint)
                oldcont = hist2D_loc_flux.GetBinContent(bint)
                if cont!=0.:
                    if oldcont!=0.:
                        hist2D_loc_flux.SetBinContent(bint, oldcont+cont)
                    else:
                        hist2D_loc_flux.SetBinContent(bint, cont)

                # fill field/time bins, only if not already has a value 
                cont2 = h2.GetBinContent(bint)
                cont4 = h4.GetBinContent(bint)
                if hist2D_loc.GetBinContent(bint)==0.:
                    hist2D_loc.SetBinContent(bint, cont2)
                if hist2D_loc_field.GetBinContent(bint)==0.:
                    hist2D_loc_field.SetBinContent(bint, cont4)
        inRoot.Close()
    
    ch.Merge(name)
    outRoot = r.TFile(name, 'update')
    r.gStyle.SetPadRightMargin(0.5)
    r.gStyle.SetOptStat(0)
    prep2D(hist2D_l_pitch, 'L value', 'pitch [deg]', '#sum{flux}',False)
    prep2D(hist2D_l_pitch_en, 'L value', 'pitch [deg]', '#entries',False)
    prep2D(hist2D_loc, 'longitude', 'latitude', 'time',False)
    prep2D(hist2D_loc_flux, 'longitude', 'latitude', '#LT flux #GT',False)
    prep2D(hist2D_loc_field, 'longitude', 'latitude', 'field [nT]',False)

    hist2D_l_pitch.Write()
    hist2D_l_pitch_en.Write()
    hist2D_loc.Write()
    hist2D_loc_flux.Write()
    hist2D_loc_field.Write()
    outRoot.Close()
    return True

def takeSqrt(hist2D_old,histDiv):

    # hist2D_old.Sumw2()
    # rebin hist of flux to hist of entries
    # if xbins != histDiv.GetNbinsX() or ybins != histDiv.GetNbinsY():
    print("rebin histogram.")
    minX=histDiv.GetXaxis().GetBinLowEdge(histDiv.GetXaxis().GetFirst())
    maxX=histDiv.GetXaxis().GetBinUpEdge(histDiv.GetXaxis().GetLast())
    minY=histDiv.GetYaxis().GetBinLowEdge(histDiv.GetYaxis().GetFirst())
    maxY=histDiv.GetYaxis().GetBinUpEdge(histDiv.GetYaxis().GetLast())
    #print("min bin x : ", minX)
    #print("max bin x : ", maxX)
    #print("min bin y : ", minY)
    #print("max bin y : ", maxY)
    hist2D_rebin = r.TH2D(hist2D_old.GetName()+str("_rebin"), hist2D_old.GetName()+str("_rebin"), histDiv.GetNbinsX(), minX, maxX, histDiv.GetNbinsY(), minY, maxY )
    
    for binx in range(0,hist2D_old.GetNbinsX()+1):
        for biny in range(0,hist2D_old.GetNbinsY()+1):
            bint = hist2D_old.GetBin(binx,biny)
            bincy = hist2D_old.GetYaxis().GetBinCenter(bint)
            bincx = hist2D_old.GetXaxis().GetBinCenter(bint)
            cont = hist2D_old.GetBinContent(bint)
            if cont!=0:
                print(bincx, bincy, cont)
                hist2D_rebin.Fill(bincx, bincy, cont)

    # normalise histogram
    #print("new hist bins x :  ",hist2D_rebin.GetNbinsX())
    #print("norm hist bins x : ",histDiv.GetNbinsX())
    #print("new hist bins y :  ",hist2D_rebin.GetNbinsY())
    #print("norm hist bins y : ",histDiv.GetNbinsY())
    hist2D_rebin.Divide(histDiv)
    xbins_r = hist2D_rebin.GetNbinsX()
    ybins_r = hist2D_rebin.GetNbinsY()
    
    # take sqrt of entries to retrieve RMS
    for binx in range(0,xbins_r+1):
        for biny in range(0,ybins_r+1):
            bint = hist2D_rebin.GetBin(binx,biny)
            cont = hist2D_rebin.GetBinContent(bint)
            hist2D_rebin.SetBinContent(bint,np.sqrt(cont))

    return hist2D_rebin

def divideHists(hist2D,histDiv):

    hist2D.Sumw2()
    
    # rebin hist of flux to hist of entries
    minX=histDiv.GetXaxis().GetBinLowEdge(histDiv.GetXaxis().GetFirst())
    maxX=histDiv.GetXaxis().GetBinUpEdge(histDiv.GetXaxis().GetLast())
    minY=histDiv.GetYaxis().GetBinLowEdge(histDiv.GetYaxis().GetFirst())
    maxY=histDiv.GetYaxis().GetBinUpEdge(histDiv.GetYaxis().GetLast())
    #print("min bin x : ", minX)
    #print("max bin x : ", maxX)
    #print("min bin y : ", minY)
    #print("max bin y : ", maxY)
    hist2D_rebin = r.TH2D(hist2D.GetName()+str("_rebin"), hist2D.GetName()+str("_rebin"), histDiv.GetNbinsX(), minX, maxX, histDiv.GetNbinsY(), minY, maxY )

    for binx in range(0,hist2D.GetNbinsX()+1):
        for biny in range(0,hist2D.GetNbinsY()+1):
            bint = hist2D.GetBin(binx,biny)
            bincx = hist2D.GetXaxis().GetBinCenter(bint)
            bincy = hist2D.GetYaxis().GetBinCenter(bint)
            cont = hist2D.GetBinContent(bint)
            if cont != 0:
                hist2D_rebin.Fill(bincx, bincy, cont)
            
    # normalise histogram
    hist2D_rebin.Divide(histDiv)

    return hist2D_rebin
