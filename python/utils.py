import os,sys
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

def home():
    return os.getcwd()

def merge(name, listOfFiles):
    ch=r.TChain("tree")  # creates a chain to process a Tree called "T"
    hist2D_l_pitch=r.TH2D("hist2D_l_pitch","hist2D_l_pitch",18,1,10,9,10.,190.)
    hist2D_loc=r.TH2D("hist2D_loc","hist2D_loc",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_flux=r.TH2D("hist2D_loc_flux","hist2D_loc_flux",361,-180.5,180.5,181,-90.5,90.5)
    hist2D_loc_field=r.TH2D("hist2D_loc_field","hist2D_loc_field",361,-180.5,180.5,181,-90.5,90.5)

    for inFile in listOfFiles:
        ch.Add(inFile)
        # open file, get histogram, add without sum over entries
        inRoot = r.TFile( inFile , 'read' )
        h1=r.TH2D()
        h2=r.TH2D()
        h3=r.TH2D()
        h4=r.TH2D()
        r.gDirectory.GetObject('hist2D_l_pitch', h1)
        r.gDirectory.GetObject('hist2D_loc', h2)
        r.gDirectory.GetObject('hist2D_loc_flux', h3)
        r.gDirectory.GetObject('hist2D_loc_field', h4)

        # fill pitch/L histogram with average
        for binx in range(0,h1.GetNbinsX()+1):
            for biny in range(0,h1.GetNbinsY()+1):
                bint = h1.GetBin(binx,biny)
                cont = h1.GetBinContent(bint)
                oldcont = hist2D_l_pitch.GetBinContent(bint)
                if cont!=0.:
                    if oldcont!=0.:
                        hist2D_l_pitch.SetBinContent(bint, (oldcont+cont)/2.)
                    else:
                        hist2D_l_pitch.SetBinContent(bint, cont)
        
        # loop through all lon/lat
        for binx in range(0,h3.GetNbinsX()+1):
            for biny in range(0,h3.GetNbinsY()+1):
                bint = h3.GetBin(binx,biny)
                # fill flux over location histogram with average
                cont = h3.GetBinContent(bint)
                oldcont = hist2D_loc_flux.GetBinContent(bint)
                if cont!=0.:
                    if oldcont!=0.:
                        hist2D_loc_flux.SetBinContent(bint, (oldcont+cont)/2.)
                    else:
                        hist2D_loc_flux.SetBinContent(bint, cont)

                # fill field/time bins, only if not already has a value 
                cont2 = h2.GetBinContent(bint)
                cont4 = h4.GetBinContent(bint)
                if hist2D_loc.GetBinContent(bint)==0.:
                    hist2D_loc.SetBinContent(bint, cont2)
                if hist2D_loc_field.GetBinContent(bint)==0.:
                    hist2D_loc_field.SetBinContent(bint, cont4)


    ch.Merge(name)
    outRoot = r.TFile(name, 'update')
    hist2D_l_pitch.Write()
    hist2D_loc.Write()
    hist2D_loc_flux.Write()
    hist2D_loc_field.Write()

    return outRoot
