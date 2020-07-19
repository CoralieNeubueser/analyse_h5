import os,sys
import numpy as np
import ROOT as r
from utils import *

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
