import ROOT

"""
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(True)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetPadTopMargin(0.06);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.12);
    ROOT.gStyle.SetPadRightMargin(.15)
    ROOT.gStyle.SetLabelColor(1, "XYZ");
    ROOT.gStyle.SetLabelFont(42, "XYZ");
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
    ROOT.gStyle.SetLabelSize(0.05, "XYZ");
    ROOT.gStyle.SetTitleSize(0.05, "XYZ");
    ROOT.gStyle.SetTitleOffset(1.0, "X");
    ROOT.gStyle.SetTitleOffset(1.1, "Y");
    ROOT.gStyle.SetTitleOffset(1.0, "Z");
    ROOT.gStyle.SetAxisColor(1, "XYZ");
    ROOT.gStyle.SetStripDecimals(True);
    ROOT.gStyle.SetTickLength(0.03, "XYZ");
    ROOT.gStyle.SetNdivisions(510, "XYZ");
    ROOT.gStyle.SetPadTickX(0);
    ROOT.gStyle.SetPadTickY(0);
    ROOT.gStyle.SetMarkerStyle(20);
    ROOT.gStyle.SetHistLineColor(1);
    ROOT.gStyle.SetHistLineStyle(1);
    ROOT.gStyle.SetHistLineWidth(3);
    ROOT.gStyle.SetFrameBorderMode(0);
    ROOT.gStyle.SetFrameBorderSize(1);
    ROOT.gStyle.SetFrameFillColor(0);
    ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineColor(1);
    ROOT.gStyle.SetFrameLineStyle(1);
    ROOT.gStyle.SetFrameLineWidth(1);
    ROOT.gStyle.SetPalette(55);
    ROOT.gStyle.SetNumberContours(100);

import numpy as np
def customPalette(zeropoint = 0.5):
	Number = 3;
	Red    = np.array([0  ,  100,  110],dtype=float)/255.
	Green  = np.array([0  ,  255,  0], dtype=float)/255.
	Blue   = np.array([99 ,  100,  2], dtype=float)/255.
	Length = np.array([0.0,  zeropoint, 1.0], dtype=float)
	nb=100;
	ROOT.TColor.CreateGradientColorTable(Number,Length,Red,Green,Blue,nb)

def drawMultipleGrid(hists,outname,limits=[],setLogY=False,setLogZ=False, ncols = 3, width=1500,height=1100):
	c = ROOT.TCanvas("c", "c", width,height)
	nhists = len(hists)
	nrows = (nhists-1)/ncols+1
	c.Divide(ncols,nrows)
	
	if len(limits) == 2:
		limits = [limits]*nhists
	if len(limits) == ncols and len(limits[0]) == 2:
		limits = limits*nrows

	for pad in range(len(hists)):
		p = c.cd(pad +1)
		if setLogY: p.SetLogy()
		if setLogZ: p.SetLogz()
		if limits: hists[pad].GetZaxis().SetRangeUser(limits[pad][0], limits[pad][1])
		hists[pad].Draw("colz")

	c.SaveAs(outname)

def saveHists(file,prefix="",filter=""):
    if prefix: prefix += '_'
    customROOTstyle()
    ROOT.gROOT.SetBatch(True)
    hists1d = ["TH1D", "TH1F", "TH1"]
    hists2d = ["TH2D", "TH2F", "TH2"]
    histObjectNames = hists1d + hists2d
    for key in file.GetListOfKeys():
        if key.IsFolder():
            dir = file.Get(key.GetName())
            saveHists(dir,prefix=prefix + key.GetName(), filter=filter)
        if key.GetClassName() in histObjectNames and filter in prefix:
            hist = file.Get(key.GetName())
            drawoptions = ""
            if key.GetClassName() in hists2d:
                drawoptions = "colz"
            drawHist(hist,prefix + key.GetName() + ".png", drawoptions = drawoptions)

def drawHist(hist,name,width=500,height=500, drawoptions=""):
    customROOTstyle()
    c = ROOT.TCanvas("c","c",width,height)
    #hist.SetLineWidth(2)
    hist.Draw(drawoptions)
    c.SaveAs(name)

def drawMultipleSame(hists,labels,filename,colors=[], width = 500, height = 500, norm = False, xtitle = "", ytitle = "", rebin = 0, leg="top",logy=False):
    customROOTstyle()
    hist_max = 0
    if not colors:
        colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
        colors = colors[:len(hists)]
    for h in hists:
        if rebin:
           h.Rebin(rebin)
        if norm:
           h.Scale(1./h.Integral())
        if h.GetMaximum() > hist_max:
            hist_max = h.GetMaximum()

    canv = ROOT.TCanvas("c","c",width,height)
    if logy:
    	 canv.SetLogy()
    first = True

    x1 = ROOT.gStyle.GetPadLeftMargin();
    x2 = 1 - ROOT.gStyle.GetPadRightMargin();
    if leg == "top":
       y2 = 1 - ROOT.gStyle.GetPadTopMargin();
       y1 = y2*.9
       y2 = ROOT.gStyle.SetPadTopMargin(y2);
    if leg == "bot":
       y1 = ROOT.gStyle.GetPadBottomMargin();
       y2 = y1 + 0.1
       ROOT.gStyle.SetPadBottomMargin(y1);

    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetNColumns(len(hists))
    for h,l,c in zip(hists,labels,colors):
        h.SetMaximum(1.2 * hist_max)
        h.SetTitle(l)
        h.SetLineColor(c)
        h.SetLineWidth(3)
        h.GetYaxis().SetTitleOffset(1.5)
        #h.SetOptStat(0)
        if first:
            if xtitle: h.GetXaxis().SetTitle(xtitle)
            if ytitle: h.GetYaxis().SetTitle(ytitle)
            h.Draw()
            first = False
        else:
            h.Draw("same")

        leg.AddEntry(h,l)


    leg.Draw()
    canv.SaveAs(filename)
