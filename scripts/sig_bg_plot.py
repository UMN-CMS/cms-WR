import ROOT
import itertools
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
import math

customROOTstyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHatchesLineWidth(2)
colors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
fillstyle = [1001, 3345, 3354, 3006]

col = ["Signal", "TT", "DY"]

ele_names = [a + ' ' + b for a,b in itertools.product(["ee"], col)]
mu_names = [a + ' ' + b for a,b in itertools.product(["#mu#mu"], col)]
print ele_names, mu_names

mass = np.ndarray(25)
ele_results = np.ndarray((3,25))
mu_results = np.ndarray((3,25))
ele_sigma = np.ndarray((3,25))
mu_sigma = np.ndarray((3,25))
ele_linei=0
mu_linei=0
with open("datacards/summary.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		line = line.split()
		CHANNEL = line[0]
		MASS, XS,SIG_NEVENTS, SIG_RATE, TT_RATE, DY_RATE, TT_SF, DY_SF, SIG_N, TT_N, DY_N, SIG_ALPHA, TT_ALPHA, DY_ALPHA = map(float, line[1:])
		if CHANNEL == "ee":
			mass[ele_linei] = int(MASS)
			ele_results[:,ele_linei] = SIG_NEVENTS, (TT_N+1)*TT_ALPHA, (DY_N+1)*DY_ALPHA
			ele_sigma[:,ele_linei] = math.sqrt(SIG_N + 1) * SIG_ALPHA, math.sqrt(TT_N + 1)  * TT_ALPHA, math.sqrt(DY_N + 1)  * DY_ALPHA

			ele_linei += 1
		elif CHANNEL ==  "mumu":
			mu_results[:,mu_linei] = SIG_NEVENTS, (TT_N+1)*TT_ALPHA, (DY_N+1)*DY_ALPHA
			mu_sigma[:,mu_linei] = math.sqrt(SIG_N + 1) * SIG_ALPHA, math.sqrt(TT_N + 1)  * TT_ALPHA, math.sqrt(DY_N + 1)  * DY_ALPHA

			mu_linei += 1

h = ROOT.TH1F("h","",1,800, 6000)
h.GetYaxis().SetTitleOffset(1.3)
h.SetMinimum(.001)


def draw(ch, res, sigma, names):
	c = ROOT.TCanvas("c","c",800,800)
	c.SetLeftMargin(.15)
	c.SetRightMargin(.05)
	c.SetLogy()
	h.Draw()
	h.SetXTitle("W_R Mass [GeV]")
	h.SetYTitle("Events")
	h.SetMaximum(res.max()*2)

	inset = ROOT.TPad("inset","",.57, .57, .93, .93)
	inset.Draw()
	inset.cd()
	ins = ROOT.TH1F("ins","",1,700, 6000)
	ins.SetMaximum(10)
	ins.Draw()

	c.cd()

	leg = ROOT.TLegend(.3, .8, .5, .94)
#leg.SetNColumns(3)
	leg.SetTextFont(42)
	leg.SetTextSize(0.032)
	leg.SetBorderSize(0)
	leg.SetFillStyle(0)

	graphs = []
	for n,(col,err) in enumerate(zip(res,sigma)):
		graphs.append( ROOT.TGraphErrors(len(mass), mass, col,ROOT.nullptr,err))
		g = graphs[-1]
		g.SetFillColor(colors[n] - 9)
		g.SetFillStyle(fillstyle[n])
		g.SetLineColor(colors[n])
		g.SetLineWidth(3)
		c.cd()
		g.Draw("3same")
		g.Draw("LXsame")
		leg.AddEntry(g, names[n], "LF")
		inset.cd()
		g.Draw("LXsame")

	c.cd()
	leg.Draw("same")
	c.RedrawAxis()
	inset.RedrawAxis()
	c.SaveAs("plots/%s_events.png" % ch)
	c.SaveAs("plots/%s_events.pdf" % ch)
	c.SaveAs("plots/%s_events.C" % ch)

draw("ele", ele_results, ele_sigma, ele_names)
draw("mu", mu_results, mu_sigma, mu_names)
