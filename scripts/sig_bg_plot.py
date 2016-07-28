import ROOT
import itertools
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle

customROOTstyle()
ROOT.gROOT.SetBatch(True)
colors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kBlack]

col = ["Signal", "TT", "DY"]

ele_names = [a + ' ' + b for a,b in itertools.product(["ee"], col)]
mu_names = [a + ' ' + b for a,b in itertools.product(["#mu#mu"], col)]
print ele_names, mu_names

mass = np.ndarray(25)
ele_results = np.ndarray((3,25))
mu_results = np.ndarray((3,25))
linei=0
with open("datacards/sig_bg.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		line = map(float,line.split())
		mass[linei] = line[0]

		ele_results[:,linei] = line[1:4]
		mu_results[:,linei]   = line[4:]
		linei += 1

h = ROOT.TH1F("h","",1,700, 6200)
h.GetYaxis().SetTitleOffset(1.3)
h.SetMinimum(.0001)


#electrons
c = ROOT.TCanvas("c","c",800,800)
c.SetLeftMargin(.15)
c.SetLogy()
h.Draw()
h.SetXTitle("W_R Mass [GeV]")
h.SetYTitle("Events")
h.SetMaximum(ele_results.max()*2)

inset = ROOT.TPad("inset","",.17, .15, .55, .55)
inset.Draw()
inset.cd()
ins = ROOT.TH1F("ins","",1,700, 6200)
ins.SetMaximum(10)
ins.Draw()

c.cd()

leg = ROOT.TLegend(.5, .7, .8, .9)
#leg.SetNColumns(3)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
for n,col in enumerate(ele_results):
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	c.cd()
	g.Draw("Lsame")
	leg.AddEntry(g, ele_names[n], "L")
	inset.cd()
	g.Draw("Lsame")

c.cd()
leg.Draw("same")
c.SaveAs("plots/ele_events.png")
c.SaveAs("plots/ele_events.pdf")
c.SaveAs("plots/ele_events.C")

#muons
c = ROOT.TCanvas("c","c",800,800)
c.SetLeftMargin(.15)
c.SetLogy()
h.Draw()
h.SetXTitle("W_R Mass [GeV]")
h.SetYTitle("Events")
h.SetMaximum(mu_results.max()*2)

inset = ROOT.TPad("inset","",.17, .15, .55, .55)
inset.Draw()
inset.cd()
ins = ROOT.TH1F("ins","",1,700, 6200)
ins.SetMaximum(10)
ins.Draw()

c.cd()

leg = ROOT.TLegend(.5, .7, .80, .9)
#leg.SetNColumns(3)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
for n,col in enumerate(mu_results):
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	c.cd()
	g.Draw("Lsame")
	leg.AddEntry(g, mu_names[n], "L")
	inset.cd()
	g.Draw("Lsame")

c.cd()
leg.Draw("same")
c.SaveAs("plots/mu_events.png")
c.SaveAs("plots/mu_events.pdf")
c.SaveAs("plots/mu_events.C")
