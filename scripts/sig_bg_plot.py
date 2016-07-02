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
with open("sig_bg.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		line = map(float,line.split())
		mass[linei] = line[0]

		ele_results[:,linei] = line[1:4]
		mu_results[:,linei]   = line[4:]
		linei += 1

h = ROOT.TH1F("h","",1,700, 6200)

h.SetMinimum(.0001)
h.SetMaximum(ele_results.max()*2)
h.GetYaxis().SetTitleOffset(1.2)

#electrons
c = ROOT.TCanvas("c","c",800,800)
c.SetLeftMargin(1.5)
c.SetLogy()
h.Draw()
h.SetXTitle("W_R Mass [GeV]")
h.SetYTitle("Events")

leg = ROOT.TLegend(.5, .7, .85, .9)
#leg.SetNColumns(3)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
for n,col in enumerate(ele_results):
	print ele_names[n]
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	g.Draw("Lsame")
	leg.AddEntry(g, ele_names[n], "L")

leg.Draw("same")
c.SaveAs("plots/ele_events.png")
c.SaveAs("plots/ele_events.pdf")
c.SaveAs("plots/ele_events.C")
c.SetLogy(0)
h.SetMaximum(ele_results.max()*1.1)
h.SetMaximum(10)
c.SaveAs("plots/ele_events_lin.png")
c.SaveAs("plots/ele_events_lin.pdf")
c.SaveAs("plots/ele_events_lin.C")


h.SetMaximum(mu_results.max()*2)
h.Draw()

leg = ROOT.TLegend(.5, .7, .85, .9)
#leg.SetNColumns(3)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

for n,col in enumerate(mu_results):
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	g.Draw("Lsame")
	leg.AddEntry(g, mu_names[n], "L")

leg.Draw("same")
c.SetLogy()
c.SaveAs("plots/mu_events.png")
c.SaveAs("plots/mu_events.pdf")
c.SaveAs("plots/mu_events.C")
c.SetLogy(0)
h.SetMaximum(mu_results.max()*1.1)
h.SetMaximum(10)
c.SaveAs("plots/mu_events_lin.png")
c.SaveAs("plots/mu_events_lin.pdf")
c.SaveAs("plots/mu_events_lin.C")
