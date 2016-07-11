import ROOT
import itertools
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle

customROOTstyle()
ROOT.gROOT.SetBatch(True)
colors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kBlack]

col = ["BB", "EB", "EE", "Global"]

ele_names = [a + ' ' + b for a,b in itertools.product(["Electron"], col)]
mu_names = [a + ' ' + b for a,b in itertools.product(["Muon"], col)]
mass = np.ndarray(25)
ele_results = np.ndarray((4,25))
ele_nt_results = np.ndarray((4,25))
mu_results = np.ndarray((4,25))
mu_nt_results = np.ndarray((4,25))
linei=0
with open("lepton_eff.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		line = map(float,line.split())
		mass[linei] = line[0]

		ele_results[:,linei] = line[1:5]
		mu_results[:,linei] = line[5:]
		linei += 1

linei=0
with open("lepton_eff_notrigger.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		line = map(float,line.split())

		print line
		ele_nt_results[:,linei] = line[1:5]
		mu_nt_results[:,linei] = line[5:]
		linei += 1


h = ROOT.TH1F("h","",1,700, 6200)

h.SetMinimum(.95)
h.SetMaximum(1.02)

#electrons
c = ROOT.TCanvas("c","c",800,800)
c.SetLeftMargin(.2)
h.Draw()
h.SetXTitle("W_{R} Mass [GeV]")
h.SetYTitle("Lepton Eff/Lepton Eff (notrigger)")
h.SetTitleOffset(1.5,"Y")

leg = ROOT.TLegend(.3, .8, .85, .9)
leg.SetNColumns(2)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
for n,col in enumerate(ele_results/ele_nt_results):
	print ele_names[n]
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	g.Draw("Lsame")
	leg.AddEntry(g, ele_names[n], "L")

leg.Draw("same")
c.SaveAs("plots/trig_ele_eff.png")
c.SaveAs("plots/trig_ele_eff.pdf")
c.SaveAs("plots/trig_ele_eff.C")

h.Draw()

#muons
leg = ROOT.TLegend(.3, .8, .85, .9)
leg.SetNColumns(2)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
for n,col in enumerate(mu_results/mu_nt_results):
	print mu_names[n]
	graphs.append( ROOT.TGraph(len(mass), mass, col))
	g = graphs[-1]
	g.SetLineColor(colors[n])
	g.SetLineWidth(3)
	g.Draw("Lsame")
	leg.AddEntry(g, mu_names[n], "L")

leg.Draw("same")
c.SaveAs("plots/trig_mu_eff.png")
c.SaveAs("plots/trig_mu_eff.pdf")
c.SaveAs("plots/trig_mu_eff.C")

h.Draw()
