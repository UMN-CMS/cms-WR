import ROOT
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle

customROOTstyle()
ROOT.gROOT.SetBatch(True)

#declare numpy arrays which hold WR mass points, and the efficiency of WR events at each mass point to pass the entire online and offline selection including mass window cuts
#mass = np.ndarray(25)
#ele_results = np.ndarray(25)
#mu_results = np.ndarray(25)
mass = np.ndarray(2)
ele_results = np.ndarray(2)
mu_results = np.ndarray(2)
linei=0

with open("effWithWindowCuts.txt","r") as f:
	for line in f:
		if line[0] == '#': continue
		print line
		line = map(float,line.split())
		mass[linei] = line[0]
		ele_results[linei] = line[1]
		mu_results[linei] = line[2]
		linei += 1



print 'mass list contains', mass
print 'ele_results list contains', ele_results
print 'mu_results list contains', mu_results
#h = ROOT.TH1F("h","",1,700, 6200)
#
#
#h.SetMinimum(0)
#h.SetMaximum(1)
#
###electrons
#c = ROOT.TCanvas("c","c",800,800)
#h.Draw()
#h.SetXTitle("M_{LLJJ} [GeV]")
#h.SetYTitle("Full Selection Eff*Accept")
#
#leg = ROOT.TLegend(.3, .8, .85, .9)
##delete if plot is ok    leg.SetNColumns(2)
#leg.SetTextFont(42)
#leg.SetTextSize(0.032)
#leg.SetBorderSize(0)
#
#graphs = []
#graphs.append( ROOT.TGraph(len(mass), mass, ele_results))
#g = graphs[-1]
#g.SetLineColor(ROOT.kBlack)
#g.SetLineWidth(3)
#g.Draw("L")
#leg.AddEntry(g, "Global", "L")
#
#leg.Draw("same")
#c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.png")
#c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.pdf")
#c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.C")
#
#
#h.Draw()
#
###muons
#leg = ROOT.TLegend(.3, .8, .85, .9)
##leg.SetNColumns(2)
#leg.SetTextFont(42)
#leg.SetTextSize(0.032)
#leg.SetBorderSize(0)
#
#graphs.append( ROOT.TGraph(len(mass), mass, mu_results))
#g = graphs[-1]
#g.SetLineColor(ROOT.kBlack)
#g.SetLineWidth(3)
#g.Draw("L")
#leg.AddEntry(g, "Global", "L")
#
#
#leg.Draw("same")
#c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.png")
#c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.pdf")
#c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.C")
