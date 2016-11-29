from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
from itertools import product
from collections import defaultdict
import numpy as np
import csv


def colors():
	return None
customROOTstyle()

csv.register_dialect("space", delimiter=" ", skipinitialspace=True)

chs = "ee","mumu"
process = "signal","DYAMC","TT"

chs_n = "ee", "#mu#mu"
process_n = "Signal", "DY", "T#bar{T}"

#masses = np.array( [ 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800,
#		3000, 3200, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5600,
#		5800, 6000 ], dtype=float)
masses = np.array( [ 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800,
		3000, 3200, 3600, 3800, 4000], dtype=float)
n_masses = len(masses)

keys = [ "_".join([a,b]) for  b,a in product(process,chs)]
#make a dictionary mapping keys to graph names
names  = dict( zip( keys, [ " ".join([a,b]) for  b,a in product(process_n,chs_n)] ))
print names
means = defaultdict(lambda: np.zeros((n_masses)))
stats = defaultdict(lambda: np.zeros((n_masses)))
systs = defaultdict(lambda: np.zeros((n_masses)))
jet_systs = defaultdict(lambda: np.zeros((n_masses)))
lep_systs = defaultdict(lambda: np.zeros((n_masses)))

# store column in result
def get_column(filename, column, result):
	with open(filename,'r') as f:
		reader = csv.DictReader(f,dialect="space")
		for row in reader:
			#strip take first two fields that are separted by "_"
			key = "_".join(row['key'].split("_")[:2])
			mass = int(row["mass"])
			if mass > 4000:
				continue
			val = float(row[column])
			mass_i = np.where(masses==mass)[0][0]
			result[key][mass_i] = val

get_column("all_unc.txt", "mean", means)
get_column("all_unc.txt", "syst_std", systs)
get_column("all_unc.txt", "syst_stat_err", stats)

get_column("jet_unc.txt", "syst_std", jet_systs)
get_column("ele_unc.txt", "syst_std", lep_systs)
get_column("mu_unc.txt",  "syst_std", lep_systs)
get_column("emu_unc.txt", "syst_std", lep_systs)

for unc in [systs, stats, jet_systs, lep_systs]:
	for key in unc:
		unc[key] /= means[key]/100.

import ROOT

c = ROOT.TCanvas("c","c",600,600)
c.SetTopMargin(.15)
c.SetLeftMargin(.15)
c.SetRightMargin(0.05)
c.SetGridx()
c.SetGridy()
c.SetTicks(1,1)

h = ROOT.TH1F("h","",1,750,4050)
h.SetMinimum(0)
ROOT.TGaxis.SetMaxDigits(3)
h.GetYaxis().SetTitleOffset(1.3)
h.SetXTitle("W_{R} Mass Hypothesis")
h.SetYTitle("Relative Uncertainty")

unc_names = ["Stat.", "Total Syst.", "Jet Syst.", "Lepton Syst"]
unc_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1]  
unc_markers = [21, 20, 24, 25] 
unc_plot_info = zip(unc_names, unc_colors, unc_markers)
uncs = [stats, systs, jet_systs, lep_systs]


graphs = []
for key in keys:

	h.SetYTitle(names[key] + " Relative Uncertainty[%]")
	h.SetMaximum(2.5 if "signal" in key else 100)
	h.Draw()
	leg = ROOT.TLegend(.15,.9,.95,.95)
	leg.SetNColumns(len(uncs))
	leg.SetBorderSize(0)
	leg.SetFillStyle(0)
	for (name,color,marker), unc in zip(unc_plot_info, uncs):
		g = ROOT.TGraph(n_masses, masses, unc[key])
		g.SetMarkerColor(color)
		g.SetMarkerStyle(marker)
		g.SetLineWidth(2)
		graphs.append(g)
		g.Draw("Psame")
		leg.AddEntry(g,name,"p")

	leg.Draw()
	c.SaveAs("plots/unc_" + key + ".png")
