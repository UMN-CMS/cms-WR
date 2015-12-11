import ROOT
import numpy as np
from array import array
import sys

from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
from collections import defaultdict
customROOTstyle()

c = ROOT.TCanvas()
#c.SetLogz()
min_limit = 1000

crossections = {800:3.65, 1200:0.663, 1400:0.389, 1600:0.177, 2000:0.0707,
		2400:0.0248, 2600:0.015, 2800:0.00913, 3000:0.00576, 3200:0.0034,
		3600:0.00154, 3800:0.00119, 4000:0.000801, 4200:0.000473, 4400:0.000375,
		4600:0.00019, 4800:0.000152, 5000:0.0000912, 5200:0.0000665,
		5600:0.0000254, 5800:0.0000202, 6000:0.0000144,}
#args = "_".join(sys.argv[1:5])
#slope1 = float(sys.argv[1])
#intercept1 = float(sys.argv[2])
#slope2 = float(sys.argv[3])
#intercept2 = float(sys.argv[4])
#files = [open(name) for name in sys.argv[5:]]
files = [open(name) for name in sys.argv[1:]]

def dictoflists(): return {"low":[], "hi":[], "limit":[]}
hists = defaultdict(dictoflists)
for f in files:
	for line in f:
		l = line.split()
		if not l: continue
		MWR = l[0]
		MNR = l[1]
		try:
			MWRi = int(MWR)
		except ValueError:
			MWR = MWR[8:]

		low = int(float(l[2]))
		hi = int(float(l[3]))
		limit = float(l[-3])
		hists[(MWR,MNR)]["low"].append(low)
		hists[(MWR,MNR)]["hi"].append(hi)
		hists[(MWR,MNR)]["limit"].append(limit)

rms ={ 800:192.103227286, 1200:229.876231807, 1400:247.285490466,
		1600:269.256149708, 2000:302.149975396, 2400:352.353484885,
		2600:374.7755583, 2800:404.408429095, 3000:429.047814501,
		3200:452.764326089, 3600:508.20332794, 3800:540.205301701,
		4000:574.609631266, 4200:598.511463743, 4400:633.835317423,
		4600:668.892901662, 4800:705.95270766, 5000:750.166249925,
		5200:777.740102556, 5600:865.823245838, 5800:907.185734652,
		6000:946.564336157}

choice ={ 800:100, 1200:200, 1400:200, 1600:200, 2000:300, 2400:400, 2600:400,
		2800:500, 3000:500, 3200:600, 3600:700, 3800:700, 4000:900, 4200:1100,
		4400:1100, 4600:1100, 4800:1300, 5000:1500, 5200:1400, 5600:1700,
		5800:1900, 6000:2100}
choice ={ 800:(700,900), 1200:(1000,1400), 1400:(1200,1600), 1600:(1400,1800),
		2000:(1700,2300), 2400:(2000,2800), 2600:(2200,3000), 2800:(2300,3300),
		3000:(2500,3500), 3200:(2600,3800), 3600:(2900,4300), 3800:(3100,4500),
		4000:(3100,4900), 4200:(3100,5300), 4400:(3300,5500), 4600:(3500,5700),
		4800:(3500,6100), 5000:(3500,6500), 5200:(3800,6600), 5600:(3900,7300),
		5800:(3900,7700), 6000:(3900,8100)}
#choice ={ 800:(700,900), 1200:(1000,1400), 1400:(1200,1600), 1600:(1400,1800),
#		2000:(1700,2300), 2400:(2000,2800), 2600:(2200,3000), 2800:(2300,3300),
#		3000:(2500,3500), 3200:(2600,3800), 3600:(2900,4300), 3800:(3100,4500),
#		4000:(3100,4900), 4200:(3100,5300), 4400:(3300,5500), 4600:(3500,5700),
#		4800:(3500,6100), 5000:(3500,6500), 5200:(3800,6600), 5600:(3900,7300),
#		5800:(3900,7700), 6000:(3900,8100)}

#for mass in choice:
#	low = round(slope1*mass + intercept1, -2)
#	hi = round(slope2*mass + intercept2, -2)
#	hi = 2*mass - low
#	if mass > 4600:
#		low = round(slope1*4600 + intercept1, -2)
#	choice[mass] = (low, hi)
#
#cs = choice.keys()
#cs.sort()
#for ci in cs:
#	print ci, choice[ci]

results_opt = []
results_min = []

for MWR,MNR in hists:
	h = hists[(MWR,MNR)]
	#arg = np.argmin(h["limit"])
	#print MWR,MNR, h["low"][arg], h["hi"][arg], h["limit"][arg]
	#low = int(MWR) - round(choice[int(MWR)],-2)
	#hi = int(MWR) + round(choice[int(MWR)],-2)
#	low, hi = choice[int(MWR)]
#	print MWR,low,hi
#	for i in range(len(h["low"])):
#		if low == h["low"][i] and hi == h["hi"][i]:
#			arg = i
#			break
#	print i,arg
#	results_opt.append( [int(MWR),  h["low"][arg], h["hi"][arg], h["limit"][arg] ])
#
	MWRi = float(MWR)
	arg = np.argmin(h["limit"])
	results_min.append( [MWRi,  h["low"][arg], h["hi"][arg], h["limit"][arg] ])
	print MWR,MNR, h["low"][arg], h["hi"][arg], h["limit"][arg]

	min_low = min(h["low"])
	max_low = max(h["low"])
	min_hi = min(h["hi"])
	max_hi = max(h["hi"])
	n_low = len(set(h["low"]))
	n_hi = len(set(h["hi"]))
	name = "WR_M_%s_Nu_M_%s" % (MWR, MNR)
	binsize = int(.02 * MWRi)
	binsize = .02
	#hist = ROOT.TH2F(name, name,n_low, min_low/MWRi,(max_low+binsize)/MWRi,n_hi,(min_hi)/MWRi,(max_hi+binsize)/MWRi)
	#hist = ROOT.TH2F(name, name, n_low, min_low/MWRi,(max_low+binsize)/MWRi,n_hi,(min_hi)/MWRi,(max_hi+binsize)/MWRi)
	hist = ROOT.TH2F(name, name, 31, .59-.2, 1.01, 21, .99, 1.41)
	hist.SetXTitle("Min")
	hist.SetYTitle("Max")
	#hist.GetZaxis().SetMoreLogLabels()
	hist.SetMinimum(min(h["limit"]) - 1e-5)
	hist.SetMaximum(min(h["limit"]) * 1.3)
	for i in range(len(h["low"])):
		hist.Fill(h["low"][i]/MWRi, h["hi"][i]/MWRi, h["limit"][i])
	hist.Draw("colz")
	c.SaveAs("plots/mass_window_opt2/window_plot" + name + ".png")

#results_opt.sort()
results_min.sort()
print results_min
#print results_opt

#for res,name in [(results_opt,"opt"), (results_min, "min")]:
for res,name in [(results_min, "min")]:
	npresults = np.array(res, dtype="float64")
	min_mass  = np.copy(npresults.T[0])
	min_low   = np.copy(npresults.T[1])/min_mass
	min_hi    = np.copy(npresults.T[2])/min_mass
	min_limit = np.copy(npresults.T[3])*.01
	npcrosssection = np.array( [ crossections[m] for m in min_mass], dtype="float64")
	
	lowmass = min(min_mass)
	himass = max(min_mass)
	h = ROOT.TH1F("h","Mass Window",1,lowmass,himass)
	h.Fill(np.average((lowmass,himass)))
	h.SetMinimum(.5)
	h.SetMaximum(1.5)

	low_graph = ROOT.TGraph(len(min_mass),  min_mass, min_low)
	hi_graph = ROOT.TGraph(len(min_mass), min_mass, min_hi)
	limit_graph = ROOT.TGraph(len(min_mass), min_mass, min_limit)
	xs_graph = ROOT.TGraph(len(min_mass), min_mass, npcrosssection)
	
	low_graph.SetLineColor(ROOT.kRed)
	h.SetLineColor(ROOT.kBlack)
	hi_graph.SetLineColor(ROOT.kBlue)
	limit_graph.SetLineColor(ROOT.kBlack)
	xs_graph.SetLineColor(ROOT.kRed)
	low_graph.SetLineWidth(2)
	h.SetLineWidth(2)
	hi_graph.SetLineWidth(2)
	limit_graph.SetLineWidth(2)
	xs_graph.SetLineWidth(2)
	
	leg1 = ROOT.TLegend(.15, .15, .5, .25)
	leg1.SetNColumns(3)
	leg1.AddEntry(low_graph, "Low/Mass", "l")
	leg1.AddEntry(hi_graph, "High/Mass", "l")
	leg1.AddEntry(h, "Mass", "l")
	leg2 = ROOT.TLegend(.55, .75, .9, .85)
	leg2.SetNColumns(2)
	leg2.AddEntry(limit_graph, "XS limit", "l")
	leg2.AddEntry(xs_graph, "WR XS", "l")
	c = ROOT.TCanvas("c","c",1200,800)
	c.Divide(1,2)

	c.cd(1)
	#hi_graph.SetTitle("Mass Window")
	h.Draw()
	hi_graph.Draw("PLsame")
	low_graph.Draw("PLsame")
	#low_graph.Fit("pol1","","",800,4000)
	#hi_graph.Fit("pol1","","",800,4000)
	leg1.Draw()
	
	c.cd(2)
	limit_graph.SetTitle("Resulting limit (if XS = .01)")
	limit_graph.Draw("APL")
	xs_graph.Draw("PLsame")
	leg2.Draw()
	#c.SaveAs("mass_window_" + args + "_" + name + ".png")
	c.SaveAs("plots/mass_window_opt2/mass_window_" + name + ".png")
 
#npresmin = np.array(results_min, dtype="float64")
#npresopt = np.array(results_opt, dtype="float64")
#min_mass  = np.copy(npresmin.T[0])
#
#c = ROOT.TCanvas()
#h = ROOT.TH1F("h","diff (opt - min)/min",1,lowmass,himass)
#h.SetMaximum(.2)
#h.SetMinimum(0)
#h.Draw()
#limit_diff = np.copy((-npresmin.T[3] + npresopt.T[3])/npresmin.T[3])
#graph = ROOT.TGraph(len(min_mass),  min_mass, limit_diff)
#graph.Draw("PLsame")
#c.SaveAs("plots/mass_window_opt2mass_window_diff_" + args + ".png")

