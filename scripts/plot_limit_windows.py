import ROOT

from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
from collections import defaultdict
customROOTstyle()

c = ROOT.TCanvas()
c.SetLogz()
min_limit = 1000

f = open("windows-all.txt")

def dictoflists(): return {"low":[], "hi":[], "limit":[]}
hists = defaultdict(dictoflists)
for line in f:
   l = line.split()
   if not l: continue
   MWR = l[0]
   MNR = l[1]
   low = int(l[2])
   hi = int(l[3])
   limit = float(l[-3])
   hists[(MWR,MNR)]["low"].append(low)
   hists[(MWR,MNR)]["hi"].append(hi)
   hists[(MWR,MNR)]["limit"].append(limit)

for MWR,MNR in hists:
	print MWR,MNR
	print 
	h = hists[(MWR,MNR)]
	min_low = min(h["low"])
	max_low = max(h["low"])
	min_hi = min(h["hi"])
	max_hi = max(h["hi"])
	n_low = len(set(h["low"]))
	n_hi = len(set(h["hi"]))
	name = "WR_M_%s_Nu_M_%s" % (MWR, MNR)
	hist = ROOT.TH2F(name, name,n_low,min_low,max_low+100,n_hi,min_hi,max_hi+100)
	hist.SetXTitle("Min")
	hist.SetYTitle("Max")
	hist.GetZaxis().SetMoreLogLabels()
	hist.SetMinimum(min(h["limit"]) - 1e-5)
	for i in range(len(h["low"])):
		hist.Fill(h["low"][i], h["hi"][i], h["limit"][i])
	hist.Draw("colz")
	c.SaveAs("window_plot" + name + ".png")
 

