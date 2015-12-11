import ROOT as r
import re
import numpy as np
import subprocess
from random import gauss
import ExoAnalysis.cmsWR.backgroundFits as bgfits

#def getRates(histos, norms, lumi, mass_range):
#	if len(histos) != len(norms):
#		print "LENGTH MISMATCH"
#		return
#	low,hi = mass_range
#	ret = []
#	for i in range(len(histos)):
#		lowbin = histos[i].FindBin(low)
#		hibin = histos[i].FindBin(hi)
#		ret.append( (histos[i].GetName(), histos[i].Integrate(lowbin,hibin)*norms[i]*lumi) )
#	return [ (histos[i].GetName(), histos[i].GetEntries()*norms[i]*lumi) for i in range(len(histos))]
#
#def sumRateTuple(rates):
#	return sum([r[1] for r in rates])

def makeDataCardSingleBin(outfile, bin_name, nObs, signal_tuple, background_tuples):

	nBGs = len(background_tuples)
	ns = "  ".join([str(i) for i in range(nBGs+1)])
	sig_name = signal_tuple[0]
	signal_rate = "%.2f" % (signal_tuple[1])

	bg_names = ""
	bg_rates = ""
	for bg_name, bg_rate in background_tuples:
		bg_names += bg_name + "  "
		bg_rates   += "%.2e" % bg_rate + "  "
	names = sig_name + "  " + bg_names
	rates = signal_rate + '  ' + bg_rates
	with open(outfile, "w") as out:
		out.write("imax 1  number of channels\n")
		out.write("jmax %d  number of backgrounds\n" % nBGs)
		out.write("bin " + bin_name + "\n")
		out.write("observation %d\n" % nObs)
		out.write("------------\n")
		out.write("bin" + ("    " + bin_name)* (nBGs + 1) + "\n")
		out.write("process  " + names + "\n")
		out.write("process  " + ns + "\n")
		out.write("rate  " + rates + "\n")
		out.write("------------\n")
	return signal_rate, tuple(bg_rates.split())

def getMassHisto(dirname, treename, process, binsize=100):
	histname, filename = process
	infile = r.TFile.Open(filename)
	tree = infile.Get(dirname + '/' + treename)
	r.gROOT.cd()
	nbins = 10000/int(binsize)
	h = r.TH1F(histname, histname, nbins, 0, nbins * binsize)
	tree.Draw("fourObjectMass>>" + histname, "" ,"goff")
	return h

def runCombine(command):
		output = subprocess.check_output(command)
		p = re.compile(r'mean   expected limit: r < ([0-9.]*) \+/- ([0-9.]*)')
		mean,meanError = p.search(output).groups()
		p = re.compile(r'68% expected band : ([0-9.]*) < r < ([0-9.]*)')
		onesig_minus,onesig_plus = p.search(output).groups()
		p = re.compile(r'95% expected band : ([0-9.]*) < r < ([0-9.]*)')
		twosig_minus,twosig_plus = p.search(output).groups()
		return (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)

#def getHists(dirname, treename, processes):
#	 return [getMassHisto(dirname,  treename,  p) for p in processes]

WR_eejj_norm = 0.0707/50000
#DYJets_eejj_norm = 6025.2/9052671
#TTonly_eejj_norm = 831.76/96834559
#WZ_eejj_norm = 66.1/991232
#ZZ_eejj_norm = 15.4/996168
#bg_norms = [DYJets_eejj_norm, TTonly_eejj_norm, WZ_eejj_norm, ZZ_eejj_norm]

#processes = [("DYJets_eejj", "data/analyzed_DYJets_Madgraph_M_50_25ns_eejj_signal_region.root"),
#		("TTOnly_eejj", "data/analyzed_TTOnly_PowhegPythia_25ns_eejj_signal_region_reMiniAOD.root"),
#		("WZ_eejj", "data/analyzed_WZ_25ns_eejj_signal_region.root"),
#		("ZZ_eejj", "data/analyzed_ZZ_25ns_eejj_signal_region.root")
#		]

#backgrounds = getHists(dirname, treename, processes)



#lumi = 1900
#print "#lumi\tbg\tsig\tobs\tobs/bg\tMethod\tr(sig/expect sig)"
#bg_tuples = getRates(backgrounds, bg_norms, lumi)
#signal_tuple = getRates([signal], [WR_eejj_norm], lumi)[0]
#total_bg = sumRateTuple(bg_tuples)
#for obs in np.linspace(total_bg, (total_bg + signal_tuple[1])*1.2, 10, dtype=int):
#	data = r.TH1F("data_eejj","data_eejj",100,600,4600)
#	[data.Fill(gauss(2000, 500)) for i in range(obs)]
#
#	datacard = "WReejj_LUMI%d_OBS%d" % (lumi, obs)
#	datacard_file = "data/" + datacard + ".txt"
#	makeDataCardSingleBin(datacard_file, "eejj", data, signal_tuple, bg_tuples)
#	#for method in ["BayesianSimple", "HybridNew"]:
#	for method in ["BayesianSimple"]:
#		output = subprocess.check_output(["combine","-M",method,"-S","0",datacard_file,"-n", datacard])
#		combinefile = r.TFile.Open("higgsCombine" + datacard + "." + method + ".mH120.root")
#		tree = combinefile.Get("limit")
#		tree.GetEntry(0)
#		print "%d\t%d\t%d\t%d\t%.3f\t%s\t%.4f" %(lumi, int(total_bg), signal_tuple[1],
#				obs, obs/float(total_bg), method, tree.limit)
#		subprocess.check_output(["mv", "higgsCombine" + datacard + "." + method + ".mH120.root", "data/"])
	

