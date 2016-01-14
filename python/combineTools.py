import ROOT as r
import re
import numpy as np
import subprocess
from random import gauss
import ExoAnalysis.cmsWR.backgroundFits as bgfits

##
# @brief creates a datacard for combine for signal and background
#
# @param outfile String of filename
# @param bin_name Name for bin
# @param nObs number of observed
# @param signal_tuple (sig_name, signal_rate)
# @param background_tuples list of (name, rate) for backgrounds
#
# @return signal_rate, tuple of bg rates
def makeDataCardSingleBin(outfile, bin_name, nObs, signal_tuple, background_tuples, systematics = []):

	nBGs = len(background_tuples)
	nSystematics = len(systematics)
	ns = "  ".join([str(i) for i in range(nBGs+1)])
	sig_name = signal_tuple[0]
	signal_rate = "%.4g" % (signal_tuple[1])

	bg_names = ""
	bg_rates = ""
	name_lookup = {sig_name:0}
	for i,(bg_name, bg_rate) in enumerate(background_tuples):
		name_lookup[bg_name] = i+1
		bg_names += bg_name + "  "
		bg_rates   += "%.4g" % bg_rate + "  "
	names = sig_name + "  " + bg_names
	rates = signal_rate + '  ' + bg_rates
	with open(outfile, "w") as out:
		out.write("imax 1  number of channels\n")
		out.write("jmax %d  number of backgrounds\n" % nBGs)
		if systematics:
			out.write("kmax %d  number of nuisance parameters\n" % nSystematics)
		out.write("bin " + bin_name + "\n")
		out.write("observation %d\n" % nObs)
		out.write("------------\n")
		out.write("bin" + ("    " + bin_name)* (nBGs + 1) + "\n")
		out.write("process  " + names + "\n")
		out.write("process  " + ns + "\n")
		out.write("rate  " + rates + "\n")
		out.write("------------\n")
		if systematics:
			for name, syst_type, channels in systematics:
				channel_list = ['-']*(nBGs + 1)
				for chan_name, err in channels:
					channel_list[name_lookup[chan_name]] = str(err)
				out.write(name + ' ' + syst_type + ' ' + ' '.join(channel_list) + '\n')
	return signal_rate, tuple(bg_rates.split())

##
# @brief gets fourObjectMass histogram from tree in root file
#
# @param dirname name for TDirectory that contains TTree
# @param treename name of TTree
# @param process tuple of (histname, filename). histname is name of output histogram. filename which contains TTree
# @param binsize
#
# @return 
def getMassHisto(dirname, treename, process, binsize=100):
	histname, filename = process
	infile = r.TFile.Open(filename)
	tree = infile.Get(dirname + '/' + treename)
	r.gROOT.cd()
	nbins = 10000/int(binsize)
	h = r.TH1F(histname, histname, nbins, 0, nbins * binsize)
	tree.Draw("fourObjectMass>>" + histname, "" ,"goff")
	return h

##
# @brief calls and parses command for `combine'
#
# @param command a list to pass to subprocess.check_output
#
# @return (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)
def runCombine(command):
	try:
		output = subprocess.check_output(command)
		p = re.compile(r'mean   expected limit: r < ([0-9.]*) \+/- ([0-9.]*)')
		mean,meanError = p.search(output).groups()
		p = re.compile(r'68% expected band : ([0-9.]*) < r < ([0-9.]*)')
		onesig_minus,onesig_plus = p.search(output).groups()
		p = re.compile(r'95% expected band : ([0-9.]*) < r < ([0-9.]*)')
		twosig_minus,twosig_plus = p.search(output).groups()
		return (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)
	except:
		print output

#def getHists(dirname, treename, processes):
#	 return [getMassHisto(dirname,  treename,  p) for p in processes]

