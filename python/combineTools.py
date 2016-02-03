import ROOT as r
import re
import numpy as np
import subprocess
from random import gauss
import ExoAnalysis.cmsWR.backgroundFits as bgfits
import ExoAnalysis.cmsWR.cross_sections as xs
import math

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
		out.write("observation %.4g\n" % nObs)
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
# @param MWR Mass of WR
# @param binsize
#
# @return TH1F of fourobjectmass
def getEEJJMassHisto(MWR, binsize=100):
	filename = "data/all_analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-%d_Nu_M-%d.root" % (MWR, MWR/2)
	dirname =  "unmatchedSignalRecoAnalyzerFive"
	treename = "signalRecoObjectsWithAllCuts"
	histname = "WR_eejj_" + str(MWR)
	infile = r.TFile.Open(filename)
	tree = infile.Get(dirname + '/' + treename)
	r.gROOT.cd()
	nbins = 10000/int(binsize)
	h = r.TH1F(histname, histname, nbins, 0, nbins * binsize)
	tree.Draw("fourObjectMass>>" + histname, "" ,"goff")
	nevents = 50000.
	h.Scale( xs.WR_eejj[MWR]/nevents )
	return h

def getMuMuJJMassHisto(MWR, binsize=100):
	filename = "data/mumu_histos.root"
	infile = r.TFile.Open(filename)
	r.gROOT.cd()
	h = infile.Get("h_" + str(MWR)).Clone()
	orig_lumi = 2460.
	h.Scale((xs.WR_mumujj[MWR]/0.01)*(1./orig_lumi))
	return h

def getSignalMassHisto(MWR, channel, binsize=100):
	if "ee" in channel or "EE" in channel:
		return getEEJJMassHisto(MWR,binsize=binsize)
	else:
		return getMuMuJJMassHisto(MWR,binsize=binsize)

def getTTBarNEvents(MWR, channel, lumi):
	if "ee" in channel or "EE" in channel:
		filename = "data/ttBarBkgndEleEstimate.root"
		workname = "ttBarElectronEstimate"
		corr_name = "EEtoEMuCorrection"
	else:
		filename = "data/ttBarBkgndMuonEstimate.root"
		workname = "ttBarMuonEstimate"
		corr_name = "MuMutoEMuCorrection"
	f = r.TFile.Open(filename)
	work = f.Get(workname)
	pdf = work.pdf("rescaledExpPdf")
	mass = work.var("fourObjectMass")
	data = work.data("eMuRealDataSet")

	mass.setRange(str(MWR), mass_cut[MWR][0], mass_cut[MWR][1])
	argset = r.RooArgSet(mass)
	integral = pdf.createIntegral(argset, r.RooFit.NormSet(argset), r.RooFit.Range(str(MWR)))

	real_norm = work.var("realEMuDataNormalization").getVal()
	corr = work.var(corr_name).getVal()
	orig_lumi = 2488.245
	scale = real_norm*corr*lumi/orig_lumi
	nevents = integral.getVal()*scale
	f.Close()
	return  nevents

def getDYNEvents(MWR, channel, lumi):
	if "ee" in channel or "EE" in channel:
		filename = "data/DYElectronFits.root"
		workname = "DYElectronFits"
		dataname = "DY_EEDataSet_120to200"
	else:
		filename = "data/DYMuonFits.root"
		workname = "DYMuonFits"
		dataname = "DY_MuMuDataSet_120to200"

	f = r.TFile.Open(filename)
	work = f.Get(workname)
	pdf = work.pdf("rescaledExpPdf")
	mass = work.var("fourObjectMass")
	data = work.data(dataname)

	mass.setRange(str(MWR), mass_cut[MWR][0], mass_cut[MWR][1])
	argset = r.RooArgSet(mass)
	integral = pdf.createIntegral(argset, r.RooFit.NormSet(argset), r.RooFit.Range(str(MWR)))

	norm = data.sumEntries()
	orig_lumi = 2488.245
	scale = norm*lumi/orig_lumi
	nevents = integral.getVal()*scale
	f.Close()
	return nevents

def getNEvents(MWR, channel, process, lumi):
	if process == "signal":
		h = getSignalMassHisto(MWR, channel)
		low,hi = mass_cut[int(MWR)]
		lowbin = h.FindBin(low)
		hibin = h.FindBin(hi) - 1
		nevents = h.Integral(lowbin, hibin)*lumi
		return nevents
	elif process == "TTBar":
		return getTTBarNEvents(MWR, channel, lumi)
	elif process == "DY":
		return getDYNEvents(MWR, channel, lumi)
	else:
		return None
##
# @brief calls and parses command for `combine'
#
# @param command a list to pass to subprocess.check_output
#
# @return (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)
def runCombine(command):
	try:
		output = subprocess.check_output(command)
		p = re.compile(r'median expected limit: r < ([0-9.]*)')
		median = p.search(output).group(1)
		p = re.compile(r'mean   expected limit: r < ([0-9.]*) \+/- ([0-9.]*)')
		mean,meanError = p.search(output).groups()
		p = re.compile(r'68% expected band : ([0-9.]*) < r < ([0-9.]*)')
		onesig_minus,onesig_plus = p.search(output).groups()
		p = re.compile(r'95% expected band : ([0-9.]*) < r < ([0-9.]*)')
		twosig_minus,twosig_plus = p.search(output).groups()
		return median, (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)
	except:
		errfile = open("_".join(command) + ".err","w")
		errfile.write(output)
		return None

mass_cut = {
		800 :( 700,  900),
		1200:(1000, 1400),
		1400:(1200, 1600),
		1600:(1350, 1850),
		2000:(1700, 2300),
		2400:(2050, 2750),
		2600:(2200, 3000),
		2800:(2400, 3200),
		3000:(2550, 3450),
		3200:(2700, 3700),
		3600:(3000, 4150),
		3800:(3000, 4350),
		4000:(3000, 4600),
		4200:(3000, 4850),
		4400:(3000, 5050),
		4600:(3000, 5300),
		4800:(3000, 5500),
		5000:(3000, 5750),
		5200:(3000, 6000),
		5600:(3000, 6450),
		5800:(3000, 6650),
		6000:(3000, 6900),
		}
