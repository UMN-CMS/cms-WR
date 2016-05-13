#!/usr/bin/env python

import ROOT as r
import re
import numpy as np
import subprocess
from random import gauss
import ExoAnalysis.cmsWR.backgroundFits as bgfits
import ExoAnalysis.cmsWR.cross_sections as xs
import math
import datetime
from os import environ

datafolder = environ['LOCALRT'] + "/src/ExoAnalysis/cmsWR/data/"
configfolder = environ['LOCALRT'] + "/src/ExoAnalysis/cmsWR/configs/"
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
		out.write("observation %d\n" % round(nObs))
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

class miniTreeInterface:
	chnlName = {"ee":"EE","mumu":"MuMu"}
	def __init__(self,
			base="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/forPeterRootFiles/",
			tag = "noMllAndZptWeights", ):

		if tag: tag = "_" + tag
		self.filefmt_dict = {"base":base, "tag":tag}
		self.filefmt = "{base}/selected_tree_{mode}_signal_{channel}{chnlName}{tag}.root"

		emufile = r.TFile.Open("flavor_fits.root")
		if emufile:
			f_EE = emufile.Get("f_EE")
			f_MuMu = emufile.Get("f_MuMu")
			self.tt_emu_ratio = {
					"ee": f_EE.GetParameter(1),
					"mumu": f_MuMu.GetParameter(1),
					}
			self.tt_emu_error = {
					"ee": f_EE.GetParError(1),
					"mumu": f_MuMu.GetParError(1),
					}

	def getNEvents(self, MWR, channel, process):
		MWR = int(MWR)
		if process == "signal":
			mode = "WRto" + self.chnlName[channel] + "JJ_" + str(MWR) + "_" + str(MWR/2)
		elif "TT" in process or "DY" in process:
			mode = process
		else:
			return None

		self.filefmt_dict["channel"] = channel
		self.filefmt_dict["chnlName"] = self.chnlName[channel]
		self.filefmt_dict["mode"] = mode

		filename = self.filefmt.format(**self.filefmt_dict)
		f = r.TFile.Open(filename)
		if not f: raise IOError


		masses = r.std.vector(int)()
		f.GetObject("signal_mass",masses)
		masses = [m for m in masses]
		mass_i = masses.index(int(MWR))
		tree = f.Get("tf1")
		if "TT" in process or "DY" in process:
			num = np.array([event.FitIntegralInRange[mass_i]*event.Normalization for event in tree])
		else:
			num = np.array([event.NEventsInRange[mass_i] for event in tree])

		mean = np.mean(num)
		std = np.std(num)
		if "TT" in process:
			scale = self.tt_emu_ratio[channel]
			mean *= scale
			std *= scale

		syst = 1 + std/mean
		return mean, syst

##
# @brief calls and parses command for `combine'
#
# @param command a list to pass to subprocess.check_output
#
# @return (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus)
def runCombine(command):
	name  = "combine" + datetime.datetime.now().time().isoformat()
	out_file = open(name + ".out", "w+")
	err_file = open(name + ".err", "w")
	command = "unbuffer " + command
	command = command.split(' ')
	jobname = command[command.index('-n') + 1]
	method = command[command.index('-M') + 1]
	print method, "HybridNew" in method, jobname
	try:
		seed = command[command.index('-s') + 1]
	except ValueError:
		seed = "123456"
	try:
		mass = command[command.index('-m') + 1]
	except ValueError:
		mass = "120"
		
	try:
		if "HybridNew" in method:
			rs = []
			print "do different quantiles"
			for q in [.025, .16, 0.5, .84, .975]:
				#print q
				run_command = command + ["--expectedFromGrid=%f" % q]
				#print " ".join(run_command)
				subprocess.call(run_command, stdout=out_file, stderr=err_file)
				out_file.seek(0)
				output = out_file.read()
				p = re.compile(r'Limit: r < ([0-9.]*)')
 				matches  = p.findall(output)
				rs.append(matches[-1])
			twosig_minus, onesig_minus, median, onesig_plus, twosig_plus = tuple(rs)
		else:
			#print " ".join(command)
			subprocess.call(command, stdout=out_file, stderr=err_file)
			out_file.seek(0)
			output = out_file.read()
			if not "--toys" in command and not "-t" in command:
				p = re.compile(r'Limit: r < ([0-9.]*)')
 				matches  = p.findall(output)
 				if not matches: raise RuntimeError
				return matches[-1]
			
			
			outfile = r.TFile.Open("higgsCombine" + jobname + "." + method + ".mH" + mass + "." + seed + ".root")
			limitTree = outfile.Get("limit")
			limits = np.zeros(limitTree.GetEntries())
			for i in range(limitTree.GetEntries()):
				limitTree.GetEntry(i)
				limits[i] = limitTree.limit 
			#limitTree.Draw("limit>>tmphist")
			#h = r.gDirectory.Get("tmphist")
			q = [2.5, 16, 50, 84, 97.5]
			twosig_minus, onesig_minus, median, onesig_plus, twosig_plus = np.percentile(limits, q)

		if not all([median, onesig_minus, onesig_plus, twosig_minus, twosig_plus]):
			print "combine parse failed"
			print median, onesig_minus, onesig_plus, twosig_minus, twosig_plus
			raise TypeError
		return median, (onesig_minus, onesig_plus), (twosig_minus, twosig_plus)
	except Exception as e:
		raise e
		return None

mass_cut =  {mass:(low,hi) for mass,low,hi in [ map(int,s.split()) for s in open(configfolder + "mass_cuts.txt",'r').read().split('\n') if s and s[0] != "#"]}

import sys
if __name__ == '__main__':
	ID = sys.argv[1]
	command = " ".join(sys.argv[2:])
	print command
	result = runCombine(command)
	print ("COMBINE", ID, result)