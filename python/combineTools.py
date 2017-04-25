#!/usr/bin/env python

import ROOT as r
import re
import numpy as np
import subprocess
from collections import defaultdict
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
def makeDataCardSingleBin(outfile, bin_name, nObs, signal_tuple, background_tuples, systematics):

	nBGs = len(background_tuples)
	ns = "  ".join([str(i) for i in range(nBGs+1)])
	sig_name = signal_tuple[0]
	signal_rate = "%.4g" % (signal_tuple[1])

	bg_names = ""
	bg_rates = ""
	name_lookup = {sig_name:0}
	for i, (bg_name, bg_rate) in enumerate(background_tuples):
		name_lookup[bg_name] = i+1
		bg_names += bg_name + "  "
		bg_rates   += "%.4g" % bg_rate + "  "
	names = sig_name + "  " + bg_names
	rates = signal_rate + '  ' + bg_rates
	with open(outfile, "w") as out:
		out.write("imax 1  number of channels\n")
		out.write("jmax %d  number of backgrounds\n" % nBGs)
		if systematics:
			out.write("kmax %d  number of nuisance parameters\n" % systematics.n_nuisance)
		out.write("bin " + bin_name + "\n")
		out.write("observation %d\n" % round(nObs))
		out.write("------------\n")
		out.write("bin" + ("    " + bin_name)* (nBGs + 1) + "\n")
		out.write("process  " + names + "\n")
		out.write("process  " + ns + "\n")
		out.write("rate  " + rates + "\n")
		out.write("------------\n")
		if systematics:
			out.write(str(systematics))
	return signal_rate, tuple(bg_rates.split())

#getOtherMassHisto allows an individual root file containing a TTree to be opened and read
#this currently draws unweighted evts as it is used to determine the total efficiency of
#online + offline + mass window cuts in the signal region, assuming 50k generated evts in the dataset
def getOtherMassHisto(cutstring, treename, process, binsize=100):
	histname, filename = process
	infile = r.TFile.Open(filename)
	tree = infile.Get(treename)
	r.gROOT.cd()
	nbins = 10000/int(binsize)
	h = r.TH1F(histname, histname, nbins, 0, nbins * binsize)
	if "TT" not in histname:
		tree.Draw("WR_mass>>" + histname, "1*("+cutstring+")", "goff")

	#if the histo corresponds to TT, then weight the events with the appropriate SF
	if "TT_ee" in histname:
		tree.Draw("WR_mass>>" + histname, "0.4194*("+cutstring+")", "goff")
	if "TT_mumu" in histname:
		tree.Draw("WR_mass>>" + histname, "0.6563*("+cutstring+")", "goff")

	return h


##end getOtherMassHisto



class AnalysisResultsInterface:
	chnlName = {"ee":"EE", "mumu":"MuMu", "emu":"EMu"}
	def __init__(self,
			base="./",
			tag = "",
			SFfilename = configfolder + "/MCScaleFactors.txt",
			makeplots = True, #default is False
			):

		self.makeplots = makeplots
		if tag: tag = "_" + tag
		self.filefmt_dict = {"base":base, "tag":tag}
		self.filefmt = "{base}/selected_tree_{mode}_{minitreename}{chnlName}{tag}.root"

		self.SF = defaultdict(lambda : defaultdict(dict))
		with open(SFfilename) as f:
			for line in f:
				if line[0] == "#": continue
				mode, ch, sf, stat, syst  = line.split()
				ch = ch.lower()
				self.SF[mode][ch]["SF"] = float(sf)
				self.SF[mode][ch]["unc"] = math.sqrt(float(stat)**2 + float(syst)**2)


		self.masses = []
		self.results = {}

	def getUncertainty(self, mode, channel):
		return 1 + self.SF[mode][channel]["unc"]/self.SF[mode][channel]["SF"]

	def getNEvents(self, MWR, channel, process, systematics, scale = 1.0):
		""" returns mean, syst, stat """
		key = channel + "_" + process
		if "signal" == process:
			key += "_" + str(MWR)

		MWR = int(MWR)
		if key not in self.results:
			f = self.OpenFile(channel, process, MWR)
			if not self.masses: self.GetMasses(f)
			self.ProcessFile(key, f)

		mass_i = self.masses.index(MWR)

		syst_mean = self.results[key]["syst"]["mean"][mass_i]*scale
		tmp_syst  = self.results[key]["syst"]["std"] [mass_i]*scale
		tmp_stat  = self.results[key]["stat"]        [mass_i]*scale
		central_value      = self.results[key]["central"]["weighted"][mass_i]*scale
		central_unweighted = self.results[key]["central"]["unweighted"][mass_i]


		var = tmp_syst**2 + tmp_stat**2
		if syst_mean > 0:
			mean = syst_mean
			alpha = var/mean
			N = mean**2/var
			if N < 1:
				N = 0
				mean = math.sqrt(var)
				alpha = mean
		elif syst_mean < 0:
			mean = math.sqrt(var)
			alpha = mean
			N = 1
		else:
			N=0
			alpha = self.results[key]["syst"]["mean"][0]*scale/( self.results[key]["central"]["unweighted"][0] + 1)
			mean = .0001

		raw = [ syst_mean, mean, central_value, tmp_syst, tmp_stat, var, central_unweighted, ]
		#lump the WR pdf uncertainty in with the lumi uncertainty, but don't change the lumi uncertainty for DY or TT or OTHER
		if process in ["DYAMC", "TT"]:
			systematics.add(process, "lumi", 1.027)
			raw += [1.027]
		else:
			systematics.add(process, "lumi", 1.037)
			raw += [1.037]

		
		systematics.add(process, process + "_unc", (N,alpha))
		raw +=[N,alpha]
		
		if process in ["DYAMC", "TT"]:
			systematics.add(process,process + "_SF", self.getUncertainty(process, channel))
			raw += [self.getUncertainty(process, channel)]
		else:
			raw += [0.0]

		if "_800" == key[-4:]:
			key = key[:-3] + "0800"
		raw = ["%0.4f" % x for x in raw]
		raw = ["raw", key, "%04d" % MWR] + raw
		print " ".join(raw)

		return mean

	def getMeanStd(self, tree, fromFit=False):
		means = np.zeros(len(self.masses))
		stds = np.zeros(len(self.masses))
		if fromFit:
			draw_str = "FitIntegralInRange[%d]*Normalization"
		else:
			draw_str = "NEventsInRange[%d]"

		if self.makeplots: c = r.TCanvas("c", "c", 600, 600)
		for mass_i in range(len(self.masses)):
			ms = str(self.masses[mass_i])
			if "signal" in self.currentkey and self.currentkey.split("_")[2] != ms: continue

			tree.Draw(draw_str % mass_i, "", "goff")
			h = r.gDirectory.Get("htemp")
			means[mass_i] = h.GetMean()
			stds[mass_i] = h.GetStdDev()

			if self.makeplots:
				r.gStyle.SetOptStat(1001110)
				c.SetLogy()
				h.SetTitle(self.currentkey + tree.GetName() + " Mass " + ms)
				h.Draw()
				c.SaveAs("plots/" + self.currentkey  + tree.GetName() + "_mass" + ms + ".png")

		return means, stds

	def ProcessFile(self, key, f):
		if key in self.results: return

		self.currentkey = key
		tree = f.Get("syst_tree")
		if not tree or tree.GetEntries() == 0:
			tree = f.Get("central_value_tree")
		syst_means, syst_stds = self.getMeanStd(tree) 

		central_tree = f.Get("central_value_tree")
		central_tree.GetEntry(0)
		try:
			central_value = np.array(central_tree.NEventsInRange)
		except AttributeError:
			central_value = np.zeros(len(self.masses))
		try:
			central_error = np.array(central_tree.ErrorEventsInRange)
		except AttributeError:
			central_error = np.zeros(len(self.masses))
		try:
			central_unweighted = np.array(central_tree.UnweightedNEventsInRange)
		except AttributeError:
			central_unweighted = np.zeros(len(self.masses))

		if "TT" in key or "DYAMC" in key:
			if "TT" in key: mode = "TT"
			elif "DYAMC" in key: mode = "DYAMC"
			#elif "OTHER" in key: mode = "OTHER"
			if "ee" in key: channel = "ee"
			elif "mumu" in key: channel = "mumu"
			scale = self.SF[mode][channel]["SF"]
			syst_means *= scale
			syst_stds *= scale
			central_value *= scale
			central_error *= scale

		self.results[key] = {
				"syst": {
					"mean":syst_means.tolist(),
					"std":syst_stds.tolist()
					},
				"stat": central_error.tolist(),
				"central": {
					"weighted": central_value.tolist(),
					"unweighted": central_unweighted.tolist(),
					}
				}


	def OpenFile(self, channel, process, MWR):
		oldtag = self.filefmt_dict["tag"]
		if process == "signal":
			mode = "WRto" + self.chnlName[channel] + "JJ_" + str(MWR) + "_" + str(MWR/2)
			if MWR == 1800: 
				mode = "WRto" + self.chnlName[channel] + "JJ_1800_1400"
			minitreename = "signal_" + channel
		elif "TT" in process:
			self.filefmt_dict["tag"] = self.chnlName[channel]
			channel = "emu"
			mode = "data"
			minitreename = "flavoursideband"
		elif "DY" or "data" in process:
			mode = process
			minitreename = "signal_" + channel
		#elif "OTHER" in process:
		#	mode = "topW"
		#	minitreename = "signal_" + channel
		else:
			return None

		self.filefmt_dict["minitreename"] = minitreename
		self.filefmt_dict["chnlName"] = self.chnlName[channel]
		self.filefmt_dict["mode"] = mode

		filename = self.filefmt.format(**self.filefmt_dict)
		self.filefmt_dict["tag"] = oldtag
		print "Opening File ", filename
		f = r.TFile.Open(filename)
		if not f: raise IOError
		return f

	def GetMasses(self, f):
		masses = r.std.vector(int)()
		f.GetObject("signal_mass", masses)
		self.masses = [m for m in masses]
	
	def setTag(self, tag):
		self.filefmt_dict["tag"] = tag

	def getMassHisto(self, MWR, channel, process):
		r.gROOT.SetBatch(True)
		key = channel + "_" + process
		if "signal" == process:
			key += "_" + str(MWR)

		MWR = int(MWR)
		f = self.OpenFile(channel, process, MWR)
		if not self.masses: self.GetMasses(f)

		tree = f.Get("Tree_Iter0")

		mass_histo = r.TH1D("mass","mass",590,600, 6500)
		tree.Draw("WR_mass>>mass","","")

		mass_histo.SetDirectory(0)
		return mass_histo

		
class Systematics:
	#for formatting. if a string format as string otherwise use fmt_spec
	#class nuisance_value:
	#	def __init__(self, value):
	#		self.value = value
	#	def __format__(self, spec):
	#		if(isinstance(self.value, str)): return str(self.value)
	#		else: return format(self.value, spec)

	def __init__(self, channel_names, nuisance_params):
		self.channel_names = channel_names
		self.rows = nuisance_params
		self.n_nuisance = len(self.rows)
		self.values = {}
	def add(self, channel, syst, value):
		key = channel+syst
		self.values[key] = value
	def __str__(self):
		maxlen = max([len(name) for name, s_type in self.rows])
		floatsize = 10
		prefix = "{name:{maxlen}} {type} {N:6}"
		fmt = " {:{floatsize}.4e}"
		s = ""
		for name, s_type in self.rows:
			line = prefix 
			N = ""
			for channel in self.channel_names:
				key = channel + name
				try:
					N, v = self.values[key]
					N = int(round(N))
					line += fmt.format(v, floatsize=floatsize)
				except KeyError as e:
					line += " " * floatsize + "-"
				except TypeError:
					v = self.values[key]
					line += fmt.format(v, floatsize=floatsize )
			s += line.format(name=name, type=s_type, N=N, maxlen=maxlen) + "\n"
		return s
	def __len__(self):
		return len(self.values)

		




##
# @brief calls and parses command for `combine'
#
# @param command a list to pass to subprocess.check_output
#
# @return (mean, meanError), (onesig_minus, onesig_plus), (twosig_minus, twosig_plus)
def runCombine(command, ID):
	name  = "combine" + ID
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
			limits = []
			for toy in limitTree:
				if not np.isnan(toy.limit):
					limits.append(toy.limit)
			#limitTree.Draw("limit>>tmphist")
			#h = r.gDirectory.Get("tmphist")
			q = [2.5, 16, 50, 84, 97.5]
			twosig_minus, onesig_minus, median, onesig_plus, twosig_plus = np.percentile(np.array(limits), q)

		if not all([median, onesig_minus, onesig_plus, twosig_minus, twosig_plus]):
			print "combine parse failed"
			print median, onesig_minus, onesig_plus, twosig_minus, twosig_plus
			raise TypeError
		return median, (onesig_minus, onesig_plus), (twosig_minus, twosig_plus)
	except Exception as e:
		raise e
		return None

mass_cut = {"ee":{},
		      "mumu":{}
		     }
with open(configfolder + "mass_cuts.txt") as mc_file:
	for line in mc_file:
		if line[0] == "#": continue
		ch, mass, low, hi = line.split()
		mass_cut[ch.lower()][int(mass)] = (int(low), int(hi))

import sys
if __name__ == '__main__':
	ID = sys.argv[1]
	command = " ".join(sys.argv[2:])
	print ID, command
	result = runCombine(command, ID)
	print ("COMBINE", ID, result)
