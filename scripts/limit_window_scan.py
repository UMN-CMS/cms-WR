import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np

import sys


#MWRs = [ 800, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 2600, 2800, 3000, 3200, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5600, 5800, 6000]
#MWRs = [ 6000]
#chnl = "mumu"
MWRs = [ WRMASS]
chnl = "CHNL"


from itertools import combinations, product
binsize = 50

if len(sys.argv) > 1:
	MWRs = [MWRs[int(i)] for i in sys.argv[1:]]

#make WR mass histos with weighted evts for bkgnds
identifier = "ee"
channeltag = "eeEE"
wrfinalstate = "EE"
if "mu" in chnl:
	identifier = "mumu"
	channeltag = "mumuMuMu"
	wrfinalstate = "MuMu"

treename = "Tree_Iter0"
dy = combineTools.getMassHistoTwo(treename, ("DY_"+identifier, "analysisCppOutputRootFiles/selected_tree_DYMadInclAndHT_signal_" + channeltag + "_withMllWeight.root"), binsize=binsize)
top = combineTools.getMassHistoTwo(treename, ("TT_"+identifier, "analysisCppOutputRootFiles/selected_tree_data_flavoursidebandEMuEE.root"), binsize=binsize)

for MWR in MWRs:
	#binsize =MWR*.02
	signal = combineTools.getMassHistoTwo(treename, ("WR_"+identifier+"jj_" + str(MWR), "analysisCppOutputRootFiles/selected_tree_WRto%sJJ_%d_%d_signal_%s.root" % (wrfinalstate, MWR, MWR/2, channeltag) ), binsize=binsize)   #store a WR mass histo with weighted evts in signal
	
	#np.linspace returns a numpy array of evenly spaced numbers which span the interval specified by the first two arguments
	#the third arg specifies the number of elements in the array
	low_mass_bins = np.linspace(.4,1, 21)*MWR
	hi_mass_bins = np.linspace(1,1.4, 21)*MWR
	for low_m, high_m in product(low_mass_bins, hi_mass_bins):
		if low_m > MWR: continue
		if high_m < MWR: continue
		if (high_m - low_m) < 200: continue
		if (high_m - low_m) > (MWR*0.5): continue
		if (high_m - low_m) > 2500: continue
		low_bin = signal.FindBin(low_m)
		hi_bin = signal.FindBin(high_m) - 1
		
		#events in signal and bkgnd histos are already weighted
		signal_tuple = ("WR_"+identifier+"jj", signal.Integral(low_bin, hi_bin))
		bg_tuples = []
		bg_tuples.append(("TTBar", top.Integral(low_bin, hi_bin) ))
		bg_tuples.append(("DY", dy.Integral(low_bin, hi_bin) ))

		#now sum the total number of bkgnd events
		nBG = sum( [ b for a,b in bg_tuples ] )

		#give the datacard file a unique name which identifies the WR mass point, lepton channel and lljj mass range
		datacard = "WR"+identifier +"jj_MASS%d_LOW%d_HIGH%d" % (MWR, low_m, high_m)
		datacard_file = "data/" + datacard + ".txt"
		
		#now make the datacard and store it in the txt file specified by datacard_file
		#last arg decides whether or not systematics are added to the datacard
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, identifier+"jj", nBG, signal_tuple, bg_tuples, False)
		method = "ProfileLikelihood"
		ntoys = "300"
		##in command defined below: "-S0" disables evaluation of systematics, "-M %s" sets the combine run mode to the string 
		##stored in method, "-n %s" sets the job name equal to the string datacard
		command = "combine -M %s -S0 %s -n %s --toys %s" %(method, datacard_file, datacard, ntoys)
		ret = combineTools.runCombine(command,datacard)
		mean, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
		
		print MWR, MWR/2, low_m, high_m, sig, ("%s "*len(bgs)) % bgs , twosig_minus, onesig_minus, mean, onesig_plus, twosig_plus
		subprocess.check_output(["rm", "higgsCombine" + datacard + "." + method + ".mH120.123456.root"])
