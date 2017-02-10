import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np

import sys


#MWRs = [ 800, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 2600, 2800, 3000, 3200, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5600, 5800, 6000]
MWRs = [ 2000]
#MWRs = [ WRMASS]



from itertools import combinations, product
binsize = 50

if len(sys.argv) > 1:
	MWRs = [MWRs[int(i)] for i in sys.argv[1:]]

#make WR mass histos with weighted evts for bkgnds
identifier = "ee"
treename = "Tree_Iter0"
dy = combineTools.getMassHistoTwo(treename, ("DY_"+identifier, "analysisCppOutputRootFiles/selected_tree_DYMadInclAndHT_signal_eeEE_withMllWeight.root"), binsize=binsize)
top = combineTools.getMassHistoTwo(treename, ("TT_"+identifier, "analysisCppOutputRootFiles/selected_tree_data_flavoursidebandEMuEE.root"), binsize=binsize)

for MWR in MWRs:
	#binsize =MWR*.02
	signal = combineTools.getMassHistoTwo(treename, ("WR_eejj_"+ str(MWR), "analysisCppOutputRootFiles/selected_tree_WRtoEEJJ_%d_%d_signal_eeEE.root" % (MWR, MWR/2) ), binsize=binsize)   #store a WR mass histo with weighted evts in signal
	
	#in the TH1F stored in signal (defined immediately above), get the bin number which corresponds to MLLJJ equal to MWR
	#this should be the bin with the greatest entries
	mid_bin = signal.FindBin(MWR)
	#while signal.Integral(mid_bin - bin_range, mid_bin + bin_range)/ signal.Integral() < .98:
		#bin_range +=1
	#search = bin_range*binsize

	#np.linspace returns a numpy array of evenly spaced numbers which span the interval specified by the first two arguments
	#the third arg specifies the number of elements in the array
	low_mass_bins = np.linspace(.4,1, 21)*MWR
	hi_mass_bins = np.linspace(1,1.4, 21)*MWR
	#mass_bins = range(int(MWR - search), int(MWR + (search + binsize)), int(binsize))
	#for low_m, high_m in combinations(mass_bins, 2):
	for low_m, high_m in product(low_mass_bins, hi_mass_bins):
		if low_m > MWR: continue
		if high_m < MWR: continue
		if (high_m - low_m) < 200: continue
		if (high_m - low_m) > (MWR*0.5): continue
		if (high_m - low_m) > 2000: continue
		low_bin = signal.FindBin(low_m)
		hi_bin = signal.FindBin(high_m) - 1
		
		#events in signal and bkgnd histos are already weighted
		signal_tuple = ("WR_eejj", signal.Integral(low_bin, hi_bin))
		bg_tuples = []
		bg_tuples.append(("dy", dy.Integral(low_bin, hi_bin) ))
		bg_tuples.append(("top", top.Integral(low_bin, hi_bin) ))

		#now sum the total number of bkgnd events
		nBG = sum( [ b for a,b in bg_tuples ] )

		#give the datacard file a unique name which identifies the WR mass point, lepton channel and lljj mass range
		datacard = "WReejj_MASS%d_LOW%d_HIGH%d" % (MWR, low_m, high_m)
		datacard_file = "data/" + datacard + ".txt"
		
		#now make the datacard and store it in the txt file specified by datacard_file
		#last arg decides whether or not systematics are added to the datacard
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, "eejj", nBG, signal_tuple, bg_tuples, False)
		method = "ProfileLikelihood"
		ntoys = "200"
		#in command defined below: "-S0" disables evaluation of systematics, "-M %s" sets the combine run mode to the string 
		#stored in method, "-n %s" sets the job name equal to the string datacard
		command = "combine -M %s -S0 %s -n %s --toys %s" %(method, datacard_file, datacard, ntoys)
		ret = combineTools.runCombine(command,datacard)
		mean, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
		
		print MWR, MWR/2, low_m, high_m, sig, ("%s "*len(bgs)) % bgs , twosig_minus, onesig_minus, mean, onesig_plus, twosig_plus
		subprocess.check_output(["rm", "higgsCombine" + datacard + "." + method + ".mH120.123456.root"])
