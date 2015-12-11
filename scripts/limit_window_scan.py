import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.backgroundFits as bgFits
import numpy as np

import sys



lumi = 1900
MWRs = [ 800, 1200, 1400, 1600, 2000, 2400, 2600, 2800, 3000, 3200, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5600, 5800, 6000]
from itertools import combinations, product
binsize = 100

if len(sys.argv) > 1:
	MWRs = [MWRs[int(i)] for i in sys.argv[1:]]

for MWR in MWRs:
	binsize =MWR*.02
	dirname =  "unmatchedSignalRecoAnalyzerFive"
	treename = "signalRecoObjectsWithAllCuts"
	signal = combineTools.getMassHisto(dirname,  treename,  ("WR_eejj_"+ str(MWR), "data/all_analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-%d_Nu_M-%d.root" % (MWR, MWR/2) ), binsize=binsize)
	
	mid_bin = signal.FindBin(MWR)
	#while signal.Integral(mid_bin - bin_range, mid_bin + bin_range)/ signal.Integral() < .98:
		#bin_range +=1
	#search = bin_range*binsize
	low_mass_bins = np.linspace(.4,.8, 21)*MWR
	hi_mass_bins = np.linspace(1,1.4, 21)*MWR
	#mass_bins = range(int(MWR - search), int(MWR + (search + binsize)), int(binsize))
	#for low_m, high_m in combinations(mass_bins, 2):
	for low_m, high_m in product(low_mass_bins, hi_mass_bins):
		if low_m > MWR: continue
		if high_m < MWR: continue
		low_bin = signal.FindBin(low_m)
		hi_bin = signal.FindBin(high_m) - 1
		signal_tuple = ("WR_eejj", signal.Integral(low_bin, hi_bin) * (0.01/50000.) * lumi)
		bg_tuples = []
		#bg_tuples.append( ("TTbar", bgFits.TTbarFit.Integral(low_m, high_m)/bgFits.binWidth* (lumi /2000.) ))
		#bg_tuples.append( ("DY",    bgFits.DYFit.Integral(low_m, high_m)   /bgFits.binWidth* (lumi /2000.) ))
		bg_tuples.append( ("AllBG",  bgFits.All_ee.Integral(low_m, high_m) * (lumi /2000.) ))

		datacard = "WReejj_MASS%d_LOW%d_HIGH%d" % (MWR, low_m, high_m)
		datacard_file = "data/" + datacard + ".txt"
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, "eejj", 0, signal_tuple, bg_tuples)
		method = "ProfileLikelihood"
		command = ["combine","-M",method,"-S","0",datacard_file,"-n", datacard,"-t","1000"]
		ret = combineTools.runCombine(command)
		(mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret

		print MWR, MWR/2, low_m, high_m, sig, ("%s "*len(bgs)) % bgs , twosig_minus, onesig_minus, mean, onesig_plus, twosig_plus
		subprocess.check_output(["rm", "higgsCombine" + datacard + "." + method + ".mH120.123456.root"])
