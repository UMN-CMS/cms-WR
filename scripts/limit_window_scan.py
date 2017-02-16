import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np

import sys


#MWRs = [ 800, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 2600, 2800, 3000, 3200, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5600, 5800, 6000]
#MWRs = [ 2000]
#chnl = "ee"
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
unweightedTreeName = "central_value_tree"
dy = combineTools.getMassHistoTwo(True, treename, ("DY_"+identifier, "analysisCppOutputRootFiles/selected_tree_DYMadInclAndHT_signal_" + channeltag + "_withMllWeight.root"), binsize=binsize)
dyUnwgt = combineTools.getMassHistoTwo(False, treename, ("UnweightedDY_"+identifier, "analysisCppOutputRootFiles/selected_tree_DYMadInclAndHT_signal_" + channeltag + "_withMllWeight.root"), binsize=binsize)
top = combineTools.getMassHistoTwo(True, treename, ("TT_"+identifier, "analysisCppOutputRootFiles/selected_tree_data_flavoursidebandEMuEE.root"), binsize=binsize)
topUnwgt = combineTools.getMassHistoTwo(False, treename, ("UnweightedTT_"+identifier, "analysisCppOutputRootFiles/selected_tree_data_flavoursidebandEMuEE.root"), binsize=binsize)
scale=1.0   #reduced below 1 for low MWR signal events only

#for systematics
nuisance_params = []
nuisance_params.append(("lumi",        "lnN"))
nuisance_params.append(("TTBar_SF",       "lnN"))
nuisance_params.append(("DY_SF",       "lnN"))
nuisance_params.append(("WR_"+identifier+"jj"+"_unc",  "gmN"))
nuisance_params.append(("TTBar_unc",      "gmN"))
nuisance_params.append(("DY_unc",   "gmN"))
dySFUnc = 1.4028
ttSFUnc = 1.0255
if "mu" in chnl:
	dySFUnc = 1.4034
	ttSFUnc = 1.0402

for MWR in MWRs:
	#binsize =MWR*.02
	signal = combineTools.getMassHistoTwo(True, treename, ("WR_"+identifier+"jj_" + str(MWR), "analysisCppOutputRootFiles/selected_tree_WRto%sJJ_%d_%d_signal_%s.root" % (wrfinalstate, MWR, MWR/2, channeltag) ), binsize=binsize)   #store a WR mass histo with weighted evts in signal
	signalUnwgt = combineTools.getMassHistoTwo(False, treename, ("UnweightedWR_"+identifier+"jj_" + str(MWR), "analysisCppOutputRootFiles/selected_tree_WRto%sJJ_%d_%d_signal_%s.root" % (wrfinalstate, MWR, MWR/2, channeltag) ), binsize=binsize)   #save a WR mass histo with unweighted evts in case a mass window is tested which contains zero weighted WR evts
	
	#np.linspace returns a numpy array of evenly spaced numbers which span the interval specified by the first two arguments
	#the third arg specifies the number of elements in the array
	low_mass_bins = np.linspace(.4,1, 11)*MWR
	hi_mass_bins = np.linspace(1,1.4, 11)*MWR
	for low_m, high_m in product(low_mass_bins, hi_mass_bins):
		if low_m > MWR: continue
		if high_m < MWR: continue
		if (high_m - low_m) < 200: continue
		if (high_m - low_m) > (MWR*0.5): continue
		if (high_m - low_m) > 2500: continue
		if low_m < 600: continue
		low_bin = signal.FindBin(low_m)
		hi_bin = signal.FindBin(high_m) - 1
		
		#events in signal and bkgnd histos are already weighted
		if MWR < 1300: scale = 0.1
		signal_tuple = ("WR_"+identifier+"jj", (signal.Integral(low_bin, hi_bin))*scale)
		bg_tuples = []
		if top.Integral(low_bin, hi_bin) > 0: bg_tuples.append(("TTBar", top.Integral(low_bin, hi_bin) ))
		else: bg_tuples.append(("TTBar", 0.0001 ))
		bg_tuples.append(("DY", dy.Integral(low_bin, hi_bin) ))

		#now sum the total number of bkgnd events
		nBG = sum( [ b for a,b in bg_tuples ] )

		#now calculate the stat and syst uncertainties and write them to a Systematics class object
		systematics = combineTools.Systematics(["WR_"+identifier+"jj", "TTBar", "DY"], nuisance_params)   #make systematics class objects for signal and bkgnds
		
		#add lumi uncertainty to signal and bkgnds, and PDF uncertainty in quadrature with lumi uncertainty for signal
		systematics.add("WR_"+identifier+"jj", "lumi", 1.202)
		systematics.add("TTBar", "lumi", 1.027)
		systematics.add("DY", "lumi", 1.027)
		
		#add other uncertainties affecting bkgnd normalization
		systematics.add("TTBar", "TTBar_SF", ttSFUnc)
		systematics.add("DY", "DY_SF", dySFUnc)

		#calculate stat uncertainties and compute N and alpha for signal and bkgnds
		meanTT = top.Integral(low_bin, hi_bin)
		varTT = meanTT
		NTT = 1
		alphaTT = 1
		avgWeightWR = (signal.Integral(low_bin, hi_bin)/signalUnwgt.Integral(low_bin, hi_bin))  #needed to compute variance
		meanWR = signal.Integral(low_bin, hi_bin)*scale
		varWR = meanWR*avgWeightWR
		NWR = 1
		alphaWR = 1
		avgWeightDY = (dy.Integral(low_bin, hi_bin)/dyUnwgt.Integral(low_bin, hi_bin))  #needed to compute variance
		meanDY = dy.Integral(low_bin, hi_bin)
		varDY = meanDY*avgWeightDY
		NDY = 1
		alphaDY = 1
		if meanTT > 0:
			NTT = meanTT**2/varTT
			alphaTT = varTT/meanTT
		else: #meanTT equals 0
			NTT = 0
			alphaTT = top.Integral(600, 6500)/(1 + topUnwgt.Integral(600, 6500))

		if meanDY > 0:
			NDY = meanDY**2/varDY
			alphaDY = varDY/meanDY
		else: #meanDY equals 0
			NDY = 0
			alphaDY = dy.Integral(600, 6500)/(1 + dyUnwgt.Integral(600, 6500))

		if meanWR > 0:
			NWR = meanWR**2/varWR
			alphaWR = varWR/meanWR
		else: #meanWR equals 0
			NWR = 0
			alphaWR = signal.Integral(600, 6500)/(1 + signalUnwgt.Integral(600, 6500))


		#add uncertainties (purely statistical here) which affect signal and bkgnd shapes
		systematics.add("WR_"+identifier+"jj", "WR_"+identifier+"jj"+"_unc", (NWR,alphaWR))
		systematics.add("TTBar", "TTBar_unc", (NTT,alphaTT))
		systematics.add("DY", "DY_unc", (NDY,alphaDY))


		#give the datacard file a unique name which identifies the WR mass point, lepton channel and lljj mass range
		datacard = "WR"+identifier +"jj_MASS%d_LOW%d_HIGH%d" % (MWR, low_m, high_m)
		datacard_file = "data/" + datacard + ".txt"
		
		#now make the datacard and store it in the txt file specified by datacard_file
		#last arg decides whether or not systematics are added to the datacard
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, identifier+"jj", nBG, signal_tuple, bg_tuples, systematics)
		method = "BayesianToyMC"
		ntoys = "150"   #default is 150
		##in command defined below: "-S0" disables evaluation of systematics, "-M %s" sets the combine run mode to the string 
		##stored in method, "-n %s" sets the job name equal to the string datacard
		command = "combine -M %s -H ProfileLikelihood -S1 %s -n %s --toys %s" %(method, datacard_file, datacard, ntoys)
		ret = combineTools.runCombine(command,datacard)
		mean, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
		
		print MWR, MWR/2, low_m, high_m, sig, ("%s "*len(bgs)) % bgs , twosig_minus, onesig_minus, mean, onesig_plus, twosig_plus
		subprocess.check_output(["rm", "higgsCombine" + datacard + "." + method + ".mH120.123456.root"])
