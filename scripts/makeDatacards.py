import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.cross_sections as xs 

import argparse

parser = argparse.ArgumentParser(description='Make datacards')
parser.add_argument('--no-syst', dest='nosyst', action='store_true',
		help='do not write systematics to datacards')
parser.add_argument('--draw-plots', dest='drawplots', action='store_true',
		help='draw plots')
parser.add_argument('-d', '--dir', dest='basedir',
		default="./",
		help='base dir for analysis tree files')
parser.add_argument('-t', '--tag', dest='tag',
		default="",
		help='tag name for analysis tree files')
parser.add_argument('-o', '--outdir', dest='outdir',
		default="datacards/",
		help='where to store the datacards')


args = parser.parse_args()

minitrees = combineTools.AnalysisResultsInterface(
			base=args.basedir,
			tag =args.tag,
			makeplots=args.drawplots
			)

nuisance_params = []
nuisance_params.append(("lumi",        "lnN"))
nuisance_params.append(("TT_SF",       "lnN"))
nuisance_params.append(("DYAMC_SF",       "lnN"))

#omit the signal_unc nuisance parameter when making datacards to calculate the observed limit
nuisance_params.append(("signal_unc",  "gmN"))

nuisance_params.append(("TT_unc",      "gmN"))
nuisance_params.append(("DYAMC_unc",   "gmN"))
#nuisance_params.append(("OTHER_unc",   "gmN"))
unscale_by_xs = True	#set to true to allow 800 GeV and 1.0 TeV MWR limit jobs to finish successfully
for channel in ["ee", "mumu"]:
	sig_name = "WR_" + channel + "jj"
	MWR = []
	signal = []
	obs = []
	bg = []
	systematics_list = []
	for mass in sorted(combineTools.mass_cut[channel]):
		try:
			systematics = combineTools.Systematics(["signal", "TT", "DYAMC"], nuisance_params)
			if unscale_by_xs:
				scale =  .001/xs.WR_jj[channel][mass]
			else:
				scale = 1.0
			signalNevents = minitrees.getNEvents(mass, channel, "signal", systematics, scale = scale)

			TTBar = minitrees.getNEvents(mass, channel, "TT", systematics)
			DY = minitrees.getNEvents(mass, channel, "DYAMC", systematics)
			
			#this data value is the number of events observed in data in a specific mass window
			#only needed when making datacards to calculate the observed limit
			data = minitrees.getNEvents(mass, channel, "data", systematics)
			
			#saved for reference, in case another process needs to be added to the datacard
			#Other = minitrees.getNEvents(mass, channel, "OTHER", systematics)


			MWR.append(mass)
			
			#add the number of events observed in data to the list named signal so that the datacard
			#lists the number of obs events in data in both the observation field and the signal events field, next to the expected TT and DY events
			signal.append(signalNevents)
			#signal.append(data)

			obs.append(data)
			
			bg.append([TTBar, DY])

			if args.nosyst: systematics = None
			systematics_list.append(systematics)
		except (IOError,KeyError) as e:
			print mass, "not found"

	bg_names = ["TTBar", "DY"]   #names shown in the expected evt columns of each datacard

	for i in range(len(MWR)):
		print channel, MWR[i], signal[i]/sum(bg[i])
		signal_tuple = (sig_name, signal[i])
		bg_tuples = zip(bg_names, bg[i])
		nBG = sum(bg[i])
		nObs = obs[i]

		datacard = "WR%sjj_MASS%04d" % (channel, MWR[i])
		datacard_file = args.outdir + "/" + datacard + ".txt"
		
		#the third argument in makeDataCardSingleBin sets the number of observed events in the datacard
		#use nObs for observed limits, and nBG for expected limits
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nObs,
				signal_tuple, bg_tuples, systematics=systematics_list[i])
