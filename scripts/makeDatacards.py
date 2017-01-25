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
nuisance_params.append(("signal_unc",  "gmN"))
nuisance_params.append(("TT_unc",      "gmN"))
nuisance_params.append(("DYAMC_unc",   "gmN"))
unscale_by_xs = False	#set to true to allow 800 GeV and 1.0 TeV MWR limit jobs to finish successfully
for channel in ["ee", "mumu"]:
	sig_name = "WR_" + channel + "jj"
	MWR = []
	signal = []
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

			MWR.append(mass)
			signal.append(signalNevents)
			bg.append([TTBar, DY])

			if args.nosyst: systematics = None
			systematics_list.append(systematics)
		except (IOError,KeyError) as e:
			print mass, "not found"

	bg_names = ["TTBar", "DY"]

	for i in range(len(MWR)):
		print channel, MWR[i], signal[i]/sum(bg[i])
		signal_tuple = (sig_name, signal[i])
		bg_tuples = zip(bg_names, bg[i])
		nBG = sum(bg[i])

		datacard = "WR%sjj_MASS%04d" % (channel, MWR[i])
		datacard_file = args.outdir + "/" + datacard + ".txt"
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
				signal_tuple, bg_tuples, systematics=systematics_list[i])
