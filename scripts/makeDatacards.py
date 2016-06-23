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

minitrees = combineTools.miniTreeInterface(
			base=args.basedir,
			tag =args.tag,
			makeplots=args.drawplots
			)

nuisance_params = []
nuisance_params.append( ("lumi",         "lnN", ))
nuisance_params.append( ("TT_ratio",     "lnN", ))
nuisance_params.append( ("DY_SF",        "lnN", ))
nuisance_params.append( ("signal_syst",  "lnN", ))
nuisance_params.append( ("signal_stat",  "gmN", ))
nuisance_params.append( ("TT_syst",      "lnN", ))
nuisance_params.append( ("TT_stat",      "gmN", ))
nuisance_params.append( ("DYAMC_syst",   "lnN", ))
nuisance_params.append( ("DYAMC_stat",   "gmN", ))
unscale_by_xs = True
for channel in ["ee", "mumu"]:
	sig_name = "WR_" + channel + "jj"
	MWR = []
	signal = []
	bg = []
	systematics_list = []
	for mass in sorted(combineTools.mass_cut):
		try:
			systematics = combineTools.Systematics(["signal", "TT", "DYAMC"], nuisance_params)
			signalNevents, sig_syst, sig_stat = minitrees.getNEvents(mass, channel, "signal", systematics)
			if unscale_by_xs:
				uns = .001/xs.WR_jj[channel][mass]
				signalNevents *= uns

			TTBar, TTBar_syst, TTBar_stat = minitrees.getNEvents(mass, channel, "TT", systematics)
			DY, DY_syst, DY_stat = minitrees.getNEvents(mass, channel, "DYAMC", systematics)

			MWR.append(mass)
			signal.append(signalNevents)
			bg.append([TTBar, DY])

			if args.nosyst: systematics = None
			systematics_list.append(systematics)
		except IOError:
			print mass, "not found"

	bg_names = ["TTBar", "DY"]

	for i in range(len(MWR)):
		print MWR[i], signal[i], bg[i]
		signal_tuple = (sig_name, signal[i])
		bg_tuples = zip(bg_names, bg[i])
		nBG = sum(bg[i])

		datacard = "WR%sjj_MASS%04d" % (channel, MWR[i])
		datacard_file = args.outdir + "/" + datacard + ".txt"
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
				signal_tuple, bg_tuples, systematics=systematics_list[i])
