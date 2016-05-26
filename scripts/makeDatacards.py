import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.cross_sections as xs 

import argparse

parser = argparse.ArgumentParser(description='Make datacards')
parser.add_argument('--no-syst', dest='nosyst', action='store_true',
		help='do not write systematics to datacards')
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
			makeplots=True
			)

unscale_by_xs = True
for channel in ["ee","mumu"]:
	sig_name = "WR_" + channel + "jj"
	MWR = []
	signal = []
	bg = []
	systematics_list  = []
	for mass in sorted(combineTools.mass_cut):
		try:
			systematics = []
			signalNevents, sig_syst, sig_stat = minitrees.getNEvents(mass, channel, "signal")
			if unscale_by_xs:
				uns = .001/xs.WR_jj[channel][mass]
				signalNevents *= uns

			TTBar, TTBar_syst, TTBar_stat = minitrees.getNEvents(mass, channel, "TT")
			minitrees.setTag("_withMllWeight")
			DY, DY_syst, DY_stat = minitrees.getNEvents(mass, channel, "DYAMC")
			minitrees.setTag("")

			MWR.append(mass)
			signal.append(signalNevents)
			bg.append([TTBar, DY])

			#TODO: add Lumi uncertainty and others
			#systematics.append( ("", "", [ (sig_name, sig_syst )]))
			print "nosyst", args.nosyst
			if not args.nosyst:
				systematics.append( ("Signal_syst", "lnN", [ (sig_name, sig_syst )]))
				systematics.append( ("TTBar_syst", "lnN", [ ("TTBar", TTBar_syst)] ))
				systematics.append( ("DY_syst", "lnN", [ ("DY", DY_syst)] ))
				systematics.append( ("TTBar_stat", "lnN", [ ("TTBar", TTBar_stat)] ))
				systematics.append( ("DY_stat", "lnN", [ ("DY", DY_stat)] ))
			systematics_list.append(systematics)
		except IOError:
			print mass, "not found"

	bg_names = ["TTBar", "DY"]

	for i in range(len(MWR)):
		print MWR[i], signal[i], sum(bg[i])
		print "syst", len(systematics_list[i])
		signal_tuple = (sig_name, signal[i])
		bg_tuples = zip(bg_names, bg[i])
		nBG = sum(bg[i])

		datacard = "WR%sjj_MASS%04d" % (channel,MWR[i])
		datacard_file = args.outdir + "/" + datacard + ".txt"
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
				signal_tuple, bg_tuples, systematics=systematics_list[i])
