import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.cross_sections as xs 

minitrees = combineTools.miniTreeInterface()

for channel in ["ee","mumu"]:
	sig_name = "WR_" + channel + "jj"
	MWR = []
	signal = []
	bg = []
	systematics_list  = []
	for mass in sorted(combineTools.mass_cut):
		try:
			#TODO: add Lumi uncertainty and others
			systematics = []
			signalNevents, sig_syst = minitrees.getNEvents(mass, channel, "signal") * .001/xs.WR_jj[channel][mass]

			TTBar,TTBar_syst = minitrees.getNEvents(mass, channel, "TT")
			DY,DY_syst = minitrees.getNEvents(mass, channel, "DYAMC")

			MWR.append(mass)
			signal.append(signalNevents)
			bg.append([TTBar, DY])
			systematics.append( ("Signal_unc", "lnN", [ (sig_name, sig_syst )]))
			systematics.append( ("TTBar_unc", "lnN", [ ("TTBar", TTBar_syst)] ))
			systematics.append( ("DY_unc", "lnN", [ ("DY", DY_syst)] ))
			systematics_list.append(systematics)
		except IOError:
			print mass, "not found"

	bg_names = ["TTBar", "DY"]

	for i in range(len(MWR)):
		print MWR[i], signal[i], sum(bg[i])
		signal_tuple = (sig_name, signal[i])
		bg_tuples = zip(bg_names, bg[i])
		nBG = sum(bg[i])

		datacard = "WR%sjj_MASS%04d" % (channel,MWR[i])
		datacard_file = thisdir + "/datacards/" + datacard + ".txt"
		sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
				signal_tuple, bg_tuples, systematics=systematics_list[i])
