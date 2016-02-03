import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.cross_sections as xs 
import ExoAnalysis.cmsWR.plot as plt
import math

import sys

syst_unc = float(sys.argv[1])+1
name = sys.argv[2]

lumi = 2400.
channel = "ee"
sig_name = "WR_" + channel + "jj"
MWR = []
signal = []
bg = []
systematics_list  = []
for mass in sorted(combineTools.mass_cut):
	systematics = []
	MWR.append(mass)
	signalNevents = combineTools.getNEvents(mass, channel, "signal", lumi)
	signal.append(signalNevents)
	systematics.append( ("Signal_stat", "lnN", [ (sig_name,syst_unc )]))

	TTBar = combineTools.getNEvents(mass, channel, "TTBar", lumi)
	DY = combineTools.getNEvents(mass, channel, "DY", lumi)
	systematics.append( ("TTBar_stat", "lnN", [ ("TTBar", syst_unc)] ))
	systematics.append( ("DY_stat", "lnN", [ ("DY", syst_unc)] ))
	bg.append([TTBar, DY])
	systematics_list.append(systematics)

bg_names = ["TTBar", "DY"]

plotter = plt.limit1d()
#plotter.addTheory(xs.WR_jj[channel])

for i in range(len(MWR)):
	signal_tuple = (sig_name, signal[i])
	bg[i] = [x for x in bg[i]]
	bg_tuples = zip(bg_names, bg[i])
	nBG = sum(bg[i])

	datacard = "WR%sjj_MASS%s_%s" % (channel,MWR[i], name)
	datacard_file = "data/" + datacard + ".txt"
	#systematics =  [ ("Signal_stat", "lnN", [ (sig_name, XX)])]
	#systematics.append( "TTBar_stat", "lnN", [ ("TTBar", XX)] )
	#systematics.append( "DY_stat", "lnN", [ ("DY", XX)] )
	sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
			signal_tuple, bg_tuples, systematics=systematics_list[i])
	method = "ProfileLikelihood"
	command = ["combine","-M",method,datacard_file,"-n", datacard,"-t","1000"]
	command = ["combine","-M",method,"-H","ProfileLikelihood","-S0",datacard_file,"-n", datacard,"-t","1000000"]
	ret = combineTools.runCombine(command)
	plotter.add(MWR[i], ret)
	median, (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
	print MWR[i], median, mean, meanError, onesig_minus, onesig_plus, twosig_minus,twosig_plus, xs.WR_jj[channel][MWR[i]]

plotter.plot("plots/limWR" + channel + "jj_" + name, x_title = "M_{W_{R}} [GeV]",
		y_title="Limit/Expected ", y_limits = (1e-3,10), leg_y = .18 )
#plotter.plot("plots/limWR" + channel + ".png", x_title = "M_{W_{R}} [GeV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", y_range = (1e-3,10))
