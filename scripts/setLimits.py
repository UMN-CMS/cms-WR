import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.condorTools as condorTools
import ExoAnalysis.cmsWR.cross_sections as xs 
import ExoAnalysis.cmsWR.plot as plt
import math

import sys, os

mode = sys.argv[1]
syst_unc = float(sys.argv[2])+1
name = sys.argv[3]

lumi = 2400.
channel = "ee"
sig_name = "WR_" + channel + "jj"
MWR = []
signal = []
bg = []
systematics_list  = []
for mass in sorted(combineTools.mass_cut):
	try:
		systematics = []
		signalNevents = combineTools.getNEvents(mass, channel, "signal", lumi) * .001/xs.WR_jj[channel][mass]
		TTBar = combineTools.getNEvents(mass, channel, "TTBar", lumi)
		DY = combineTools.getNEvents(mass, channel, "DY", lumi)

		MWR.append(mass)
		signal.append(signalNevents)
		bg.append([TTBar, DY])
		systematics.append( ("Signal_unc", "lnN", [ (sig_name,syst_unc )]))
		systematics.append( ("TTBar_unc", "lnN", [ ("TTBar", syst_unc)] ))
		systematics.append( ("DY_unc", "lnN", [ ("DY", syst_unc)] ))
		systematics_list.append(systematics)
	except IOError:
		print mass, "not found"


bg_names = ["TTBar", "DY"]

plotter = plt.limit1d()
#plotter.addTheory(xs.WR_jj[channel])

thisdir = os.getcwd()

job = condorTools.Job(thisdir + "/scripts/batch_run", name, prodSpace="/local/cms/user/phansen/limits/")
for i in range(len(MWR)):
	print MWR[i], signal[i], sum(bg[i])
	signal_tuple = (sig_name, signal[i])
	bg_tuples = zip(bg_names, bg[i])
	nBG = sum(bg[i])

	datacard = "WR%sjj_MASS%s_%s" % (channel,MWR[i], name)
	datacard_file = thisdir + "/data/" + datacard + ".txt"
	#systematics =  [ ("Signal_unc", "lnN", [ (sig_name, XX)])]
	#systematics.append( "TTBar_unc", "lnN", [ ("TTBar", XX)] )
	#systematics.append( "DY_unc", "lnN", [ ("DY", XX)] )
	sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, channel + "jj", nBG,
			signal_tuple, bg_tuples, systematics=systematics_list[i])

	systematics = True
	jobname = "%s_EXPECTED" % MWR[i]
	command = "combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood -S%d %s -n %s -T 5000 " %(systematics, datacard_file, datacard)
	prefix  = thisdir + "/python/combineTools.py " + jobname
	job.addJob( prefix + " " + command, jobname)
	#jobname = "%s_OBSERVED" % MWR[i]
	#command = "combine -M BayesianToyMC -H ProfileLikelihood -S%d %s -n %s " %(systematics, datacard_file, datacard)
	#prefix  = thisdir + "/python/combineTools.py " + jobname
	#job.addJob( prefix + " " + command, jobname)
	#jobname = "%s_EXPECTED" % MWR[i]
	#command = "combine -M BayesianToyMC -H ProfileLikelihood -S%d %s -n %s --toys 100" %(systematics, datacard_file, datacard)
	#prefix  = thisdir + "/python/combineTools.py " + jobname
	#job.addJob( prefix + " " + command, jobname)
	#ret = combineTools.runCombine(command)
	#print ret
	#plotter.add(MWR[i], ret)
	#median, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
	#print MWR[i], median, onesig_minus, onesig_plus, twosig_minus,twosig_plus, xs.WR_jj[channel][MWR[i]]

#job.submit()

#plotter.plot("plots/limWR" + channel + "jj_" + name, x_title = "M_{W_{R}} [GeV]",
		#y_title="Limit/Expected ", y_limits = (1e-3,10), leg_y = .18 )
#plotter.plot("plots/limWR" + channel + ".png", x_title = "M_{W_{R}} [GeV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", y_range = (1e-3,10))
