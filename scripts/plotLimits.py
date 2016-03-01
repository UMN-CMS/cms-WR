import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt
import ExoAnalysis.cmsWR.cross_sections as xs 
import sys 

channel = "ee"
name = sys.argv[1]
file_format = sys.argv[2]

plotter = plt.limit1d(.001)
plotter.addTheory(xs.WR_jj[channel])

obs_file_format = file_format + "_OBSERVED.log"
exp_file_format = file_format + "_EXPECTED.log"
for mass in sorted(combineTools.mass_cut):
	filename =  exp_file_format % mass
	with open(filename, "r") as f:
		for line in f:
			if "COMBINE" in line:
				__, m, ret = eval(line)
				m = int(m.split('_')[0])
				plotter.add(m, ret)
	filename = obs_file_format % mass
	with open(filename, "r") as f:
		for line in f:
			if "COMBINE" in line:
				__, m, obs = eval(line)
				m = int(m.split('_')[0])
				#plotter.addObserved(m, obs)



plotter.plot("plots/limWR" + channel + "jj_" + name, x_title = "M_{W_{R}} [GeV]",
		y_title="Limit on XS(pb)", y_limits = (1e-4,1e-1), leg_y = .58 )
#plotter.plot("plots/limWR" + channel + ".png", x_title = "M_{W_{R}} [GeV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", y_range = (1e-3,10))
