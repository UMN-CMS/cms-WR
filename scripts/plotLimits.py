import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt
import ExoAnalysis.cmsWR.cross_sections as xs 
import sys 

channel = "ee"
name = sys.argv[1]
res_file = sys.argv[2]

plotter = plt.limit1d(.001)
plotter.addTheory(xs.WR_jj[channel])

with open(res_file, "r") as f:
	for line in f:
		if "COMBINE" in line:
			__, m, ret = eval(line)
			mass,mode = m.split('_')
			if mode == "EXPECTED":
				plotter.add(mass, ret)
			elif mode == "OBSERVED":
				plotter.addObserved(mass, ret)

plotter.plot("plots/limWR" + channel + "jj_" + name, x_title = "M_{W_{R}} [GeV]",
		y_title="Limit on XS(pb)", y_limits = (1e-3,1e-1), leg_y = .58 )
#plotter.plot("plots/limWR" + channel + ".png", x_title = "M_{W_{R}} [GeV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", y_range = (1e-3,10))
