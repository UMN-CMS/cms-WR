import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt
import ExoAnalysis.cmsWR.cross_sections as xs 
import os 
import sys


configfile = open("configs/2015-v1.conf")
config = dict( [ line.strip().split('=') for line in configfile])

#prodSpace = "/local/cms/user/phansen/limits/"
#prodSpace = "/afs/cern.ch/work/p/phansen/public/wr/limits/" + tag + "/"

tag,prodSpace sys.argv[1:] 
name = config["productionTAG"]

results = []
for log in os.listdir(prodSpace + name):
	if log[-4:] != ".log": continue
	with open(prodSpace + name + "/" + log) as f:
		for line in f:
			if "COMBINE" in line:
				print line
				results.append(eval(line.strip()))
results.sort()

plotters = {"ee":plt.limit1d(.001), "mumu":plt.limit1d(.001)}
plotters["ee"].addTheory(xs.WR_jj["ee"])
plotters["mumu"].addTheory(xs.WR_jj["mumu"])

for res in results:
	_, m, ret = res
	channel,mass,mode = m.split('_')
	if mode == "EXPECTED":
		plotters[channel].add(mass, ret)
	#elif mode == "OBSERVED":
		#plotters[channel].addObserved(mass, ret)

plotters["ee"].plot("plots/limWReejj" + name + tag, x_title = "M_{W_{R}} [GeV]",
	y_title="Limit on XS(pb)", y_limits = (1e-3,1e-1), leg_y = .58 )
plotters["mumu"].plot("plots/limWRmumujj" + name + tag, x_title = "M_{W_{R}} [GeV]",
	y_title="Limit on XS(pb)", y_limits = (1e-3,1e-1), leg_y = .58 )
#plotters.plot("plots/limWR" + channel + ".png", x_title = "M_{W_{R}} [GeV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", y_range = (1e-3,10))
