import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt
import ExoAnalysis.cmsWR.cross_sections as xs 
import os 
import sys
import pickle


configfile = open("configs/2015-v1.conf")
config = dict( [ line.strip().split('=') for line in configfile])

#prodSpace = "/local/cms/user/phansen/limits/"
#prodSpace = "/afs/cern.ch/work/p/phansen/public/wr/limits/" + tag + "/"

tag,prodSpace = sys.argv[1:] 

prodSpace = prodSpace + "/" + tag + "/"
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

def write_new_card(old_dc_filename, outfilename, obs):
	with open(old_dc_filename,'r') as dcf:
		dc = dcf.read()
	start = dc.find("observation")
	end = dc.find("\n",start)
	new_dc = dc.replace(dc[start:end], "observation %d" % obs)
	with open(outfilename,'w') as outf:
		outf.write(new_dc)

try:
	with open("obs_limit.p","r") as obsp:
		print "loading from pickle file"
		obs_limits = pickle.load(obsp)
except IOError:
	obs_limits = {"ee":[], "mumu":[]}
	for ch in ["ee","mumu"]:
		for obs in range(3):
			jobname =  ch + "obs" + str(obs)
			new_dc = "datacards/" + jobname + ".txt"
			write_new_card("datacards/WR{ch}jj_MASS5200.txt".format(ch=ch),new_dc, obs)
			obs_limits[ch].append( float(combineTools.runCombine("combine -M BayesianToyMC -n " +jobname + " " + new_dc, jobname)))
	with open("obs_limit.p","w") as obsp:
		pickle.dump(obs_limits, obsp)

plotters = {"ee":plt.limit1d("e",.001), "mumu":plt.limit1d("#mu",.001)}
plotters["ee"].addTheory(xs.WR_jj["ee"])
plotters["mumu"].addTheory(xs.WR_jj["mumu"])


for res in results:
	_, m, ret = res
	channel,mass,mode = m.split('_')
	if mode == "EXPECTED":
		plotters[channel].add(mass, ret)
	#elif mode == "OBSERVED":
		#plotters[channel].addObserved(mass, ret)

#full plot
plotters["ee"].plot("plots/limWReejj" + name + tag + "_full", x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_y = .56 )
plotters["mumu"].plot("plots/limWRmumujj" + name + tag + "_full", x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_y = .56 )

#zoomed plot
plotters["ee"].plot("plots/limWReejj" + name + tag, x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_x = .62, leg_h = .22, leg_y = .63, x_limits = (600,4000))
plotters["mumu"].plot("plots/limWRmumujj" + name + tag, x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_x = .62, leg_h = .22, leg_y = .63, x_limits = (600,4000))

#plot with obs lines
plotters["ee"].addObsLines(obs_limits["ee"])
plotters["mumu"].addObsLines(obs_limits["mumu"])

plotters["ee"].plot("plots/limWReejj" + name + tag + "_ObsLines", x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_y = .56 )
plotters["mumu"].plot("plots/limWRmumujj" + name + tag + "_ObsLines", x_title = "M_{W_{R}} [GeV]",
	y_title="XS(pb)", y_limits = (1e-3,1e-1), leg_y = .56 )
