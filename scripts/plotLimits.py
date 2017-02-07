import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt
from ExoAnalysis.cmsWR.PlotUtils import ave_tuple
import ExoAnalysis.cmsWR.cross_sections as xs 
import os 
import sys
import pickle


configfile = open("configs/2015-v1.conf")
config = dict( [ line.strip().split('=') for line in configfile])

tag = sys.argv[1] 

#prodSpace = "/local/cms/user/phansen/limits/"
#prodSpace = "/afs/cern.ch/work/p/phansen/public/wr/limits/" + tag + "/"
prodSpace = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/limitSetting/" + tag + "/"
name = config["productionTAG"]

results = []
for log in os.listdir(prodSpace + name):
	if log[-4:] != ".log": continue
	with open(prodSpace + name + "/" + log) as f:
		for line in f:
			if "COMBINE" in line:
				results.append(eval(line.strip()))

ee_interp = []
mumu_interp = []
for res in results:
	_, m, ret = res
	if m == "ee_3200_EXPECTED":
		ee_interp.append(ret)
	elif m == "ee_3600_EXPECTED":
		ee_interp.append(ret)
	elif m == "mumu_3200_EXPECTED":
		mumu_interp.append(ret)
	elif m == "mumu_3600_EXPECTED":
		mumu_interp.append(ret)

if len(ee_interp) > 1:
	results.append(("COMBINE", "ee_3400_EXPECTED", ave_tuple(*ee_interp)))
if len(mumu_interp) > 1:
	results.append(("COMBINE", "mumu_3400_EXPECTED", ave_tuple(*mumu_interp)))

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
			print "calculating Observed = ", obs
			jobname =  ch + "obs" + str(obs)
			new_dc = "datacards/" + jobname + ".txt"
			write_new_card("datacards/WR{ch}jj_MASS5200.txt".format(ch=ch),new_dc, obs)
			obs_limits[ch].append( float(combineTools.runCombine("combine -M BayesianToyMC -n " +jobname + " " + new_dc, jobname)))
	with open("obs_limit.p","w") as obsp:
		pickle.dump(obs_limits, obsp)

plotters = {"ee":plt.limit1d("e",.001), "mumu":plt.limit1d("#mu",.001)}
plotters["ee"].addTheory(xs.WR_jj["ee"])
plotters["mumu"].addTheory(xs.WR_jj["mumu"])

plotters2d = {
		"ee"  : plt.limit2d("data/offlineEEEfficienciesVsMassesNoGenNuFilterWithMassWindowCuts.txt",    "ee2d" , "ee",     (17, 700, 4100), (40, 50, 4050), xs=.001),
		"mumu": plt.limit2d("data/offlineMuMuEfficienciesVsMassesNoGenNuFilterWithMassWindowCuts.txt", "mumu2d", "#mu#mu", (17, 700, 4100), (40, 50, 4050), xs=.001)
		}

plotters2d["ee"].addTheory(xs.WR_jj_offdiagonal["ee"])
plotters2d["mumu"].addTheory(xs.WR_jj_offdiagonal["mumu"])

for res in results: 
	_, m, ret = res
	channel,mass,mode = m.split('_')
	if mode == "EXPECTED":
		print m, ret
		plotters[channel].add(mass, ret)
		plotters2d[channel].add(mass,ret)
	#elif mode == "OBSERVED":
		#plotters[channel].addObserved(mass, ret)


ytitle = "#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow %sqq)" 
#full plot
#plotters["ee"].plot("plots/limWReejj" + name + tag + "_full", x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "ee", y_limits = (1e-3,1e-1), leg_y = .56 )
#plotters["mumu"].plot("plots/limWRmumujj" + name + tag + "_full", x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "#mu#mu", y_limits = (1e-3,1e-1), leg_y = .56 )

#zoomed plot
#plotters["ee"].plot("plots/limWReejj" + name + tag, x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "ee", y_limits = (1e-3,1e-1), leg_x = .62, leg_h = .22, leg_y = .63, x_limits = (600,4000))
#plotters["mumu"].plot("plots/limWRmumujj" + name + tag, x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "#mu#mu", y_limits = (1e-3,1e-1), leg_x = .62, leg_h = .22, leg_y = .63, x_limits = (600,4000))

#plot with obs lines
#plotters["ee"].addObsLines(obs_limits["ee"])
#plotters["mumu"].addObsLines(obs_limits["mumu"])

#plotters["ee"].plot("plots/limWReejj" + name + tag + "_ObsLines", x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "ee", y_limits = (1e-3,1e-1), leg_y = .56 )
#plotters["mumu"].plot("plots/limWRmumujj" + name + tag + "_ObsLines", x_title = "M_{W_{R}} [GeV]",
#	y_title=ytitle % "#mu#mu", y_limits = (1e-3,1e-1), leg_y = .56 )

plotters2d["ee"].plot("plots/lim2dWReejj" + name + tag)
plotters2d["mumu"].plot("plots/lim2dWRmumujj" + name + tag)

