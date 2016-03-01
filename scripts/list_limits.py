import sys
import ExoAnalysis.cmsWR.cross_sections as xs 
import ExoAnalysis.cmsWR.combineTools as combineTools
import ExoAnalysis.cmsWR.plot as plt

infile = sys.argv[1]
header = []
MWR = []
signal = []
bg = []
with open(infile,"r") as f:
	for line in f:
		if not header:
			header = line.split()
			continue
		l = line.split()
		MWR.append(int(l[0]))
		signal.append(float(l[1]))
		bg.append(map(float,l[2:]))

sig_name = header[1]
bg_names = header[2:]

plotter = plt.limit1d(0.01)
plotter.addTheory(xs.WR_mumujj)

for i in range(len(MWR)):
	signal_tuple = (sig_name, signal[i])
	bg[i] = [x*2 for x in bg[i]]
	bg_tuples = zip(bg_names, bg[i])
	nBG = sum(bg[i])

	datacard = "WRmumujj_MASS%s" % MWR[i]
	datacard_file = "data/" + datacard + ".txt"
	sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, "mumujj", nBG, signal_tuple, bg_tuples)
	method = "ProfileLikelihood"
	command = ["combine","-M",method,"-S","0",datacard_file,"-n", datacard,"-t","1000"]
	ret = combineTools.runCombine(command)
	plotter.add(MWR[i], ret)
	median, (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
	print MWR[i], median, mean, meanError, onesig_minus, onesig_plus, twosig_minus,twosig_plus, xs.WR_mumujj[MWR[i]]

plotter.plot("limWRmumu.png", x_title = "M_{W_{R}} [TEV]", y_title="#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mu) [fb]", )

		
