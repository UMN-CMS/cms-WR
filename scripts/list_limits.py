import sys
import ExoAnalysis.cmsWR.cross_sections as xs 
import ExoAnalysis.cmsWR.combineTools as combineTools

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

for i in range(len(MWR)):
	signal_tuple = (sig_name, signal[i])
	bg_tuples = zip(bg_names, bg[i])
	nBG = sum(bg[i])

	datacard = "WRmumujj_MASS%s" % MWR[i]
	datacard_file = "data/" + datacard + ".txt"
	sig, bgs = combineTools.makeDataCardSingleBin(datacard_file, "mumujj", nBG, signal_tuple, bg_tuples)
	method = "ProfileLikelihood"
	command = ["combine","-M",method,"-S","0",datacard_file,"-n", datacard,"-t","1000"]
	ret = combineTools.runCombine(command)
	(mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = ret
	print MWR[i], mean, xs.WR_mumujj[MWR[i]]

		
