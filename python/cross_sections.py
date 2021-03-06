from os import environ
folder = environ['LOCALRT'] + "/src/ExoAnalysis/cmsWR"


WR_eejj = {}
WR_mumujj = {}

WR_jj = {"ee":WR_eejj, "mumu":WR_mumujj}

with open(folder + "/configs/datasets.dat") as datasetsfile:
	for line in datasetsfile:
		line = line.split()
		name = line[0]
		namesplit = name.split("_")
		if len(namesplit) != 3: continue
		channel, MWR, MHNu, = namesplit
		if channel == "WRtoEEJJ":
			xs = float(line[2])
			WR_eejj[int(MWR)] = xs
		elif channel == "WRtoMuMuJJ":
			xs = float(line[2])
			WR_mumujj[int(MWR)] = xs

WR_eejj_offdiagonal = {}
WR_mumujj_offdiagonal = {}
WR_jj_offdiagonal = {"ee":WR_eejj_offdiagonal, "mumu":WR_mumujj_offdiagonal}

with open(folder + "/data/crosssections.txt","r") as xsfile:
	for line in xsfile:
		line = line.strip()
		if line[0] == "#": continue
		ch, mwr, mnu, xs, xs_err = line.split()
		WR_jj_offdiagonal[ch][(int(mwr),int(mnu))] = float(xs)


#WR_eejj = {
#800 :3.65,
#1000:1.78,
#1200:0.663,
#1400:0.389,
#1600:0.177,
#1800:0.117,
#2000:0.0707,
#2200:0.045,
#2400:0.0248,
#2600:0.015,
#2800:0.00913,
#3000:0.00576,
#3200:0.0034,
#3400:0.00263,
#3600:0.00154,
#3800:0.00119,
#4000:0.000801,
#4200:0.000473,
#4400:0.000375,
#4600:0.00019,
#4800:0.000152,
#5000:9.12e-05,
#5200:6.65e-05,
#5600:2.54e-05,
#5800:2.02e-05,
#6000:1.44e-05,
#}
#
#WR_mumujj = {
#800:4.06,
#1000:1.76,
#1200:0.717,
#1400:0.336,
#1600:0.174,
#1800:0.0897,
#2000:0.0622,
#2200:0.0415,
#2400:0.0276,
#2600:0.0142,
#2800:0.00913,
#3000:0.00524,
#3200:0.00442,
#3400:0.00285,
#3600:0.0016,
#3800:0.00107,
#4000:0.000702,
#4200:0.000508,
#4400:0.000326,
#4600:0.00019,
#4800:0.000172,
#5000:9.44e-05,
#5200:6.39e-05,
#5400:4.18e-05,
#5600:2.73e-05,
#5800:1.87e-05,
#6000:1.37e-05,
#}
