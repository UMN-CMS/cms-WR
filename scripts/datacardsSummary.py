import ExoAnalysis.cmsWR.cross_sections as xs
import os
import re
import math

datacards = os.listdir("datacards/")

p = re.compile("WR(ee|mumu)jj_MASS(\d*).txt")
filenames = []
for datacard in datacards:
	m = p.match(datacard)
	try:
		ch, mass = m.groups()
	except:
		continue
	filenames.append(( ch, mass, 'datacards/' + datacard))

print "#CHANNEL MASS XS(fb) SIG_NEVENTS SIG_RATE TT_RATE DY_RATE SIG_SD TT_SD DY_SD TT_SF DY_SF SIG_N TT_N DY_N SIG_ALPHA TT_ALPHA DY_ALPHA"
for ch,mass,fn in sorted(filenames):
	with open(fn) as f:
		unc = {}
		for line in f:
			if "rate" in line:
				__, sig, tt, dy = line.split()
			elif "SF" in line:
				line = line.split()
				name = line[0]
				unc[name] = float([x for x in line[2:] if '-' != x][0])
			elif "gmN" in line:
				line = line.split()
				name = line[0]
				N = int(line[2])
				unc[name] = float([x for x in line[3:] if '-' != x][0])
				unc[name + "N"] = N

		print ch,mass, xs.WR_jj[ch][int(mass)]/.001, "%.3f" % (float(xs.WR_jj[ch][int(mass)]/.001)*float(sig)), 
		print sig, tt, dy,
		print "%.4f" % (math.sqrt((unc["signal_uncN"]+1))*unc["signal_unc"]),
		print "%.4f" % (math.sqrt((unc["TT_uncN"]    +1))*unc["TT_unc"]    ),
		print "%.4f" % (math.sqrt((unc["DYAMC_uncN"] +1))*unc["DYAMC_unc"] ),
		print unc["TT_SF"],
		print unc["DYAMC_SF"],
		print unc["signal_uncN"],
		print unc["TT_uncN"],
		print unc["DYAMC_uncN"],
		print unc["signal_unc"],
		print unc["TT_unc"],
		print unc["DYAMC_unc"],
		print



	
