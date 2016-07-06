import ExoAnalysis.cmsWR.cross_sections as xs
import os
import re

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

print "#CHANNEL MASS XS(fb) SIG_RATE TT_RATE DY_RATE SIG_N TT_N DY_N SIG_ALPHA TT_ALPHA DY_ALPHA"
for ch,mass,fn in sorted(filenames):
	with open(fn) as f:
		print ch,mass, xs.WR_jj[ch][int(mass)]/.001,
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
		print sig, tt, dy,
		print unc["TT_SF"],
		print unc["DYAMC_SF"],
		print unc["signal_uncN"],
		print unc["TT_uncN"],
		print unc["DYAMC_uncN"],
		print unc["signal_unc"],
		print unc["TT_unc"],
		print unc["DYAMC_unc"],
		print



	
