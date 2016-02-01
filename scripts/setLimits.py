import ExoAnalysis.cmsWR.combineTools as combineTools

lumi = 2400.
for mass in sorted(combineTools.mass_cut):
	for channel in ["ee","mumu"]:
		for process in ["signal","TTBar","DY"]:
			print mass, channel, process, combineTools.getNEvents(mass,channel, process, lumi) 
