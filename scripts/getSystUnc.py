from ExoAnalysis.cmsWR.combineTools import AnalysisResultsInterface

onlyjet = AnalysisResultsInterface(base="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysOnlyJetScaleSyst/")
onlyele = AnalysisResultsInterface(base="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysSmearEleScaleSyst/")
onlymu  = AnalysisResultsInterface(base="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysOnlyMuScaleIdSyst/")
allsyst = AnalysisResultsInterface(base="/afs/cern.ch/work/p/phansen/public/wr/4000toys/")


for interface,name,channels in zip([onlyjet, onlyele, onlymu, allsyst], ["Jet", "Electron", "Muon", "All"], [["ee","mumu"],["ee"],["mumu"],["ee","mumu"]]):
	for channel in channels:
		for mode in ["signal", "DYAMC"]:
			for mass in [800, 1000, 2000, 2800, 3200, 3800]:
				f = interface.OpenFile(channel, mode, mass)
				if not  interface.masses: interface.GetMasses(f)
				key = channel + "_" + mode
				if "signal" == mode:
					key += "_" + str(mass)

				interface.ProcessFile(key, f)
				print name, channel, mode, mass, interface.results[key]["syst"]["std"][interface.masses.index(mass)]



