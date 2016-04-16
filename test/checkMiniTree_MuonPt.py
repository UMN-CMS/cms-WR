import ROOT

configfile = open("configs/2015-v1.conf")
config = dict( [ line.strip().split('=') for line in configfile])

datasets = {"DYJets_amctnlo":"", "DYJets_madgraph":"", "DYToEE_powheg":""}
for dataset in datasets:
	print dataset
	f = ROOT.TFile.Open("root://eoscms//eos/cms/store//user/shervin/ntuples/%s%s/unmerged-allRange.root" % (dataset, config["productionTAG"]))
	tree = f.Get("miniTree_dytagandprobe/t")
	for event in tree:
		for mu_p4 in event.muons_p4:
			if mu_p4.Pt() > 1000.0:
				print event.event, mu_p4.Pt()

