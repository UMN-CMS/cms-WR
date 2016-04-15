import ROOT
import subprocess

prefix = "root://eoscms//eos/cms/store/user/shervin/ntuples/"

eos = "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select"

command = eos + " ls " + prefix
datasets = subprocess.check_output(command.split())

configfile = open("configs/2015-v1.conf")
config = dict( [ line.strip().split('=') for line in configfile])

tree_channels = [
		"miniTree_flavoursideband",
		"miniTree_lowdileptonsideband",
		"miniTree_signal_ee",
		"miniTree_signal_mumu",
		"miniTree_dytagandprobe",
		]


print "# Dataset", "%s "*len(tree_channels) % tuple(tree_channels)
for dataset in datasets.split():
	if config["productionTAG"] in dataset:
		f = ROOT.TFile.Open(prefix + dataset + '/unmerged-allRange.root')
		if f:
			print dataset,
			for tree_channel in tree_channels:
				t = f.Get(tree_channel + "/t")
				if t:
					print t.GetEntries(),
				else: 
					print "NULL",
			print 
		else:
			print "#",dataset, "dataset not found"
