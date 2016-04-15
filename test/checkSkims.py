import ROOT
import subprocess
import sys
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()

import csv 
config = csv.DictReader(open("configs/datasets.dat","r"), delimiter='\t')

handleMuons = Handle('std::vector<pat::Muon>')
muonTag = "slimmedMuons"
handleVertex = Handle('std::vector<reco::Vertex>')
vertexTag = "offlineSlimmedPrimaryVertices"
handleGen = Handle('std::vector<reco::GenParticle>')
genTag = "prunedGenParticles"

datasets = {"DYJets_amctnlo":"", "DYJets_madgraph":"", "DYToEE_powheg":""}
for line in config:
	if line["#name"] in datasets:
		datasets[line["#name"]] = line["skimmed_dataset"]

for dataset, eosdir in datasets.items():
	print dataset
	subprocess.call(("makefilelist.sh " + dataset + " " + eosdir).split())

	files = []
	with open("filelist/" + dataset + ".list","r") as filelist:
		for line in filelist:
			files.append(line.strip())
	events = Events(files)
	print "N Entries in ",dataset, ":",  events.size()
	print "%10s | %8s | %8s | %8s | %8s | %8s | %8s | %6s" %("event", "muon pt", "inner pt", "outer pt", "global pt", "max gen pt", "isHighPt", "isGood")
	muonformat =  "%10d | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f | %8d | %6d"
	for event in events:
		event.getByLabel(muonTag, handleMuons)
		muons = handleMuons.product()
		event.getByLabel(vertexTag, handleVertex)
		vertices = handleVertex.product()
		for mu in muons:
			if mu.isHighPtMuon(vertices.at(0)) and mu.pt() > 1000:
				event.getByLabel(genTag, handleGen)
				genParticles = handleGen.product()
				genPT = [ p.pt() for p in genParticles if abs(p.pdgId()) == 13]
				if genPT:
					maxgenPT = max(genPT)
				else:
					maxgenPT = 0.0
					#if p.status() == 1 and abs(p.pdgId()) == 13 and abs(p.eta() - mu.eta()) < .1:

				innerpt = mu.innerTrack().pt() if mu.isAValidMuonTrack(mu.InnerTrack) else 0.
				outerpt = mu.outerTrack().pt() if mu.isAValidMuonTrack(mu.OuterTrack) else 0.
				globalpt = mu.globalTrack().pt() if mu.isAValidMuonTrack(mu.CombinedTrack) else 0.

				print muonformat % (event.object().id().event(), mu.pt(), innerpt, outerpt, globalpt, maxgenPT, mu.isHighPtMuon(vertices.at(0)), mu.isGood("All"))
