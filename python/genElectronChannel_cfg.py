import FWCore.ParameterSet.Config as cms

process = cms.Process("GENEEJJ")

### define the filters and producers to be used at gen level

##filters to select gen jets and leptons without any selection requirements
process.bareGenJet = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("slimmedGenJets"),
		src = cms.InputTag("ak4GenJets"),
		cut = cms.string("")
		)

process.bareGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenJet"),
		minNumber = cms.uint32(2)
		)

process.bareGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("prunedGenParticles"),
		src = cms.InputTag("genParticles"),
		cut = cms.string("pdgId == 11 || pdgId == -11")
		)

process.bareGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenEle"),
		minNumber = cms.uint32(2)
		)

process.genAnalyzerOne = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsNoCuts"),
		genLeptonCollection = cms.InputTag("bareGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("bareGenJet","","GENEEJJ")
		)


##filters on pt and eta of gen leptons and jets
process.ptEtaRestrictedGenJet = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("slimmedGenJets"),
		src = cms.InputTag("bareGenJet","","GENEEJJ"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

process.ptEtaRestrictedGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedGenJet"),
		minNumber = cms.uint32(2)
		)

process.ptEtaRestrictedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareGenEle"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

process.ptEtaRestrictedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle"),
		minNumber = cms.uint32(2)
		)

process.ptEtaRestrictedLeadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		cut = cms.string("pt>60")
		)

process.ptEtaRestrictedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

process.genAnalyzerTwo = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

##filters on dilepton mass, and the producer necessary to facilitate this filter
process.genDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadingGenEle ptEtaRestrictedSubleadingGenEle"),
		role = cms.string("leading subleading"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

process.genDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

process.genAnalyzerThree = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaAndDileptonMassCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

process.runGenAnalysis = cms.Path(
		process.bareGenEle
		*process.bareGenEleFilter
		*process.bareGenJet
		*process.bareGenJetFilter
		*process.genAnalyzerOne
		
		*process.ptEtaRestrictedGenJet
		*process.ptEtaRestrictedGenJetFilter
		*process.ptEtaRestrictedSubleadingGenEle
		*process.ptEtaRestrictedSubleadingGenEleFilter
		*process.ptEtaRestrictedLeadingGenEle
		*process.ptEtaRestrictedLeadingGenEleFilter
		*process.genAnalyzerTwo

		*process.genDiElectronCandidate
		*process.genDiElectronCandidateFilter
		*process.genAnalyzerThree

		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analysis_genElectronChannel.root')
)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		'file:/uscms/home/skalafut/nobackup/WR_starting2015/GEN-SIM_13TeV/7A8730D6-9E86-E411-A442-002590747DF0.root',

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)






