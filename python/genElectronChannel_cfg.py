import FWCore.ParameterSet.Config as cms

process = cms.Process("GENEEJJ")

## load the filters, producers, and sequences defined in the config file fragment
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

## run these analyzers after selecting gen electrons, jets, and quarks without any cuts applied
process.genAnalyzerOne = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsNoCuts"),
		genLeptonCollection = cms.InputTag("bareGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("bareGenJet","","GENEEJJ")
		)

process.matchedGenAnalyzerOne = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedGenObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle","","GENEEJJ"),
		genSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle","","GENEEJJ"),
		genQuarkCollection = cms.InputTag("bareMatchedGenQuark","","GENEEJJ")
		)


## run these analyzers after applying pt and eta cuts on the leptons and jets
process.genAnalyzerTwo = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

process.matchedGenAnalyzerTwo = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ")
		)


## run these analyzers after applying pt, eta, and dilepton mass cuts
process.genAnalyzerThree = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaAndDileptonMassCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

process.matchedGenAnalyzerThree = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaAndDileptonMassCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ")
		)

## analyzers to run with dR(lepton,hadron) > someValue cut applied
process.matchedGenAnalyzerFour = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaDileptonMassAndDrCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ")
		)

## analyzers to run with dR(lepton,hadron) > someValue cut and four obj inv mass cut applied
process.matchedGenAnalyzerFive = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(True),
		minFourObjMass = cms.double(600),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ")
		)



#process.runUnmatchedGenAnalysis = cms.Path(
#		process.bareGenParticleSeq
#		*process.genAnalyzerOne
#		*process.ptEtaConstrainedSeq
#		*process.genAnalyzerTwo
#		*process.genDiElectronCandidateSeq
#		*process.genAnalyzerThree
#		)

process.runMatchedGenAnalysis = cms.Path(
		process.bareMatchedGenParticleSeq
		*process.matchedGenAnalyzerOne
		*process.matchedPtEtaConstrainedSeq
		*process.matchedGenAnalyzerTwo
		*process.genMatchedDiElectronCandidateSeq
		*process.matchedGenAnalyzerThree
		*process.matchedGenAnalyzerFour
		*process.matchedGenAnalyzerFive
		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analysis_genElectronChannel_using_prunedGenParticles.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		#'file:/uscms/home/skalafut/nobackup/WR_starting2015/GEN-SIM_13TeV/7A8730D6-9E86-E411-A442-002590747DF0.root',
		#'file:/uscms/home/skalafut/nobackup/WR_starting2015/GEN-SIM_13TeV/00774ADC-8F87-E411-B6EF-00266CF32E78.root'
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_1.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_2.root',

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)




