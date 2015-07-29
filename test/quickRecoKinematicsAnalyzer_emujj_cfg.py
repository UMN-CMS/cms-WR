import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEMuJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoEMuChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates in bkgnd evts with these analyzers
#use the matchedAnalyzer plugin to work with separate electron and muon collections
process.recoAnalyzerOne = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("emuBareRecoLepton"),
		subleadingLeptonCollection = cms.InputTag("emuBareRecoLeptonTwo"),
		quarkCollection = cms.InputTag("emuBareRecoJet"),
		matchingLeadingLeptonCollection = cms.InputTag(""),
		matchingSubleadingLeptonCollection = cms.InputTag(""),
		matchingQuarkCollection = cms.InputTag("")
		)
#
#process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaCuts"),
#		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
#		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
#		doDileptonMassCut = cms.bool(False),
#		minDileptonMass = cms.double(-1)
#		)
#
#process.recoAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaAndDileptonMassCuts"),
#		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
#		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)
#
#process.recoAnalyzerFour = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaDileptonMassAndDrCuts"),
#		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
#		jetsCollection = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)
#
#process.recoAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
#		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
#		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)


#################################
#Paths
process.unmatchedBkgndRecoPath = cms.Path(
		process.emuBareRecoParticleSeq
		*process.recoAnalyzerOne
		#*process.ptEtaRestrictedSeq
		#*process.recoAnalyzerTwo
		#*process.recoDileptonCandidateSeq
		#*process.recoAnalyzerThree
		#*process.recoDrSeparationSeq
		#*process.recoAnalyzerFour
		#*process.recoFourObjMassSeq
		#*process.recoAnalyzerFive
		)
process.schedule = cms.Schedule(process.unmatchedBkgndRecoPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analyzed_SMPL_skim_low_mass_region_emujj.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
       'file:overwrite.root'
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


