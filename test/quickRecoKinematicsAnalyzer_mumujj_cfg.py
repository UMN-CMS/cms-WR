import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOMuMuJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoMuonChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates in bkgnd evts with these analyzers
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("mumuBareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
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
		process.mumuBareRecoParticleSeq
		*process.mumuBareRecoDrSeparationSeq
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
		#fileName = cms.string('analyzed_SMPL_skim_low_mass_region_mumujj.root')
		fileName = cms.string('/eos/uscms/store/user/skalafut/DYJetsSkimmedEEJJ_13TeV_25ns_mumujj_skim_low_mass_region/analyzed_DYJets_skim_low_mass_region_mumujj.root')



)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)


inputFiles = cms.untracked.vstring()
inputFiles.extend([

	])

	
process.source = cms.Source( "PoolSource",
    fileNames = inputFiles, 
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


