import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEMuJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoEMuChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Filters
process.trigSelector = cms.EDFilter("triggerFilter",
		checkThisHltPath = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1")
		)


#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates using these modules
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("recoObjectsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("emuBareRecoLeptonOne"),
		leptonsTwoCollection = cms.InputTag("emuBareRecoLeptonTwo"),
		jetsCollection = cms.InputTag("emuBareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		#checkThisHltPath = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2")
	
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
		*process.emuBareRecoDrSeparationSeq
		*process.trigSelector
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
		#fileName = cms.string('/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_TTJets_skim_low_mass_region_emujj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/MuonEG/analyzed_lowMassRegionSkim/analyzed_MuonEG_skim_low_mass_region_emujj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/SingleMuon/analyzed_lowMassRegionSkim/analyzed_SingleMuon_skim_low_mass_region_emujj.root')
		fileName = cms.string('/eos/uscms/store/user/skalafut/SingleElectron/analyzed_lowMassRegionSkim/analyzed_SingleElectron_skim_low_mass_region_emujj.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

inputFiles = cms.untracked.vstring()
inputFiles.extend([
	   #################################################################
	   #TTJets low mass skimmed evts


	   ##end TTJets evts which passed low mass skim requirements

	   ##MuonEG 50ns real data evts which passed emujj low mass skim requirements
	   #'file:/eos/uscms/store/user/skalafut/MuonEG/realData_MuonEG_13TeV_50ns_emujj_signalAndLowMassRegionSkim_atFNALLPC_redo/150802_220304/0000/realData_emuLowMassRegionSkim_1.root'

	   ##SingleMuon 50ns real data evts which passed emujj low mass skim requirements
	   #'file:/eos/uscms/store/user/skalafut/SingleMuon/realData_SingleMuon_13TeV_50ns_emujj_signalAndLowMassRegionSkim_atFNALLPC_redo/150802_220553/0000/realData_emuLowMassRegionSkim_1.root',
	

	   ##SingleElectron 50ns real data evts which passed emujj low mass skim requirements
	   'file:/eos/uscms/store/user/skalafut/SingleElectron/realData_SingleElectron_13TeV_50ns_emujj_signalAndLowMassRegionSkim_atFNALLPC_redo/150802_220520/0000/realData_emuLowMassRegionSkim_1.root',
	
	   ])
#end inputFiles
 
process.source = cms.Source( "PoolSource",
    fileNames = inputFiles, 
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


