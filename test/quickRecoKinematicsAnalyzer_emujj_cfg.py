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
		#checkThisHltPath = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1")
	
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
		fileName = cms.string('/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_TTJets_50ns_skim_low_mass_region_emujj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/MuonEG/analyzed_lowMassRegionSkim/analyzed_MuonEG_skim_low_mass_region_emujj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/SingleMuon/analyzed_lowMassRegionSkim/analyzed_SingleMuon_skim_low_mass_region_emujj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/SingleElectron/analyzed_lowMassRegionSkim/analyzed_SingleElectron_skim_low_mass_region_emujj.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

inputFiles = cms.untracked.vstring()
inputFiles.extend([
	##TTJets to EMuJJ 50ns low mass skim evts
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_1.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_10.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_11.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_12.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_13.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_14.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_15.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_16.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_17.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_18.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_19.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_2.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_20.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_21.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_22.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_23.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_24.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_25.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_26.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_27.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_28.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_29.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_3.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_30.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_31.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_32.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_33.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_34.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_35.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_36.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_37.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_38.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_39.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_4.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_40.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_41.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_42.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_43.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_44.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_45.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_46.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_47.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_48.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_49.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_5.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_50.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_51.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_52.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_54.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_55.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_56.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_57.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_58.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_59.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_6.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_60.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_61.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_62.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_63.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_64.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_7.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_8.root',
       'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EMuJJ_atFNALLPC/150803_125729/0000/miniAODEMuSkimLowMassRegion_9.root',
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


