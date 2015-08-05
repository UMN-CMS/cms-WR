import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Filters
process.trigSelector = cms.EDFilter("triggerFilter",
		checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2")
		)


#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates in bkgnd evts with these analyzers
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		#checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2")
		)

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
		process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
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
		fileName = cms.string('/eos/uscms/store/user/skalafut/DoubleEG/analyzed_lowMassRegionSkim/analyzed_DoubleEG_skim_low_mass_region_eejj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/SingleElectron/analyzed_lowMassRegionSkim/analyzed_SingleElectron_skim_low_mass_region_eejj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_TTJets_50ns_skim_low_mass_region_eejj.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

inputFiles = cms.untracked.vstring()
inputFiles.extend([
	   
	   #################################################################
	   #TTJets to EEJJ 50ns low mass skimmed evts
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_1.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_10.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_11.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_12.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_13.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_14.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_15.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_16.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_17.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_18.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_2.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_20.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_21.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_22.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_23.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_24.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_25.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_26.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_27.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_29.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_3.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_30.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_31.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_32.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_33.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_34.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_35.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_36.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_37.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_38.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_39.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_4.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_40.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_42.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_43.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_44.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_45.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_46.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_47.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_48.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_49.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_5.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_50.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_51.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_52.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_53.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_54.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_55.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_56.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_57.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_58.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_59.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_6.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_60.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_61.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_62.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_63.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_64.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_7.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_8.root',
       #'file:/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTJets_13TeV_50ns_skim_low_mass_and_signal_regions_EEJJ_atFNALLPC/150803_125706/0000/miniAODEleEleSkimLowMassRegion_9.root',

	   ##DoubleEG 50ns low mass skim evts
	   'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronLowMassRegionSkim_1.root'

	   ##SingleElectron 50ns low mass skim evts
	   #'file:/eos/uscms/store/user/skalafut/SingleElectron/realData_SingleElectron_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_143037/0000/realData_electronLowMassRegionSkim_1.root'


	])
	#end inputFiles

process.source = cms.Source( "PoolSource",
    #fileNames = cms.untracked.vstring(
    #   #'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronLowMassRegionSkim_1.root'
    #   #'file:/eos/uscms/store/user/skalafut/SingleElectron/realData_SingleElectron_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_143037/0000/realData_electronLowMassRegionSkim_1.root'
    #   
	#   
	#    ),
	fileNames = inputFiles,
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


