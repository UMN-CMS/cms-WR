import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Analyzers

## analyze the kinematic distributions of the reco jets and leptons in WR->ENu->EEJJ evts with these analyzers

process.unmatchedSignalRecoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJet"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

process.unmatchedSignalRecoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsWithPtEtaCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

process.unmatchedSignalRecoAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsWithPtEtaAndDileptonMassCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.unmatchedSignalRecoAnalyzerFour = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("testUnmatchedSignalRecoObjectsWithPtEtaDileptonMassAndDrCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.unmatchedSignalRecoAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("testUnmatchedSignalRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

#process.investigateTrigger = cms.EDAnalyzer('triggerAnalyzer',
#		treeName = cms.string("electronChannelSignalHLT"),
#		commaSeparatedHltPaths = cms.string(),
#		trigResultsColl = cms.InputTag("TriggerResults","","HLT"),
#		trigObjectStandAloneColl = cms.InputTag("selectedPatTrigger")
#		)


#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.bareRecoParticleSeq
		*process.unmatchedSignalRecoAnalyzerOne
		*process.ptEtaRestrictedSeq
		*process.unmatchedSignalRecoAnalyzerTwo
		*process.recoDileptonCandidateSeq
		*process.unmatchedSignalRecoAnalyzerThree
		*process.recoDrSeparationSeq
		*process.unmatchedSignalRecoAnalyzerFour
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		)

#process.unmatchedSignalTriggerStudyPath = cms.Path(
#		process.investigateTrigger
#		)


process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analyzed_electronChannel_doubleEG_realData_midJuly_2015.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		##skims of real data! 13 TeV, 50ns  doubleEG July 2015
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_1.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_10.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_100.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_101.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_102.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_103.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_104.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_105.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_106.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_107.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_108.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_109.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_11.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_110.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_111.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_112.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_113.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_114.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_115.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_116.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_117.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_118.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_119.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_12.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_120.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_121.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_122.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_123.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_13.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_14.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_15.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_16.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_17.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_18.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_19.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_2.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_20.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_21.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_22.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_23.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_24.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_25.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_26.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_27.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_28.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_29.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_3.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_30.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_31.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_32.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_33.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_34.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_35.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_36.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_37.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_38.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_39.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_4.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_40.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_41.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_42.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_43.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_44.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_45.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_46.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_47.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_48.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_49.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_5.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_50.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_51.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_52.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_53.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_54.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_55.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_56.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_57.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_58.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_59.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_6.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_60.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_61.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_62.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_63.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_64.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_65.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_66.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_67.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_68.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_69.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_7.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_70.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_71.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_72.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_73.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_74.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_75.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_76.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_77.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_78.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_79.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_8.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_80.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_81.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_82.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_83.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_84.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_85.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_86.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_87.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_88.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_89.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_9.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_90.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_91.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_92.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_93.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_94.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_95.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_96.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_97.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_98.root',
        'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_13TeV_50ns_eejj_signalRegionSkim/150716_211704/0000/realData_electronSignalRegionSkim_99.root',       

		###WR->ENu->EEJJ 13TeV, 40 PU, MWR=2.6 TeV, MNu=1.3 TeV miniAOD files
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_9.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_2.root',

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


