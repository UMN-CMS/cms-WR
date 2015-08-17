import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Filters
process.trigSelector = cms.EDFilter("triggerFilter",
		checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2"),
		alsoCheckThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v2"),

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
		)

process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsWithPtEtaCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

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
		#process.trigSelector
		process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		*process.recoAnalyzerOne
		*process.ptEtaRestrictedSeq
		*process.recoAnalyzerTwo
		#*process.recoDileptonCandidateSeq
		#*process.recoAnalyzerThree
		#*process.recoDrSeparationSeq
		#*process.recoAnalyzerFour
		#*process.recoFourObjMassSeq
		#*process.recoAnalyzerFive
		)
process.schedule = cms.Schedule(process.unmatchedBkgndRecoPath)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('/eos/uscms/store/user/skalafut/DoubleEG/analyzed_DoubleEG_skim_low_mass_region_eejj.root')
		#fileName = cms.string('/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_TTJets_50ns_skim_low_mass_region_eejj.root')
		fileName = cms.string('/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_WJets_skim_low_mass_region_eejj.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

inputFiles = cms.untracked.vstring()
inputFiles.extend([
	
	   #################################################################
	   #WJets to EEJJ 50ns low mass skimmed evts
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_1.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_10.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_100.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_101.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_103.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_104.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_105.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_106.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_107.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_108.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_109.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_11.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_110.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_111.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_112.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_113.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_114.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_115.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_116.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_117.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_118.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_119.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_12.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_120.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_121.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_122.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_123.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_124.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_125.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_126.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_127.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_128.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_129.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_13.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_130.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_131.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_132.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_133.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_134.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_135.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_136.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_137.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_138.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_139.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_14.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_140.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_141.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_142.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_143.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_144.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_145.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_146.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_147.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_148.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_149.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_15.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_150.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_151.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_152.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_153.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_154.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_155.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_156.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_157.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_158.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_159.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_16.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_160.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_161.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_162.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_163.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_164.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_165.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_166.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_167.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_168.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_169.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_17.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_170.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_171.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_172.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_173.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_174.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_175.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_176.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_177.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_178.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_179.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_18.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_180.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_181.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_182.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_183.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_184.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_185.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_186.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_187.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_188.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_189.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_19.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_190.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_191.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_192.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_193.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_194.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_195.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_196.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_197.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_198.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_199.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_2.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_20.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_200.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_201.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_202.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_203.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_204.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_205.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_206.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_207.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_208.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_209.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_21.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_210.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_211.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_212.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_213.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_214.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_215.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_216.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_217.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_218.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_219.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_22.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_220.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_221.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_222.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_223.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_224.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_225.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_226.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_227.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_228.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_229.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_23.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_230.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_231.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_232.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_233.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_234.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_235.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_236.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_237.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_238.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_239.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_24.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_240.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_241.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_242.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_243.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_244.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_245.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_246.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_247.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_248.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_249.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_25.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_250.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_251.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_252.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_253.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_254.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_255.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_256.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_257.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_258.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_259.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_26.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_260.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_261.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_262.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_263.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_264.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_265.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_266.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_267.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_268.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_269.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_27.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_270.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_271.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_272.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_273.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_274.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_275.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_276.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_277.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_278.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_279.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_28.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_280.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_281.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_282.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_283.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_284.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_285.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_286.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_287.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_288.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_289.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_29.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_290.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_291.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_292.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_293.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_294.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_295.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_296.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_297.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_298.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_299.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_3.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_30.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_300.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_301.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_302.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_303.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_304.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_305.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_306.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_307.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_308.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_309.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_31.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_310.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_311.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_312.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_313.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_314.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_315.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_316.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_317.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_318.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_319.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_32.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_320.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_33.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_34.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_35.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_36.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_37.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_38.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_39.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_4.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_40.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_41.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_42.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_43.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_44.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_45.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_46.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_47.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_48.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_49.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_5.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_50.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_51.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_52.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_53.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_54.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_55.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_56.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_57.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_58.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_59.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_6.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_60.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_61.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_62.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_63.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_64.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_65.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_66.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_67.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_68.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_69.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_7.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_70.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_71.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_72.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_73.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_74.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_75.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_76.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_77.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_78.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_79.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_8.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_80.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_81.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_82.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_83.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_84.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_85.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_86.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_87.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_88.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_89.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_9.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_90.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_91.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_92.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_93.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_94.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_95.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_96.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_97.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_98.root',
       'file:/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToEEJJ50nsSkim_13TeV_skim_low_mass_and_signal_regions_atFNALLPC/150807_132811/0000/miniAODEleEleSkimLowMassRegion_99.root',
   


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
	   #'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronLowMassRegionSkim_1.root'

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


