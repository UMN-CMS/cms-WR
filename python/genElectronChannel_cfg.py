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

process.matchedGenAnalyzerOne = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedGenObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle","","GENEEJJ"),
		subleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle","","GENEEJJ"),
		quarkCollection = cms.InputTag("bareMatchedGenQuark","","GENEEJJ"),
		matchingLeadingLeptonCollection = cms.InputTag("","",""),
		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
		matchingQuarkCollection = cms.InputTag("","","")
		)


## run these analyzers after applying pt and eta cuts on the leptons and jets
process.genAnalyzerTwo = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

process.matchedGenAnalyzerTwo = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ"),
		matchingLeadingLeptonCollection = cms.InputTag("","",""),
		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
		matchingQuarkCollection = cms.InputTag("","","")
		)


## run these analyzers after applying pt, eta, and dilepton mass cuts
process.genAnalyzerThree = cms.EDAnalyzer('genAnalyzer',
		treeName = cms.string("genObjectsWithPtEtaAndDileptonMassCuts"),
		genLeptonCollection = cms.InputTag("ptEtaRestrictedSubleadingGenEle","","GENEEJJ"),
		genJetCollection = cms.InputTag("ptEtaRestrictedGenJet","","GENEEJJ")
		)

process.matchedGenAnalyzerThree = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaAndDileptonMassCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ"),
		matchingLeadingLeptonCollection = cms.InputTag("","",""),
		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
		matchingQuarkCollection = cms.InputTag("","","")
		)

## analyzers to run with dR(lepton,hadron) > someValue cut applied
process.matchedGenAnalyzerFour = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaDileptonMassAndDrCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ"),
		matchingLeadingLeptonCollection = cms.InputTag("","",""),
		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
		matchingQuarkCollection = cms.InputTag("","","")
		)

## analyzers to run with dR(lepton,hadron) > someValue cut and four obj inv mass cut applied
process.matchedGenAnalyzerFive = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedGenObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(True),
		minFourObjMass = cms.double(600),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		saveGenMatched = cms.bool(False),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle","","GENEEJJ"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle","","GENEEJJ"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedGenQuark","","GENEEJJ"),
		matchingLeadingLeptonCollection = cms.InputTag("","",""),
		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
		matchingQuarkCollection = cms.InputTag("","","")
		)



process.printMuOrTauDecayChain = cms.Path(
		#process.hasGenNuMuOrTau
		#*process.hasGenNuMuOrTauFilter
		process.printParticleTree
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

process.schedule = cms.Schedule(process.runMatchedGenAnalysis)

process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('analysis_genElectronChannel_using_miniAOD.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_1300_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_520_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_2080_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_1400_MNu_700_using_GEN_SIM.root')
		fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_50_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_100_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_200_using_GEN_SIM.root')
		#fileName = cms.string('analysis_genElectronChannel_MWR_2600_MNu_300_using_GEN_SIM.root')





)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)


process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		##WR->e + Nu_e --> eejj miniAOD files 
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

		##WR->e + Nu_e ->eejj GEN-SIM files
		##MWR 2600 MNu 1300
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_2.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_27.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_28.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_29.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_30.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_31.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_32.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_33.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_34.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_35.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_36.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_37.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_38.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_39.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_40.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_41.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_42.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_43.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_44.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_45.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_46.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_47.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_48.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_49.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_5.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_50.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_51.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_52.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_53.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_54.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_55.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_56.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_57.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_58.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_59.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_60.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_61.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_62.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_63.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_64.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_65.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_66.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_67.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_68.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_69.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_7.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_70.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_71.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_72.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_73.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_74.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_75.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_76.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_77.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_78.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_79.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_8.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_80.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_1300_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-1300_13TeV_GEN-SIM/150519_203237/0000/WR_M-2600_ToENu_M-1300_13TeV_GEN_SIM_9.root',

		##MWR 2600 MNu 520
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_2.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_27.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_28.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_29.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_30.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_31.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_32.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_33.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_34.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_35.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_36.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_37.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_38.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_39.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_40.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_41.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_42.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_43.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_44.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_45.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_46.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_47.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_48.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_49.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_5.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_50.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_51.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_52.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_53.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_54.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_55.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_56.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_57.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_58.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_59.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_60.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_61.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_62.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_63.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_64.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_65.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_66.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_67.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_68.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_69.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_7.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_70.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_71.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_72.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_73.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_74.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_75.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_76.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_77.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_78.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_79.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_8.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_80.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000/WR_M-2600_ToENu_M-520_13TeV_GEN_SIM_9.root',

		##MWR 2600 MNu 2080
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_2.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_27.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_28.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_29.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_30.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_31.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_32.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_33.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_34.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_35.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_36.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_37.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_38.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_39.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_40.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_41.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_42.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_43.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_44.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_45.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_46.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_47.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_48.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_49.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_5.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_50.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_51.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_52.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_53.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_54.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_55.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_56.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_57.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_58.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_59.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_60.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_61.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_62.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_63.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_64.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_65.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_66.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_67.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_68.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_69.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_7.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_70.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_71.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_72.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_73.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_74.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_75.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_76.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_77.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_78.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_79.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_8.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_80.root',
		#'file:/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000/WR_M-2600_ToENu_M-2080_13TeV_GEN_SIM_9.root',

		##MWR 1400 MNu 700
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_2.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_27.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_28.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_29.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_30.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_31.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_32.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_33.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_34.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_35.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_36.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_37.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_38.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_39.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_40.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_41.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_42.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_43.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_44.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_45.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_46.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_47.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_48.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_49.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_5.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_50.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_51.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_52.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_53.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_54.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_55.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_56.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_57.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_58.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_59.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_60.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_61.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_62.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_63.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_64.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_65.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_66.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_67.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_68.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_69.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_7.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_70.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_71.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_72.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_73.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_74.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_75.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_76.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_77.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_78.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_79.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_8.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_80.root',
		#'file:/eos/uscms/store/user/skalafut/WR_1400_ToENu_700_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-1400_MNu-700_13TeV_GEN-SIM/150519_203412/0000/WR_M-1400_ToENu_M-700_13TeV_GEN_SIM_9.root',

		##MWR 2600 MNu 50


		##MWR 2600 MNu 100


		##MWR 2600 MNu 200


		##MWR 2600 MNu 300


		##WR->LNu->LLJJ GEN-SIM files
		#'file:/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/WR_GEN_SIM.root'

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



