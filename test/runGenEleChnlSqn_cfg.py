import FWCore.ParameterSet.Config as cms

process = cms.Process("GENEEJJ")

## load the filters, producers, and sequences defined in the config file fragment
process.load('ExoAnalysis.cmsWR.linkGenSequencesForTwoDimLimits_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()

#process.matchedGenAnalyzerOne = cms.EDAnalyzer('matchedAnalyzer',
#		treeName = cms.string("matchedGenObjectsNoCuts"),
#		doDeltaRcut = cms.bool(False),
#		doFourObjMassCut = cms.bool(False),
#		minFourObjMass = cms.double(-1),
#		minDeltaRforLeptonJetExclusion = cms.double(-1),
#		saveGenMatched = cms.bool(False),
#		leadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle","","GENEEJJ"),
#		subleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle","","GENEEJJ"),
#		quarkCollection = cms.InputTag("bareMatchedGenQuark","","GENEEJJ"),
#		matchingLeadingLeptonCollection = cms.InputTag("","",""),
#		matchingSubleadingLeptonCollection = cms.InputTag("","",""),
#		matchingQuarkCollection = cms.InputTag("","","")
#		)



process.runMatchedGenAnalysis = cms.Path(
		process.eeChnlCutSeq
		)

process.schedule = cms.Schedule(process.runMatchedGenAnalysis)

process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('analysis_genElectronChannel_MWR_800_MNu_400_using_miniAOD.root')
		fileName = cms.string(options.output)

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(options.files),
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)



