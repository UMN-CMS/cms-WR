import FWCore.ParameterSet.Config as cms

process = cms.Process("GenAndRecoDYJetsHTandInclusive")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.genDYJetsElectronChannelModules_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('ExoAnalysis.cmsWR.JEC_cff')
process.load('ExoAnalysis.cmsWR.recoElectronChannelSidebandUnmatchedModules_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v4', '')

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")


#################################
#Filters
process.genHTFilterOne = cms.EDFilter('htFilter',
		cutThreshold = cms.double(100.),
		isLowerBound = cms.bool(False),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

process.genHTFilterTwo = cms.EDFilter('htFilter',
		cutThreshold = cms.double(200.),
		isLowerBound = cms.bool(False),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)


#################################
#Analyzers

## analyze the kinematic distributions of electrons and jets at gen level with these analyzers
process.genMatchedAnalyzerOne = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("genGluonsAndQuarksNoCuts"),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

process.genMatchedAnalyzerTwo = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("genElectronsNoCuts"),
		inputCollection = cms.InputTag("dyJetsBareMatchedGenEle")
		)

process.recoAnalyzerOne = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("recoJetsNoCuts"),
		inputCollection = cms.InputTag("bareRecoJet")
		)

process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoJetsAndLeptonsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		jetsCollection = cms.InputTag("bareRecoJet"),
		leptonsCollection = cms.InputTag("zeeCheckLepton")
		)


#################################
#Paths
process.studyDYJetsDecay = cms.Path(
		#process.printDYJetsParticleTree
		process.dyJetsBareMatchedGenGluon
		*process.dyJetsBareMatchedGenQuark
		*process.dyJetsMergeGenMatchedPartons #merge gen quark and gen gluon collections
		*process.genHTFilterOne
		#*process.genHTFilterTwo
		##select electrons from Z then calculate Z pT
		#*process.dyJetsBareMatchedGenEle
		#*process.dyJetsBareMatchedGenEleFilter
		#*process.genMatchedAnalyzerOne
		#*process.genMatchedAnalyzerTwo
		###RECO sumPt of two leading jets
		#process.jecSequence
		#*process.patJetCorrFactorsReapplyJEC
		#+process.patJetsReapplyJEC
		#*process.bareRecoJet
		#*process.bareRecoJetFilter
		#*process.recoAnalyzerOne
		#+process.egmGsfElectronIDSequence
		#*process.HEEPIDSequence
		#*process.checkZeeSeq
		#*process.recoAnalyzerTwo
		)


process.schedule = cms.Schedule(process.studyDYJetsDecay)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('recoHadronAndLeptonKinematicsGenDYJetsMadHT600toInf.root')
		#fileName = cms.string('recoHadronAndLeptonKinematicsGenDYJetsMadHT400to600.root')
		#fileName = cms.string('recoHadronAndLeptonKinematicsGenDYJetsMadHT200to400.root')
		#fileName = cms.string('recoHadronAndLeptonKinematicsGenDYJetsMadHT100to200.root')
		#fileName = cms.string('genHadronKinematicsGenDYJetsMadIncl_HT100to200.root')
		fileName = cms.string('discard.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT600toInfMiniAOD0.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT600toInfMiniAOD1.root'
	#	),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT400to600MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT400to600MiniAOD0.root'
	#	),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT200to400MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT200to400MiniAOD0.root'
	#	),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT100to200MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT100to200MiniAOD2.root'
	#	),
	fileNames = cms.untracked.vstring(
		#'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadInclusiveMiniAOD1.root',
		#'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadInclusiveMiniAOD2.root',
		#'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadInclusiveMiniAOD3.root',
		'root://xrootd-cms.infn.it//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/04A311B8-5D6D-E511-A052-00266CFCC490.root'
		),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


