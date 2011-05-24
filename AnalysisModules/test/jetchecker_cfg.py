import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring("/store/data/Run2011A/DoubleElectron/RECO/May7ReReco-v1/0001/8A8D04CE-7C7A-E011-9C3B-00261894393C.root")
    
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo-alex_r0x.root')
    )

process.demo = cms.EDAnalyzer('JetChecker',
    maxDeltaVRs = cms.vdouble(5,6,7)
)


process.p = cms.Path(process.demo)
