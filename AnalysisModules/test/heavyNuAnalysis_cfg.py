import FWCore.ParameterSet.Config as cms
import sys

#isMC=sys.modules['__main__'].isMC
#isSIGNAL=sys.modules['__main__'].isSIGNAL
#process = sys.modules['__main__'].process

Training=False
isMC=True
isSIGNAL=False

process = cms.Process("PAT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
#process.source = cms.Source("PoolSource",
#    fileNames=cms.untracked.vstring('file:/local/cms/user/dahmes/wr2010/run2010B/oct7/mySkim/mySkim-oct7_005001.root')
#    fileNames=cms.untracked.vstring('/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/148/829/9E1A7E3C-89E0-DF11-9CDF-000423D98B6C.root')
#    fileNames=cms.untracked.vstring( "file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_1-reco-pool.root" )
##    fileNames=cms.untracked.vstring('/store/mc/S10/TTbar/GEN-SIM-RECO/START36_V9_S09-v1/0055/044B121D-9678-DF11-B39B-001E0B5F95B2.root')
#)
#process.load('HeavyNu.AnalysisModules.in_cff')

#if not isMC:
#    from HeavyNu.AnalysisModules.goodRunList_cfi import lumisToProcess
#    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC flag is SET <===================="
#    process.GlobalTag.globaltag = cms.string('START36_V10::All')
    process.GlobalTag.globaltag = cms.string('START38_V12::All')
else:
    print "=================> MC flag is NOT SET <===================="
    process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
##
if not isMC:
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'])

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)
process.hNu = cms.EDFilter("HeavyNu",
                           DoLog        = cms.bool( True ),
                           BtagName     = cms.string('jetProbabilityBJetTags'),
                           minMu1pt     = cms.double(60.),
                           minMu2pt     = cms.double(15.),
                           minJetPt     = cms.double(40),
                           maxJetAbsEta = cms.double(2.5),
                           minMuMuMass  = cms.double(200),
                           min4objMass  = cms.double(520),
                           minMuonJetdR = cms.double(0.3), # (0.5)
                           muonTrackIsoLimitGeV = cms.double(10.0),

                           isSignal     = cms.bool(False),
                           mNuRnormalization = cms.double(1000.0)
                           )

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
if isSIGNAL:
    process.hNu.isSignal = cms.bool(True)

process.myAnalysisPath = cms.Path(
    process.patDefaultSequence *
    process.hNu 
)

if not isMC:
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('heavynu_candevents.root'),
    # save only events passing the full path
                                      SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('myAnalysisPath')),
                                  
                                      #outputCommands = cms.untracked.vstring("keep *")
                                      )
    process.out_step = cms.EndPath(process.output)
    process.schedule = cms.Schedule(process.myAnalysisPath,
                                    process.out_step)
