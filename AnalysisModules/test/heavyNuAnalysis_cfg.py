import FWCore.ParameterSet.Config as cms

#isMC=False
isMC=True

process = cms.Process("PAT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
#process.source = cms.Source("PoolSource",      
##    fileNames=cms.untracked.vstring( "file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_1-reco-pool.root" )
#    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/S10/TTbar/GEN-SIM-RECO/START36_V9_S09-v1/0055/044B121D-9678-DF11-B39B-001E0B5F95B2.root')
#)
process.load('HeavyNu.AnalysisModules.in_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC FLAG IS SET <===================="
    process.GlobalTag.globaltag = cms.string('START36_V10::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_36X_V12::All')

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
       minMu1pt     = cms.double(60.),
       minMu2pt     = cms.double(15.),
       minJetPt     = cms.double(40),
       minMuMuMass  = cms.double(140),
       minMuonJetdR = cms.double(0.3), # (0.5)
)

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
