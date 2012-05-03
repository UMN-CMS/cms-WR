import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronsSKIMHeavyNu")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 100

from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/736/C87655EF-FB83-E111-A922-003048D3733E.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GR_R_50_V13::All'


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.1.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('pA'),
   ),
   outputCommands = cms.untracked.vstring("keep *")
)


#### Output needs to be created before working with PAT objects ###
#process.outB = cms.OutputModule( "PoolOutputModule",
#   fileName = cms.untracked.string("poolout.2.root"),
#   maxSize = cms.untracked.int32(3000000),
#   SelectEvents = cms.untracked.PSet(
#      SelectEvents = cms.vstring('p'),
#   ),
#   outputCommands = cms.untracked.vstring("keep *")
#)

###########################################################################
###  P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s  ###
###########################################################################
# cleaning: scraping filter
process.scrapingFilter = cms.EDFilter("FilterOutScraping",
                                                                            applyfilter = cms.untracked.bool(True),
                                                                            debugOn = cms.untracked.bool(False),
                                                                            numtrack = cms.untracked.uint32(10),
                                                                            thresh = cms.untracked.double(0.25)
                                                                            )



# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(3.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )


#------------------
#Load PAT sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#
# the HCAL Noise Filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

#
# this is to check event content
process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")

#
# pack all pat-related modules here
process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons") , cms.InputTag("patElectrons"))
process.patCandidates      = cms.Sequence( process.makePatMuons +
                                           process.makePatElectrons +
                                           process.makePatJets +
                                           process.patCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

# this does not filter event, only reduces the collection
process.ElectronsAbove20 = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("patElectrons"),
                                       cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(20) )
                                       )
#here's the real filter
process.oneElectronAbove20 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove20"),
                                       minNumber = cms.uint32(1)
                                       )
#here's the real filter
process.twoElectronsAbove20 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove20"),
                                       minNumber = cms.uint32(2)
                                       )

# this does not filter event, only reduces the collection
process.ElectronsAbove25 = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("ElectronsAbove20"),
                                       cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(25) )
                                       )
#here's the real filter
process.oneElectronAbove25 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove25"),
                                       minNumber = cms.uint32(1)
                                       )
#here's the real filter
process.twoElectronsAbove25 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove25"),
                                       minNumber = cms.uint32(2)
                                       )


# this does not filter event, only reduces the collection
process.ElectronsAbove30 = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("ElectronsAbove25"),
                                       cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(30) )
                                       )
#here's the real filter
process.oneElectronAbove30 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove30"),
                                       minNumber = cms.uint32(1)
                                       )
#here's the real filter
process.twoElectronsAbove30 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove30"),
                                       minNumber = cms.uint32(2)
                                       )

import HeavyNu.MuFilter.EleFilter_cfi
process.eleFilter20     = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter20.minElePt  = cms.double( 20. )
process.eleFilter20.minSCEt   = cms.double( 20. )

process.eleFilter20mass = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter20mass.minElePt  = cms.double( 20. )
process.eleFilter20mass.minSCEt   = cms.double( 20. )
process.eleFilter20mass.massMin   = cms.double( 50. )
process.eleFilter20mass.massMax   = cms.double( 130. )

process.eleFilter25 = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter25.minElePt  = cms.double( 25. )
process.eleFilter25.minSCEt   = cms.double( 25. )

process.eleFilter25mass = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter25mass.minElePt  = cms.double( 25. )
process.eleFilter25mass.minSCEt   = cms.double( 25. )
process.eleFilter25mass.massMin = cms.double( 50. )
process.eleFilter25mass.massMax = cms.double( 130. )

process.eleFilter30 = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter30.minElePt  = cms.double( 30. )
process.eleFilter30.minSCEt   = cms.double( 30. )

process.eleFilter30mass = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.eleFilter30mass.minElePt  = cms.double( 30. )
process.eleFilter30mass.minSCEt   = cms.double( 30. )
process.eleFilter30mass.massMin = cms.double( 50. )
process.eleFilter30mass.massMax = cms.double( 130. )





process.pA = cms.Path( 
    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *

    process.patDefaultSequence *

    # select electron candidates above 20 (no filtering)
    process.ElectronsAbove20 *
    # select electron candidates above 25 (no filtering)
    process.ElectronsAbove25 *
    # select electron candidates above 30 (no filtering)
    process.ElectronsAbove30 *

    # filter on presence of 1electron>20
    process.oneElectronAbove20  * 
    # filter on presence of 1electron>25
    process.oneElectronAbove25  * 
    # filter on presence of 1electron>30
    process.oneElectronAbove30  * 

    # filter on presence of 2electron(s)>20
    process.twoElectronsAbove20  * 
    # filter on presence of electron(s)>25
    process.twoElectronsAbove25 *
    # filter on presence of electron(s)>30
    process.twoElectronsAbove30

    # * process.dumpEvContent
)


process.pB = cms.Path( 
    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *

    process.patDefaultSequence *

    # select electron candidates above 20 (no filtering)
    process.ElectronsAbove20 *
    # select electron candidates above 25 (no filtering)
    process.ElectronsAbove25 *
    # select electron candidates above 30 (no filtering)
    process.ElectronsAbove30 *

    # filter on presence of 1electron>20
    process.oneElectronAbove20  * 

    process.eleFilter20 *
    process.eleFilter20mass *

    # filter on presence of 1electron>25
    process.oneElectronAbove25  * 

    process.eleFilter25 *
    process.eleFilter25mass *


    # filter on presence of 1electron>30
    process.oneElectronAbove30  * 

    process.eleFilter30 *
    process.eleFilter30mass 

    # * process.dumpEvContent
)


process.outpath = cms.EndPath(process.out)
#process.outpath = cms.EndPath(process.out2)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])
