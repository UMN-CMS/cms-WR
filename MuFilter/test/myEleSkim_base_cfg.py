import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronsSKIMHeavyNu")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

#from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring('file:input.root')
                            #fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/736/C87655EF-FB83-E111-A922-003048D3733E.root')
                            fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/data/Run2012B/DoublePhotonHighPt/AOD/13Jul2012-v1/00000/00FD2120-45D9-E111-863D-848F69FD28AA.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

#process.GlobalTag.globaltag = 'GR_R_50_V13::All'
process.GlobalTag.globaltag = 'FT_53_V6_AN1::All'


### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("myEleSkim_base.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('pA'),
   ),
   outputCommands = cms.untracked.vstring("keep *")
)


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



###########################################################################
###  P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s  ###
###########################################################################
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
process.patCandidates      = cms.Sequence(
    #process.makePatMuons +
    process.makePatElectrons +
    #process.makePatJets +
    process.patCandidateSummary
    )
process.patDefaultSequence = cms.Sequence( process.patCandidates )


# this does not filter event, only reduces the collection
process.ElectronsAbove35 = cms.EDFilter("CandViewSelector",
                                        src = cms.InputTag("patElectrons"),
                                        cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(35) )
                                        )
#here's the real filter
process.twoElectronsAbove35 = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("ElectronsAbove35"),
                                           minNumber = cms.uint32(2)
                                           )


process.pA = cms.Path( 
    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *
    
    process.patDefaultSequence *
    
    # select electron candidates above 35 (no filtering)
    process.ElectronsAbove35 *
    
    # filter on presence of electron(s)>35
    process.twoElectronsAbove35
    
    # * process.dumpEvContent
    )


process.outpath = cms.EndPath(process.out)


from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])
