import FWCore.ParameterSet.Config as cms

process = cms.Process("triggerEfficienciesEle")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

from PhysicsTools.PatAlgos.tools.metTools import *
#from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource",      
                            fileNames=cms.untracked.vstring(
    #'file:/data/franzoni/data/data-Run2012A-DoubleElectron/Run2012A-DoubleElectron-AOD-PromptReco-v1-000-191-833-36C97ED9-978C-E111-B635-001D09F23F2A.root',
    #'file:/data/franzoni/data/data-Run2012A-DoubleElectron/Run2012A-DoubleElectron-AOD-PromptReco-v1-000-191-833-5ADFF104-9C8C-E111-A37B-003048D3C982.root',
    #'file:/data/franzoni/data/data-Run2012A-DoubleElectron/Run2012A-DoubleElectron-AOD-PromptReco-v1-000-191-833-7CAB7DB8-968C-E111-B3B3-BCAEC532971D.root',
    'file:/data/franzoni/data/data-Run2012A-DoubleElectron/data-Run2012A-DoubleElectron-Run2012A-DoubleElectron-AOD-PromptReco-v1-000-193-621-84C87551-9B9A-E111-AF22-5404A63886A5.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag = 'GR_R_50_V13::All'


### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.1.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('pA'),
   ),
   outputCommands = cms.untracked.vstring("keep *")
)



###########################################################################
###  P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s  ###
###########################################################################
#------------------
#Load PAT sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.load("PhysicsTools.PatAlgos.tools.pfTools")
from PhysicsTools.PatAlgos.tools.jetTools import *
#process.patDefaultSequence.remove( process.patTaus )    

isMC=False
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                                jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))


process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTrigger#PATTriggerProducer
process.patTrigger.addL1Algos       = cms.bool(True)
process.patTriggerEvent.condGtTag   = cms.InputTag( 'conditionsInEdm' )

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

#
# this is to check event content
process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")


#
# pack all pat-related modules here
process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons") , cms.InputTag("patElectrons"))

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from RecoJets.JetProducers.kt4PFJets_cfi import *

# compute FastJet rho to correct jets
process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets","rho")


# compute FastJet rho to correct isolation
process.kt6PFJetsForIsolation = kt4PFJets.clone(
    rParam = 0.6,
    doRhoFastjet = True
    )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


inputJetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
# Add PF jets
switchJetCollection(process,
                    cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = inputJetCorrLabel,
                    doType1MET   = True,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True
                    )

process.patJets.addTagInfos = True
process.patJets.tagInfoSources  = cms.VInputTag( cms.InputTag("secondaryVertexTagInfosAOD") )



from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
process.HLTTagAndProbe = cms.EDFilter("TriggerResultsFilter",
                                      hltResults              = cms.InputTag('TriggerResults','','HLT'),   # HLT results   - set to empty to ignore HLT
                                      l1tResults              = cms.InputTag(''),                 # L1 GT results - set to empty to ignore L1
                                      l1tIgnoreMask           = cms.bool(False),                  # use L1 mask
                                      l1techIgnorePrescales   = cms.bool(False),                  # read L1 technical bits from PSB#9, bypassing the prescales
                                      daqPartitions           = cms.uint32(0x01),                 # used by the definition of the L1 mask
                                      throw                   = cms.bool(False),                  # if HLT path not in the table, crash/ignore according to true/false
                                      triggerConditions       = cms.vstring(
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5'
    )
                                      )


#
#  the analysis comes in here
#
import HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi

process.triggerEfficiencies                     = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
process.triggerEfficiencies.plotFolderName      = cms.string('inclusive')
process.triggerEfficiencies.numOfflJets         = cms.int32(  0 )
    
process.triggerEfficienciesOneJ                 = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
process.triggerEfficienciesOneJ.plotFolderName  = cms.string('OneJet')
process.triggerEfficienciesOneJ.numOfflJets     = cms.int32(  1 )
    
process.triggerEfficienciesTwoJ                 = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
process.triggerEfficienciesTwoJ.plotFolderName  = cms.string('TwoJet')
process.triggerEfficienciesTwoJ.numOfflJets     = cms.int32(  2 )
    
#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("triggereffele.plots.root"),
                                   )

# this is the path for the L1 study 
process.pA = cms.Path( 

    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *
    
    # if seeding HLT path not passed, ignore event
    process.HLTTagAndProbe *
    
    # re-compute rho in order to have cone at 0.6 and be as close as possible to heavyNuAnalysis_cfg.py
    process.kt6PFJets *
    process.kt6PFJetsForIsolation *
    
    process.patDefaultSequence *
    
    # process.dumpEvContent *

    # cerate triggere objects and triggerEvent
    # documentation here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTrigger#PATTriggerProducer
     process.patTriggerDefaultSequence *
        
    process.triggerEfficiencies *
    process.triggerEfficienciesOneJ *
    process.triggerEfficienciesTwoJ
    
    )

#process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])
