import FWCore.ParameterSet.Config as cms

# Selection of electrons from the slimmedElectrons collection (PAT)
""" \addtogroup electronSkim_Group electronSkim sequences
@{
"""

### select leading electron \ingroup electronSkim_Group
wRleadingElectronPresel = cms.EDFilter("PATElectronRefSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>45"),
                                 )

### select subleading electron
wRsubleadingElectronPresel = cms.EDFilter("PATElectronRefSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>30"),
                                 )
#wRpreSelectedElectrons = cms.EDProducer("CandViewMerger",

### create di-electron pair in signal region
wRdiElectronCandidatePresel = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectronPresel wRsubleadingElectronPresel"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                        # the cut on the pt of the daughter is to respect the order of leading and subleading:
                        # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                       cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
                                   )

### filter: at least one di-electron candidate in signal region
wRdiElectronCandidateFilterPresel = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiElectronCandidatePresel"),
                                           minNumber = cms.uint32(1)
                                           )

############################################################

### @}
