import FWCore.ParameterSet.Config as cms

# Selection of muons from the slimmedMuons collection (PAT)
""" \addtogroup muonSkim_Group muonSkim sequences
@{
"""

tunePMuons = cms.EDProducer("TunePMuonProducer",
		src = cms.InputTag("slimmedMuons")
		)


### select leading muon \ingroup muonSkim_Group
wRleadingMuonPresel = cms.EDFilter("PATMuonRefSelector",
                             src = cms.InputTag("tunePMuons"),
                             cut = cms.string("pt>45"),
                         )

### select subleading muon
wRsubleadingMuonPresel = cms.EDFilter("PATMuonRefSelector",
                                src = cms.InputTag("tunePMuons"),
                                cut = cms.string("pt>30"),
                            )


### create di-muon pair in signal region
wRdiMuonCandidatePresel = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("wRleadingMuonPresel wRsubleadingMuonPresel"),
                                   role = cms.string("leading subleading"),
                                   checkCharge = cms.bool(False),
                                   # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                   # if both muons have pt>60 GeV there will be two di-muon candidates with inversed order
                                   #cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                   cut = cms.string("mass > 0  && daughter(0).pt>daughter(1).pt"),
                                   
									   )

### filter: at least one di-muon candidate in signal region
wRdiMuonCandidateFilterPresel = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiMuonCandidatePresel"),
                                           minNumber = cms.uint32(1)
                                           )

############################################################




### create an object from two jets and two muons in the evt, and cut on its mass
wRdiMuonDijetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("muwRdiJetCandidate wRdiMuonCandidate"),
		role = cms.string("dijet dilepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 500")
		)

### filter: require at least one LLJJ object in the evt
wRdiMuonDijetCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("wRdiMuonDijetCandidate"),
		minNumber = cms.uint32(1)
		)



### create di-muon pair in sideband region
wRdiMuonSidebandCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingMuonPresel wRsubleadingMuonPresel"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("0< mass < 200 && daughter(0).pt>daughter(1).pt")
                                       )


### @}
