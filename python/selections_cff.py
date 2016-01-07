import FWCore.ParameterSet.Config as cms
### \todo add deltaR cut between leptons?
### \todo add HEEP ID for electrons
### \todo add muon ID


leadingPt=60.
subleadingPt=40.
maxEtaLeptons=2.4
jetPt = 40.

### need the reference for this selection, with link to the presentation and twiki page
jetID=" (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4)"



############################################################ Jets
from ExoAnalysis.cmsWR.JEC_cff import *

#patJetsReapplyJEC
wRJets = cms.EDFilter("PATJetRefSelector",
                      src = cms.InputTag("patJetsReapplyJEC"),
                      cut = cms.string( ("pt>%f") % (jetPt)),
		)

### select loose-ID jets
wRlooseJet = cms.EDFilter("PATJetRefSelector",
                          src = cms.InputTag("wRJets"),
                          cut = cms.string(jetID),
                          )

wRJetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("wRJets"),
                                 minNumber = cms.uint32(2)
                             )


############################################################ Electrons


# muons considered here are only those far from any jet with Pt>ptJets 
wRIsolatedElectrons = cms.EDFilter( "DeltaROverlapExclusionSelector",
                                   src = cms.InputTag("slimmedElectrons"),
                                   overlap = cms.InputTag("wRJets"),
                                   maxDeltaR = cms.double(0.4),
                                   )

### select leading electron \ingroup electronSkim_Group
wRleadingElectron = cms.EDFilter("PATElectronRefSelector",
                                 src = cms.InputTag("wRIsolatedElectrons"),
                                 cut = cms.string( 
        (("pt>%f") % (leadingPt)) 
        + (("abs(eta)<%f") % (maxEtaLeptons))
        ),
                                 )

### select subleading electron
wRsubleadingElectron = cms.EDFilter("PATElectronRefSelector",
                                    src = cms.InputTag("wRIsolatedElectrons"),
                                    cut = cms.string( 
        (("pt>%f") % (subleadingPt)) +
        (("abs(eta)<%f") % (maxEtaLeptons))
        ),
                                    )

### create di-electron pair in signal region
wRdiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectron wRsubleadingElectron"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                        # the cut on the pt of the daughter is to respect the order of leading and subleading:
                        # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                       cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
                                   )

### filter: at least one di-electron candidate in signal region
wRdiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiElectronCandidate"),
                                           minNumber = cms.uint32(1)
                                           )


############################################################ Muons

tunePMuons = cms.EDProducer("TunePMuonProducer",
		src = cms.InputTag("slimmedMuons")
		)

# muons considered here are only those far from any jet with Pt>ptJets 
tunePIsolatedMuons = cms.EDFilter( "DeltaROverlapExclusionSelector",
                                   src = cms.InputTag("tunePMuons"),
                                   overlap = cms.InputTag("wRJets"),
                                   maxDeltaR = cms.double(0.04),
                                   )



### select leading muon \ingroup muonSkim_Group
wRleadingMuon = cms.EDFilter("PATMuonRefSelector",
                             src = cms.InputTag("tunePIsolatedMuons"),
                             cut = wRleadingElectron.cut,
                         )

### select subleading muon
wRsubleadingMuon = cms.EDFilter("PATMuonRefSelector",
                                src = cms.InputTag("tunePIsolatedMuons"),
                                cut = wRsubleadingElectron.cut,
                                )


### create di-muon pair in signal region
wRdiMuonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("wRleadingMuon wRsubleadingMuon"),
                                   role = cms.string("leading subleading"),
                                   checkCharge = cms.bool(False),
                                   # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                   # if both muons have pt>60 GeV there will be two di-muon candidates with inversed order
                                   #cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                   cut = cms.string("mass > 0  && daughter(0).pt>daughter(1).pt"),
                                   
									   )

### filter: at least one di-muon candidate in signal region
wRdiMuonCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiMuonCandidate"),
                                           minNumber = cms.uint32(1)
                                           )


############################################################ E-Mu candidate
#mixed flavour candidates
wReleMuCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRleadingElectron wRsubleadingMuon"),
                                  role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
)

wRmuEleCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRleadingMuon wRsubleadingElectron"),
                                  role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
)


############################################################# di lepton candidate
wRdiLeptonCandidate = cms.EDProducer("CandViewMerger",
                                     src = cms.VInputTag("wRdiElectronCandidate", "wRdiMuonCandidate", "wReleMuCandidate", "wRmuEleCandidate")
                                     )


############################################################# Flavour sideband filter
flavourSidebandFilter = cms.EDFilter("CompositeCandAndSelector",
                                     src = cms.InputTag("wRdiLeptonCandidate"),
                                     cut = cms.string("daughter(0).isElectron && daughter(1).isElectron"),
                                     )

lowDiLeptonSidebandFilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("wRdiElectronCandidate"),
                                     minNumber = cms.uint32(1),
                                     )

selectionSequence = cms.Sequence(
    wRJets
    )
                                     
