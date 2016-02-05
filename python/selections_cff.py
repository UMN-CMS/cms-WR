import FWCore.ParameterSet.Config as cms
### \todo add deltaR cut between leptons?
### \todo add HEEP ID for electrons
### \todo add muon ID
### \todo check composite candidate cut string for pt ordering: pt_0 >= pt_1
### \todo make sure that the miniTree is using objects with looser selection


###

leadingPt=60.
subleadingPt=50.
miniTreeLeptonPt=40.
maxEtaLeptons=2.4
jetPt = 20. # same as preselection


### need the reference for this selection, with link to the presentation and twiki page
jetID=" (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4)"
#muonIDIso=" isolationR03().sumPt/pt < 0.1"
muonIDIso=' isolationR03().sumPt/pt < 0.1 && userInt("highPtID") == 1'



############################################################ Jets
from ExoAnalysis.cmsWR.JEC_cff import *

#patJetsReapplyJEC

### select tight-ID jets
wRtightJets = cms.EDFilter("PATJetSelector",
                          src = cms.InputTag("patJetsReapplyJEC"),
                          cut = cms.string(jetID),
                          )

wRJets = cms.EDFilter("PATJetSelector",
                      src = cms.InputTag("wRtightJets"),
                      cut = cms.string( ("pt>%f") % (jetPt)),
		)

wRJECUncert = JECUnc.clone(
    src = cms.InputTag('wRJets')
    )

wRJetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("wRJets"),
                                 minNumber = cms.uint32(2)
                             )

wRdiJetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRJets wRJets"),
                                  #role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
)

jetSelectionSeq = cms.Sequence( wRtightJets * wRJets * wRJECUncert ) #* wRdiJetCandidate)

############################################################ Electrons


# muons considered here are only those far from any jet with Pt>ptJets 
wRIsolatedElectrons = cms.EDFilter( "DeltaROverlapExclusionSelector",
                                   src = cms.InputTag("slimmedElectrons"),
                                   overlap = cms.InputTag("wRJets"),
                                   maxDeltaR = cms.double(0.4),
                                   )

#wRHEEPElectron = 
wRHEEPElectron = cms.EDProducer('HEEPIDSelector',
                                electrons= cms.InputTag("slimmedElectrons"),
                                eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                                )

#### this is to change the type of the collection for the analyzer
wRHEEPElectronRefiner = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("wRHEEPElectron"),
		cut = cms.string('')
		)

### select leading electron \ingroup electronSkim_Group
wRleadingElectron = cms.EDFilter("PATElectronSelector",
                                 src = cms.InputTag("wRHEEPElectron"), #wRIsolatedElectrons"),
                                 cut = cms.string( 
        (("(pt>%f) && (abs(eta)<%f)") % (leadingPt, maxEtaLeptons))
        ),
                                 )

### select subleading electron
wRsubleadingElectron = cms.EDFilter("PATElectronSelector",
                                    src = cms.InputTag("wRHEEPElectron"), #wRIsolatedElectrons"),
                                    cut = cms.string( 
        (("(pt>%f) && (abs(eta)<%f)") % (subleadingPt, maxEtaLeptons))
        ),
                                    )

wRminiTreeElectron = cms.EDFilter("PATElectronSelector",
                                    src = cms.InputTag("wRHEEPElectron"), 
                                    cut = cms.string( 
        (("(pt>%f) && (abs(eta)<%f)") % (miniTreeLeptonPt, maxEtaLeptons))
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


electronSelectionSeq = cms.Sequence(wRIsolatedElectrons *  wRHEEPElectron * wRHEEPElectronRefiner * (wRleadingElectron + wRsubleadingElectron + wRminiTreeElectron) * wRdiElectronCandidate)

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

### add here the rochester corrections

### muon ID and isolation
tunePIDIsoMuons = cms.EDFilter("PATMuonSelector",
                               src = cms.InputTag("tunePMuons"),
                               cut = cms.string(muonIDIso),
                               )

# make a collection of TuneP muons which pass isHighPt ID
scaleCorrectedMuons = cms.EDProducer("produceScaleCorrectedMuons",
		src = cms.InputTag("tunePIDIsoMuons"),
                OutputCollectionName = cms.string("")
		)

scaleCorrectedMuonsRefiner = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("scaleCorrectedMuons"),
		cut = cms.string('')
		)

### add here the trigger matching

### select leading muon \ingroup muonSkim_Group
wRleadingMuon = cms.EDFilter("PATMuonSelector",
                             src = cms.InputTag("scaleCorrectedMuons"), #tunePIsolatedMuons"),
                             cut = wRleadingElectron.cut,
                         )

### select subleading muon
wRsubleadingMuon = cms.EDFilter("PATMuonSelector",
                                src = cms.InputTag("scaleCorrectedMuons"), #tunePIsolatedMuons"),
                                cut = wRsubleadingElectron.cut,
                                )

wRminiTreeMuon = cms.EDFilter("PATMuonSelector",
                                src = cms.InputTag("scaleCorrectedMuons"), #tunePIsolatedMuons"),
                                cut = wRminiTreeElectron.cut,
                                )

from ExoAnalysis.cmsWR.produceIdIsoSF_cff import *


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

muonSelectionSeq = cms.Sequence(tunePMuons * tunePIDIsoMuons * scaleCorrectedMuons * scaleCorrectedMuonsRefiner  * wRleadingMuon * wRsubleadingMuon * muonIdIsoSF * wRdiMuonCandidate)
#muonSelectionSeq = cms.Sequence(tunePMuons * tunePIDIsoMuons * scaleCorrectedMuons * wRleadingMuon * wRsubleadingMuon * muonIdIsoSF * wRdiMuonCandidate)
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

wRdiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiLeptonCandidate"),
                                           minNumber = cms.uint32(1)
                                           )

wRCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRdiLeptonCandidate wRdiJetCandidate"),
                             # role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0"),
)
############################################################# Flavour sideband filter
flavourSidebandSelection='(daughter(0).isElectron && daughter(1).isMuon) || (daughter(1).isElectron && daughter(0).isMuon)'
diLeptonSelection = '(mass>200)'


lowFourObjectSidebandSelector = cms.EDFilter("CandViewSelector",
                                     src = cms.InputTag("wRCandidate"),
                                     cut = cms.string('(mass< 600)')
                                     )
lowFourObjectSidebandFilter  = cms.EDFilter('CandViewCountFilter',
                                          src = cms.InputTag('lowFourObjectSidebandSelector'),
                                          minNumber = cms.uint32(1)
)


flavourSidebandSelector = cms.EDFilter("CandViewSelector",
                                     src = cms.InputTag("wRdiLeptonCandidate"),
                                     cut = cms.string(flavourSidebandSelection+" && " + diLeptonSelection),
                                     )

flavourSidebandFilter = cms.EDFilter('CandViewCountFilter',
                                     src = cms.InputTag('flavourSidebandSelector'),
                                     minNumber = cms.uint32(1)
)

## to be fixed
lowDiLeptonSidebandSelector = cms.EDFilter("CandViewSelector",
                                     src = cms.InputTag("wRdiLeptonCandidate"),
                                     cut = cms.string('(mass< 180) && (!'+flavourSidebandSelection+')' )
                                     )
lowDiLeptonSidebandFilter  = cms.EDFilter('CandViewCountFilter',
                                          src = cms.InputTag('lowDiLeptonSidebandSelector'),
                                          minNumber = cms.uint32(1)
)

####
signalRegionSelector = cms.EDFilter("CandViewSelector",
                                     src = cms.InputTag("wRdiLeptonCandidate"),
                                     cut = cms.string("(! "+flavourSidebandSelection+') && '+diLeptonSelection),
                                     )

signalRegionFilter =  cms.EDFilter('CandViewCountFilter',
                                     src = cms.InputTag('signalRegionSelector'),
                                     minNumber = cms.uint32(1),
)

signalRegionEESelector = cms.EDFilter("CandViewSelector",
                                      src = cms.InputTag('signalRegionSelector')
                                      cut = cms.string('(daughter(0).isElectron && daughter(1).isElectron)'),
                                      )
signalRegionMuMuSelector = cms.EDFilter("CandViewSelector",
                                      src = cms.InputTag('signalRegionSelector')
                                      cut = cms.string('(daughter(0).isMuon && daughter(1).isMuon)'),
                                      )
signalRegionEEFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("signalRegionEESelector"),
                                    minNumber = cms.uint32(1),
                                    )
signalRegionMuMuFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("signalRegionMuMuSelector"),
                                    minNumber = cms.uint32(1),
                                    )

signalRegionFilterSeq = cms.Sequence(~flavourSidebandFilter * ~lowDiLeptonSidebandSelector)



selectionSequence = cms.Sequence(
    ( jetSelectionSeq + electronSelectionSeq + muonSelectionSeq ) * wRIsolatedElectrons * wReleMuCandidate * wRmuEleCandidate * wRdiLeptonCandidate *  flavourSidebandSelector * lowDiLeptonSidebandSelector * signalRegionSelector * signalRegionEESelector * signalRegionMuMuSelector
        )

filterSequence = cms.Sequence(
        wRJetFilter * wRdiLeptonCandidateFilter
        )
