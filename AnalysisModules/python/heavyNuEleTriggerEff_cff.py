import FWCore.ParameterSet.Config as cms

import HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi as hNuEleTriggerEff


dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")




hNuEtriggerEff                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEff.numOfflJets         = cms.int32(  0 )

hNuEtriggerEffOneJ                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJ.numOfflJets     = cms.int32(  1 )

hNuEtriggerEffTwoJ                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJ.numOfflJets     = cms.int32(  2 )

hNuEtriggerEffnoHEEP                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffnoHEEP.numOfflJets         = cms.int32(  0 )
hNuEtriggerEffnoHEEP.heepVersion         = cms.int32(  0 )

hNuEtriggerEffOneJnoHEEP                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJnoHEEP.numOfflJets     = cms.int32(  1 )
hNuEtriggerEffOneJnoHEEP.heepVersion     = cms.int32(  0 )

hNuEtriggerEffTwoJnoHEEP                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJnoHEEP.numOfflJets     = cms.int32(  2 )
hNuEtriggerEffTwoJnoHEEP.heepVersion     = cms.int32(  0 )



hNuEtriggerEffBef                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffBef.numOfflJets         = cms.int32(  0 )
hNuEtriggerEffBef.runMax              = cms.int32(   191717)

hNuEtriggerEffOneJBef                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJBef.numOfflJets     = cms.int32(  1 )
hNuEtriggerEffOneJBef.runMax          = cms.int32(   191717)


hNuEtriggerEffTwoJBef                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJBef.numOfflJets     = cms.int32(  2 )
hNuEtriggerEffTwoJBef.runMax          = cms.int32(   191717)


hNuEtriggerEffnoHEEPBef                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffnoHEEPBef.numOfflJets         = cms.int32(  0 )
hNuEtriggerEffnoHEEPBef.heepVersion         = cms.int32(  0 )
hNuEtriggerEffnoHEEPBef.runMax              = cms.int32(   191717)


hNuEtriggerEffOneJnoHEEPBef                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJnoHEEPBef.numOfflJets     = cms.int32(  1 )
hNuEtriggerEffOneJnoHEEPBef.heepVersion     = cms.int32(  0 )
hNuEtriggerEffOneJnoHEEPBef.runMax          = cms.int32(   191717)


hNuEtriggerEffTwoJnoHEEPBef                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJnoHEEPBef.numOfflJets     = cms.int32(  2 )
hNuEtriggerEffTwoJnoHEEPBef.heepVersion     = cms.int32(  0 )
hNuEtriggerEffOneJnoHEEPBef.runMax          = cms.int32(   191717)






hNuEtriggerEffAft                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffAft.numOfflJets         = cms.int32(  0 )
hNuEtriggerEffAft.runMin       = cms.int32(   191718)

hNuEtriggerEffOneJAft                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJAft.numOfflJets     = cms.int32(  1 )
hNuEtriggerEffOneJAft.runMin       = cms.int32(   191718)


hNuEtriggerEffTwoJAft                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJAft.numOfflJets     = cms.int32(  2 )
hNuEtriggerEffTwoJAft.runMin       = cms.int32(   191718)


hNuEtriggerEffnoHEEPAft                     = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffnoHEEPAft.numOfflJets         = cms.int32(  0 )
hNuEtriggerEffnoHEEPAft.heepVersion         = cms.int32(  0 )
hNuEtriggerEffnoHEEPAft.runMin             = cms.int32(   191718)


hNuEtriggerEffOneJnoHEEPAft                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffOneJnoHEEPAft.numOfflJets     = cms.int32(  1 )
hNuEtriggerEffOneJnoHEEPAft.heepVersion     = cms.int32(  0 )
hNuEtriggerEffOneJnoHEEPAft.runMin             = cms.int32(   191718)


hNuEtriggerEffTwoJnoHEEPAft                 = hNuEleTriggerEff.HeavyNuEleTriggerEff.clone()
hNuEtriggerEffTwoJnoHEEPAft.numOfflJets     = cms.int32(  2 )
hNuEtriggerEffTwoJnoHEEPAft.heepVersion     = cms.int32(  0 )
hNuEtriggerEffTwoJnoHEEPAft.runMin             = cms.int32(   191718)





TriggerStudyElectronSequence = cms.Sequence(
    #kt6PFJetsForIsolation
    #dumpEvContent *

    hNuEtriggerEff
    * hNuEtriggerEffOneJ
    * hNuEtriggerEffTwoJ
    * hNuEtriggerEffnoHEEP
    * hNuEtriggerEffOneJnoHEEP
    * hNuEtriggerEffTwoJnoHEEP


    * hNuEtriggerEffBef
    * hNuEtriggerEffOneJBef
    * hNuEtriggerEffTwoJBef
    * hNuEtriggerEffnoHEEPBef
    * hNuEtriggerEffOneJnoHEEPBef
    * hNuEtriggerEffTwoJnoHEEPBef

    * hNuEtriggerEffAft
    * hNuEtriggerEffOneJAft
    * hNuEtriggerEffTwoJAft
    * hNuEtriggerEffnoHEEPAft
    * hNuEtriggerEffOneJnoHEEPAft
    * hNuEtriggerEffTwoJnoHEEPAft

    )
