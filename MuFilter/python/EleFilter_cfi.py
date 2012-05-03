import FWCore.ParameterSet.Config as cms

eleFilter = cms.EDFilter("EleFilter",
                         electronTag  = cms.InputTag("patElectrons"),
                         ebTag    = cms.InputTag("correctedHybridSuperClusters"),
                         eeTag    = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                         
                         # require minEleNum gsf electrons above minElePt
                         minElePt  = cms.double( 30.0 ),
                         minEleNum = cms.double( 1 ),
                         
                         # additionally,
                         # require minSCNum above minSCEt separated from gfs (selected above) by at least DR>overlap
                         minSCEt      = cms.double( 30.0 ),
                         minSCNum     = cms.double( 1 ),
                         maxSCAbsEta  = cms.double( 2.6 ),
                         
                         # count SC's only if separate from any selected gsf by at least DR>overlap
                         overlap   = cms.double( 0.1 ),

                         massMin   = cms.double( 0.   ),
                         massMax   = cms.double( 8000.),
                         
                         debug     = cms.bool(False)

                         )
