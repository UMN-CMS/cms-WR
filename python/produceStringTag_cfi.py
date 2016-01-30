import FWCore.ParameterSet.Config as cms

addStringIdentifier = cms.EDProducer('produceStringTag',
                                     outputCollectionName = cms.string("datasetIdentifier"),
                                     stringStoredInOutputCollection = cms.string('')
                                     )
