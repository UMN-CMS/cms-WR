import FWCore.ParameterSet.Config as cms

mixedFlavorSignalRegionAnalyzer = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.),
		treeName = cms.string("mixedFlavorSignalRegionTree"),
		leptonsOneCollection = cms.InputTag("wRleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRleadingMuon"),
		jetsCollection = cms.InputTag("wRtightJets"),
		ignoreJets = cms.bool(False)
		)

zToEEAnalyzer = mixedFlavorSignalRegionAnalyzer.clone(treeName = cms.string("zEECheckTree"),
		leptonsOneCollection = cms.InputTag("wRleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingElectron"),
		jetsCollection = cms.InputTag(""),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		ignoreJets = cms.bool(True)
		)

zToMuMuAnalyzer = zToEEAnalyzer.clone(treeName = cms.string("zMuMuCheckTree"),
		leptonsOneCollection = cms.InputTag("scaleCorrectedMuons"),
		leptonsTwoCollection = cms.InputTag("scaleCorrectedMuons")
		)

lowDileptonMassAnalyzerMuMu = mixedFlavorSignalRegionAnalyzer.clone(treeName = cms.string("lowMmumuTree"),
		leptonsOneCollection = cms.InputTag("wRleadingMuon"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingMuon"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

lowDileptonMassAnalyzerEE = mixedFlavorSignalRegionAnalyzer.clone(treeName = cms.string("lowMeeTree"),
		leptonsOneCollection = cms.InputTag("wRleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingElectron"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

lowFourObjectMassAnalyzerEE = lowDileptonMassAnalyzerEE.clone(treeName = cms.string("lowMeejjTree"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.)
		)

lowFourObjectMassAnalyzerMuMu = lowDileptonMassAnalyzerMuMu.clone(treeName = cms.string("lowMmumujjTree"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.)
		)

sameFlavorSignalRegionAnalyzerEE = lowFourObjectMassAnalyzerEE.clone(treeName = cms.string("eeSignalRegionTree"))
sameFlavorSignalRegionAnalyzerMuMu = lowFourObjectMassAnalyzerMuMu.clone(treeName = cms.string("mumuSignalRegionTree"))




