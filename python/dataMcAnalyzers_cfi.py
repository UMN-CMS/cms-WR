import FWCore.ParameterSet.Config as cms

mixedFlavorSignalRegionAnalyzerLeadMuSubleadEle = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.),
		treeName = cms.string("mixedFlavorSignalRegionTree"),
		leptonsOneCollection = cms.InputTag("wRsubleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRleadingMuon"),
		jetsCollection = cms.InputTag("wRJets"),
		ignoreJets = cms.bool(False),
		checkThisSelector = cms.InputTag("flavourSidebandSelector"),
		dontCheckSelector = cms.InputTag(False)
		)


mixedFlavorSignalRegionAnalyzerLeadEleSubleadMu = mixedFlavorSignalRegionAnalyzerLeadMuSubleadEle.clone(leptonsOneCollection = cms.InputTag("wRleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingMuon")
		)

zToEEAnalyzer = mixedFlavorSignalRegionAnalyzerLeadMuSubleadEle.clone(treeName = cms.string("zEECheckTree"),
		leptonsOneCollection = cms.InputTag("wRHEEPElectron"),
		leptonsTwoCollection = cms.InputTag("wRHEEPElectron"),
		jetsCollection = cms.InputTag(""),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		ignoreJets = cms.bool(True),
		checkThisSelector = cms.InputTag(""),
		dontCheckSelector = cms.bool(True)
		)

zToMuMuAnalyzer = zToEEAnalyzer.clone(treeName = cms.string("zMuMuCheckTree"),
		leptonsOneCollection = cms.InputTag("scaleCorrectedMuons"),
		leptonsTwoCollection = cms.InputTag("scaleCorrectedMuons")
		)

lowDileptonMassAnalyzerMuMu = mixedFlavorSignalRegionAnalyzerLeadMuSubleadEle.clone(treeName = cms.string("lowMmumuTree"),
		leptonsOneCollection = cms.InputTag("wRleadingMuon"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingMuon"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		checkThisSelector = cms.InputTag("lowDiLeptonSidebandSelector")
		)

lowDileptonMassAnalyzerEE = lowDileptonMassAnalyzerMuMu.clone(treeName = cms.string("lowMeeTree"),
		leptonsOneCollection = cms.InputTag("wRleadingElectron"),
		leptonsTwoCollection = cms.InputTag("wRsubleadingElectron")
		)

lowFourObjectMassAnalyzerEE = lowDileptonMassAnalyzerEE.clone(treeName = cms.string("lowMeejjTree"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.),
		checkThisSelector = cms.InputTag("lowFourObjectSidebandSelector")
		)

lowFourObjectMassAnalyzerMuMu = lowFourObjectMassAnalyzerEE.clone(treeName = cms.string("lowMmumujjTree"))

sameFlavorSignalRegionAnalyzerEE = lowFourObjectMassAnalyzerEE.clone(treeName = cms.string("eeSignalRegionTree"),
		checkThisSelector = cms.InputTag("signalRegionSelector")
		)
sameFlavorSignalRegionAnalyzerMuMu = sameFlavorSignalRegionAnalyzerEE.clone(treeName = cms.string("mumuSignalRegionTree"))




