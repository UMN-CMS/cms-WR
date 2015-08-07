import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection

"""
To use this module:

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")
"""

def loadHEEPIDSelector(process):
	dataFormat = DataFormat.MiniAOD
	switchOnVIDElectronIdProducer(process, dataFormat)
	setupAllVIDIdsInModule(process,'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff',setupVIDElectronSelection)


HEEPIDSelector = cms.EDProducer('HEEPIDSelector',
      electrons= cms.InputTag("slimmedElectrons"),
      eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51")
      )

HEEPIDFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("HEEPIDSelector"),
		minNumber = cms.uint32(1)
		)

HEEPIDSequence = cms.Sequence(HEEPIDSelector * HEEPIDFilter)
