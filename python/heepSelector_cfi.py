import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection

"""
To use this module:

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

then run the module named egmGsfElectronIDSequence before running HEEPIDSequence in the path
"""

def loadHEEPIDSelector(process):
	dataFormat = DataFormat.MiniAOD
	switchOnVIDElectronIdProducer(process, dataFormat)
	setupAllVIDIdsInModule(process,'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',setupVIDElectronSelection)



HEEPIDFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("HEEPIDSelector"),
		minNumber = cms.uint32(2)
		)

HEEPIDSequence = cms.Sequence(HEEPIDSelector * HEEPIDFilter)

#only require that one electron in the low mass sideband pass HEEP ID
#HEEPIDSidebandFilter = HEEPIDFilter.clone(minNumber = cms.uint32(1) )
#HEEPIDSidebandSequence = cms.Sequence(HEEPIDSelector * HEEPIDSidebandFilter)
