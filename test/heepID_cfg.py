import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    'root://cmsxrootd.fnal.gov///store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronSignalRegionSkim_1.root'
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "electron_ntuple_mini.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up electron ID (VID framework)
#

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#
# Configure an example module for user analysis with electrons
#

process.ntupler = cms.EDAnalyzer('heepIDAnalyzer',
                                 electronsMiniAOD    = cms.InputTag("HEEPIDSelector"),
                                 eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.HEEPIDSequence * process.ntupler)
