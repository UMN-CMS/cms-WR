import FWCore.ParameterSet.Config as cms
import os, sys, imp, re

process = cms.Process('SELECTION')

############################################################ OPTIONS
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 

options.register('isMC',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "")
options.register('test',
                0,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "define the test type: 0=data, 1=signalMC, 2=background MC, 3=local file called skim_test.root")

options.register('datasetTag',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "unique dataset identifier")
options.register('unblind',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "0=blinded, 1=unblinded")
options.register('jsonFile',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "path and name of the json file")
options.register('runHLT',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "")

#default options
options.maxEvents = -1
defaultFileOutput = "myfile.root"
options.output = defaultFileOutput
#

options.parseArguments()
if(options.test==4):
    options.files="file:/uscms/home/jchaves/nobackup/DYamc80x_1.root "
    options.maxEvents=10000
    options.isMC=1
if(options.test==3):
    options.files="file:/uscms/home/jchaves/nobackup/singleMuB_80X_1.root"
    #options.maxEvents=100
    options.isMC=0
elif(options.test==2):
	options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_1.root"
	options.maxEvents= -1
	options.isMC=1
	options.datasetTag='dyJetsAmcNloInclusive'
elif(options.test==1):
    options.files='root://xrootd-cms.infn.it//store/mc/RunIIFall15MiniAODv2/WRToNuMuToMuMuJJ_MW-1000_MNu-500_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/38CDE93A-C2B8-E511-BEE7-002590FD5A72.root'
    options.maxEvents=200
    options.isMC=1

print options

if(options.output == defaultFileOutput and options.datasetTag!=''):
    options.output = options.datasetTag + '.root'

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load('ExoAnalysis.cmsWR.produceStringTag_cfi')
process.load('ExoAnalysis.cmsWR.pileupWeight_cff')

process.addStringIdentifier.stringStoredInOutputCollection = cms.string(options.datasetTag)

### \todo set the global tag in a separate file such that it will be common to all cfg files
if(options.isMC==0):
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'
else:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.options = cms.untracked.PSet(
#    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True),
#    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.files),
                            secondaryFileNames = cms.untracked.vstring(options.secondaryFiles)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService = cms.Service('TFileService', fileName = cms.string(options.output))


if(len(options.jsonFile)>0):
    print "[INFO] Using json file"
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.jsonFile).getVLuminosityBlockRange()

############################################################ OUTPUT MODULES
# this module defines the event content of our microAOD
process.load('ExoAnalysis.cmsWR.microAOD_Output_cff')

SelectEventsPSet = cms.untracked.PSet(
    SelectEvents = cms.vstring( [ 'FlavourSideband', 'SignalRegionEE', 'SignalRegionMuMu', 'LowDiLeptonSideband', 'DYtagAndProbe' ] )
    )


#define a process attribute for outputting a file which will be changed in a clone() call below
# process.microAOD_output = cms.OutputModule("PoolOutputModule",
# 		compressionAlgorithm = cms.untracked.string('LZMA'),
# 		compressionLevel = cms.untracked.int32(4),
# 		dataset = cms.untracked.PSet(
# 			dataTier = cms.untracked.string('MINIAODSIM'),
# 			filterName = cms.untracked.string('')
# 			),
# 		dropMetaData = cms.untracked.string('ALL'),
# 		eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
# 		fastCloning = cms.untracked.bool(False),
# 		fileName = cms.untracked.string(options.output),
# 		outputCommands = process.MICROAODSIMEventContent.outputCommands + [ 'keep *_*_*_SELECTION' ],
# 		overrideInputFileSplitLevels = cms.untracked.bool(True),
# 		SelectEvents = SelectEventsPSet
# 		)


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsPuppi'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
    )

if (options.isMC==0):
    updateJetCollection.jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L2L3Residual', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!



# here the full set of sequences and hlt paths used to make the first step   
process.load('ExoAnalysis.cmsWR.selections_cff')
from ExoAnalysis.cmsWR.JEC_cff import * # \todo check if this is needed
process.load('ExoAnalysis.cmsWR.treeMaker_cff')
process.load('ExoAnalysis.cmsWR.minitree_cfi')
process.load('ExoAnalysis.cmsWR.hltFilters_cff')
process.load('ExoAnalysis.cmsWR.heepSelector_cfi')

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)

process.load('ExoAnalysis.cmsWR.dataMcAnalyzers_cfi')

##############################################
##############################################

files = {"Prompt2015":"EgammaAnalysis/ElectronTools/data/ScalesSmearings/74X_Prompt_2015",
         "76XReReco" :"EgammaAnalysis/ElectronTools/data/ScalesSmearings/76X_16DecRereco_2015_Etunc",
         "80Xapproval" : "EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV1_2016_ele"}

correctionType = "80Xapproval"

process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3'),
        ),
    calibratedPatPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3'),
        ),
    )

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("wRHEEPElectron"),
                                         cut = cms.string("pt > 5 && abs(eta)<2.5") )
    
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

process.calibratedPatElectrons.electrons = "selectedElectrons"
process.calibratedPatElectrons.isMC = True
if (options.isMC==0):
    process.calibratedPatElectrons.isMC = False

##############################################
##############################################

process.blindSeq = cms.Sequence()
#process.dumperSeq = cms.Sequence(process.MakeTTree_Muons)
process.miniTTreeSeq = cms.Sequence(process.MiniTTree)
process.fullSeq = cms.Sequence(process.egmGsfElectronIDSequence * process.addStringIdentifier * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.jecSequence * process.electronHEEPSeq * process.selectedElectrons * process.calibratedPatElectrons * process.selectionSequence * process.filterSequence)

# Temporary while new MC is produced with HLT
if (options.isMC==0):
    process.signalHltSequence = cms.Sequence(process.wRHLTFilter_data)
    process.tagAndProbeHLTFilter = cms.Sequence(process.tagAndProbeHLTFilter_data)
else:
    if (options.runHLT==1):
        process.signalHltSequence = cms.Sequence(process.wRHLTFilter_MC)
        process.tagAndProbeHLTFilter = cms.Sequence(process.tagAndProbeHLTFilter_MC)
    else:
        process.signalHltSequence = cms.Sequence()
        process.tagAndProbeHLTFilter = cms.Sequence()

process.miniTree_signal_ee   = process.MiniTTree.clone()
process.miniTree_signal_mumu = process.MiniTTree.clone()
process.miniTree_flavoursideband = process.MiniTTree.clone()
process.miniTree_lowdileptonsideband = process.MiniTTree.clone()
process.miniTree_dytagandprobe = process.MiniTTree.clone()

############################################################ PATHs definition
process.SignalRegionEE      = cms.Path(process.signalHltSequence * process.fullSeq * process.blindSeq * process.signalRegionFilter * process.signalRegionEEFilter   * process.miniTree_signal_ee)
process.SignalRegionMuMu    = cms.Path(process.signalHltSequence * process.fullSeq * process.blindSeq * process.signalRegionFilter * process.signalRegionMuMuFilter * process.miniTree_signal_mumu)
process.FlavourSideband     = cms.Path(process.signalHltSequence * process.fullSeq                   * ~process.signalRegionFilter * process.flavourSidebandFilter * process.miniTree_flavoursideband)
process.LowDiLeptonSideband = cms.Path(process.signalHltSequence * process.fullSeq                   * ~process.signalRegionFilter * process.lowDiLeptonSidebandFilter * process.miniTree_lowdileptonsideband)
#process.LowMassSideband    = cms.Path(process.signalHltSequence * process.fullSeq * process.blindSeq * process.signalRegionFilter * process.miniTree_signal)

process.DYtagAndProbe = cms.Path(process.tagAndProbeHLTFilter * process.egmGsfElectronIDSequence * process.addStringIdentifier * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.jecSequence * process.electronHEEPSeq * process.selectedElectrons * process.calibratedPatElectrons * process.selectionSequence * process.miniTree_dytagandprobe * process.zToEEAnalyzer * process.zToMuMuAnalyzer)

#process.microAODoutput_step = cms.EndPath(process.microAOD_output)

############################################################ SCHEDULE
if (options.isMC==0 and options.unblind==0):
    process.blindSeq += ~process.signalRegionFilter
    print "########################################################################"
    print "# WARNING!!! You are running on DATA, but the analysis is still BLIND! #"
    print "# The signal region path will not be run!                              #"
    print "########################################################################"

    process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.DYtagAndProbe)
else:
    process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.SignalRegionEE, process.SignalRegionMuMu, process.DYtagAndProbe) #, process.microAODoutput_step)
#    process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.SignalRegionEE, process.SignalRegionMuMu, process.DYtagAndProbe, process.microAODoutput_step)


