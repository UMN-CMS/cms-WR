#this cfg is designed to run over MINIAODSIM specifically to produce a file with one minitree which contains the nPU distribution
#the puweight distribution should be ignored
#all of the code in the analyze(edm::Event& event, edm::EventSetup&) method in miniTTree.cc which is not used to obtain
#nPU MUST be commented out before running this cfg file with cmsRun
#run this file and the associated analyzer from the cmsWR dir by executing  cmsRun test/runAnalysisForNPu_cfg.py
import FWCore.ParameterSet.Config as cms
import os, sys, imp, re

process = cms.Process('NUMPU')

############################################################ OPTIONS
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 

options.register('isMC',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "")
options.register('test',
                2,
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

#default options
options.maxEvents = -1
defaultFileOutput = "myfile.root"
options.output = defaultFileOutput
#

options.parseArguments()
if(options.test==3):
    options.files="file:skim_test.root"   
    #options.maxEvents=100
    options.isMC=1
elif(options.test==2):
	options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_1.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_2.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_3.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_4.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_5.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_6.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_7.root","file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_8.root"
	#options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_1.root"
	#options.maxEvents= -1
	options.isMC=1
	options.datasetTag='dyJetsAmcNloInclusiveNoSelectionCutsForPUweights'
	#options.output = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveWithNumPUandNoSelections.root"
elif(options.test==1):
    options.files='root://eoscms//eos/cms/store/user/shervin/WRToNuMuToMuMuJJ_MW-2600_MNu-1300_TuneCUETP8M1_13TeV-pythia8/WRtoMuMuJJ_2600_1300_SHv2/160124_160701/0000/output_1.root'
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
    process.GlobalTag.globaltag = '74X_dataRun2_v5'
else:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'


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

#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TFileService = cms.Service('TFileService', fileName = cms.string(options.output))


if(len(options.jsonFile)>0):
    print "[INFO] Using json file"
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.jsonFile).getVLuminosityBlockRange()

############################################################ MODULES
process.load('ExoAnalysis.cmsWR.minitree_cfi')

process.miniTTreeSeq = cms.Sequence(process.MiniTTree)

process.miniTree_noSelections = process.MiniTTree.clone(
		muons_src = cms.InputTag(''),
		electrons_src = cms.InputTag(''),
		jets_src = cms.InputTag(''),
		jec_unc_src = cms.InputTag(''),
		muon_IDSF_central_src = cms.InputTag(''),
		muon_IsoSF_central_src = cms.InputTag(''),
		muon_IDSF_error_src = cms.InputTag(''),
		muon_IsoSF_error_src = cms.InputTag(''),
		PUWeights_src = cms.InputTag('PUWeights','PileupWeights'),
		datasetName = cms.InputTag("addStringIdentifier", "datasetIdentifier")
		)


############################################################ PATHs definition
process.NoSelectionsForNPu = cms.Path(process.addStringIdentifier * process.PUWeightsSequence * process.miniTree_noSelections)


############################################################ SCHEDULE

process.schedule = cms.Schedule(process.NoSelectionsForNPu)


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
CMSSW_BASE=os.getenv("CMSSW_BASE")

pathPrefix=CMSSW_BASE+'/src/ExoAnalysis/cmsWR/'

process.PUWeights.PileupMCFilename = cms.string(pathPrefix + "data/MCPileup.root")
process.PUWeights.PileupDataFilename = cms.string(pathPrefix + "data/DataPileup.root")

