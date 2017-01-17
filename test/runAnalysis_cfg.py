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

#default options
options.maxEvents = -1
defaultFileOutput = "DYJets_amctnlo.root"
options.output = defaultFileOutput
#

options.parseArguments()
if(options.test==3):
    options.files="file:skim_test.root"   
    #options.maxEvents=100
    options.isMC=1
elif(options.test==0):
	#options.files=['root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunC_SHv12/MuEG_RunC_1_1_n0i.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunC_SHv12/MuEG_RunC_2_1_ZKp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_10_1_dUq.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_11_1_gsr.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_12_1_tiu.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_13_2_wrN.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_14_1_Wal.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_15_1_q4V.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_16_1_eFY.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_17_1_eB8.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_18_1_Y24.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_1_1_9K5.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_2_1_cbs.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_3_1_lCY.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_4_1_ffI.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_5_1_dp2.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_6_1_NAp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_7_1_mgH.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_8_1_1fS.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v3_SHv12/MuEG_RunD_v3_9_2_Fza.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_10_1_NsR.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_11_1_lSG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_12_1_j6L.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_13_1_6k8.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_14_1_R0H.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_15_1_qAG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_16_1_1yS.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_17_1_CH1.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_18_1_KnO.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_19_1_uxo.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_1_1_sTf.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_20_1_yMq.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_21_1_fqI.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_22_1_DIO.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_23_1_Nht.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_24_1_Xia.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_25_1_60r.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_2_1_LdK.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_3_1_rXd.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_4_1_pFp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_5_2_usY.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_6_1_vXA.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_7_1_A7n.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_8_1_VKw.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/MuEG_RunD_v4_SHv12/MuEG_RunD_v4_9_1_vXW.root']
	options.maxEvents= -1 
	options.isMC=0
	options.datasetTag='reprocessMuEGData'
elif(options.test==2):
	options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYJets_amctnlo_skimFile.root"
	#options.files=['root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_10_2_KFJ.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_11_2_PR1.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_12_2_DuM.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_13_2_yHC.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_14_2_bWL.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_15_2_rK5.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_16_2_yJp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_17_2_Gff.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_18_2_KEW.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_19_2_3FO.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_1_2_7zA.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_20_2_hc0.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_21_2_nur.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_22_1_SSm.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_2_2_89c.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_3_2_yNe.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_4_2_aL0.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_5_2_RDU.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_6_2_wPy.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_7_2_YMG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_8_2_xNp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_9_2_EMH.root']
	#options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAOD2015/DYMadInclusiveMiniAOD2.root"
	options.maxEvents= -1 
	options.isMC=1
	options.datasetTag='DYJets_amctnlo'
	#options.datasetTag='reprocessTTJetsV2MCPartOne'
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


# here the full set of sequences and hlt paths used to make the first step
process.load('ExoAnalysis.cmsWR.selections_cff')
from ExoAnalysis.cmsWR.JEC_cff import * # \todo check if this is needed
process.load('ExoAnalysis.cmsWR.treeMaker_cff')
process.load('ExoAnalysis.cmsWR.minitree_cfi')
process.load('ExoAnalysis.cmsWR.hltFilters_cff')
process.load('ExoAnalysis.cmsWR.heepSelector_cfi')
process.load('ExoAnalysis.cmsWR.genDYJetsElectronChannelModules_cff')

#################################
#Filters
process.genHTFilter = cms.EDFilter('htFilter',
		cutThreshold = cms.double(100.),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)

process.load('ExoAnalysis.cmsWR.dataMcAnalyzers_cfi')


process.blindSeq = cms.Sequence()
#process.dumperSeq = cms.Sequence(process.MakeTTree_Muons)
process.miniTTreeSeq = cms.Sequence(process.MiniTTree)
#with HTFilter process.fullSeq = cms.Sequence(process.egmGsfElectronIDSequence * process.addStringIdentifier * process.PUWeightsSequence * process.jecSequence * process.selectionSequence  * process.filterSequence * process.dyJetsBareMatchedGenGluon * process.dyJetsBareMatchedGenQuark * process.dyJetsMergeGenMatchedPartons * process.genHTFilter)
process.fullSeq = cms.Sequence(process.egmGsfElectronIDSequence * process.addStringIdentifier * process.PUWeightsSequence * process.jecSequence * process.selectionSequence  * process.filterSequence)


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

#with HTFilter process.DYtagAndProbe = cms.Path(process.dyJetsBareMatchedGenGluon * process.dyJetsBareMatchedGenQuark * process.dyJetsMergeGenMatchedPartons * process.genHTFilter * process.tagAndProbeHLTFilter * process.egmGsfElectronIDSequence * process.addStringIdentifier * process.PUWeightsSequence * process.jecSequence * process.selectionSequence * process.miniTree_dytagandprobe * process.zToEEAnalyzer * process.zToMuMuAnalyzer)
process.DYtagAndProbe = cms.Path(process.tagAndProbeHLTFilter * process.egmGsfElectronIDSequence * process.addStringIdentifier * process.PUWeightsSequence * process.jecSequence * process.selectionSequence * process.miniTree_dytagandprobe * process.zToEEAnalyzer * process.zToMuMuAnalyzer)


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
	#process.schedule = cms.Schedule(process.FlavourSideband) #, process.microAODoutput_step)
#    process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.SignalRegionEE, process.SignalRegionMuMu, process.DYtagAndProbe, process.microAODoutput_step)


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
CMSSW_BASE=os.getenv("CMSSW_BASE")

pathPrefix=CMSSW_BASE+'/src/ExoAnalysis/cmsWR/'

process.PUWeights.PileupMCFilename = cms.string(pathPrefix + "data/MCPileup.root")
process.PUWeights.PileupDataFilename = cms.string(pathPrefix + "data/DataPileup.root")

