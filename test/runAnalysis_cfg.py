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
defaultFileOutput = "reprocessTTJetsV2MCPartOne.root"
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
	#original options.files="file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/puReweightingFiles/dyJetsAmcNloInclusiveM50MiniAodSim_1.root"
	options.files=['root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_11_2_n9j.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_10_2_K9s.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_12_2_Jjo.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_13_2_WGG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_14_2_8So.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_15_2_hhg.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_16_2_ltB.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_17_2_GWQ.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_18_2_GU1.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_19_2_ZXr.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_1_2_Ebe.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_20_2_dK5.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_21_2_zMB.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_22_2_xxx.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_23_2_t1T.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_24_2_6VU.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_25_2_i13.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_26_2_vM5.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_27_2_5dF.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_28_2_Eko.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_29_2_WsC.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_2_2_AFa.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_30_2_8cp.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_31_2_Q9v.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_32_2_4I4.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_33_1_Juk.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_34_2_Yo5.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_35_2_7wE.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_36_2_C8E.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_37_2_6Ck.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_38_1_wSb.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_39_2_oFG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_3_2_qre.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_40_2_iUh.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_41_2_vYb.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_42_2_4Kn.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_43_2_ebG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_44_2_TAg.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_45_2_1aF.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_46_2_78K.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_47_2_AUE.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_48_2_q4p.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_49_2_3xe.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_4_2_CVo.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_50_2_Fnx.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_51_2_Rwz.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_52_2_lma.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_53_2_guy.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_54_2_y7Q.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_55_2_CLA.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_56_2_GyH.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_57_2_bmj.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_58_2_x8V.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_59_2_qRW.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_5_2_TmO.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_60_2_6uE.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_6_2_UB3.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_7_2_jq3.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_8_2_nwG.root','root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/TTJets_DiLept_v2_SHv12/TTJets_DiLept_v2_9_2_C3A.root']
	options.maxEvents= -1 
	options.isMC=1
	#original options.datasetTag='dyJetsAmcNloInclusive'
	options.datasetTag='reprocessTTJetsV2MCPartOne'
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

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)

process.load('ExoAnalysis.cmsWR.dataMcAnalyzers_cfi')


process.blindSeq = cms.Sequence()
#process.dumperSeq = cms.Sequence(process.MakeTTree_Muons)
process.miniTTreeSeq = cms.Sequence(process.MiniTTree)
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
    #process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.SignalRegionEE, process.SignalRegionMuMu, process.DYtagAndProbe) #, process.microAODoutput_step)
	process.schedule = cms.Schedule(process.FlavourSideband) #, process.microAODoutput_step)
#    process.schedule = cms.Schedule(process.FlavourSideband, process.LowDiLeptonSideband, process.SignalRegionEE, process.SignalRegionMuMu, process.DYtagAndProbe, process.microAODoutput_step)


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
CMSSW_BASE=os.getenv("CMSSW_BASE")

pathPrefix=CMSSW_BASE+'/src/ExoAnalysis/cmsWR/'

process.PUWeights.PileupMCFilename = cms.string(pathPrefix + "data/MCPileup.root")
process.PUWeights.PileupDataFilename = cms.string(pathPrefix + "data/DataPileup.root")

