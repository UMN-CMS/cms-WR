import os,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('ExoAnalysis.cmsWR.treeMaker_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:/afs/cern.ch/user/j/jchavesb/work/WR_skims/skim_sideband_2600.root'),
                            )

os.system('das_client.py --query="file dataset=/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/jchavesb-crab_WR_TTJets_Skim-8a626735db4a89bf68c59f2da6d5e98a/USER instance=prod/phys03" --limit=0 > pappy.txt')

#sec_file = open('pappy.txt', 'r')
#sec_file = open('dyskim.txt', 'r')
sec_file = open('ttskim.txt', 'r')
mysecfilelist = []
for line in sec_file:
    # add as many files as you wish this way
    if 'skim' in line:
        mysecfilelist.append('file:'+line.strip())
#process.source.fileNames = mysecfilelist

process.TFileService = cms.Service('TFileService', fileName = cms.string('skim_ttree.root'))

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.triggerFilter = hltHighLevel.clone()
process.triggerFilter.HLTPaths = ['HLT_Mu45_eta2p1_v1','HLT_Mu50_v1']
process.triggerFilter.andOr = True # = OR
#for name, path in process.paths.items():
 #   if not name.startswith('eventCleaning'):
  #      path.insert(0, process.triggerFilter)
process.ptrig = cms.Path(process.triggerFilter)

process.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ptrig'))

process.p = cms.Path(process.wRtunePMuons + process.wRsubleadingMuon + process.wRlooseJet + process.MakeTTree_Muons)
