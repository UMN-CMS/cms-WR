import FWCore.ParameterSet.Config as cms
import glob

l = glob.glob('/uscms/home/jchaves/eos/ZRToNuMuMu_MZ0800_MNR0400/ZRToNuMuMu_MZR0800_MNR0400/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/160715_140213/0000/mini*.root')
fs = [x.replace('/uscms','file:/uscms') for x in l]
process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:/uscms/home/jchaves/nobackup/DYamc80x_2.root','file:/uscms/home/jchaves/nobackup/DYamc80x_3.root'),
                            )
outfile = 'puAna.root'

process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))


process.ana = cms.EDAnalyzer('PUAnalyzer',                             
                             )
                             

process.p = cms.Path(process.ana)
