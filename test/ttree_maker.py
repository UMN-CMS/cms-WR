import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/10000/1AEFBE02-4C02-E511-B796-0025905A60BE.root'),
                            #fileNames=cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/WRToNuMuToMuMuJJ_MW-2600_MNu-1300_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/102CF7B3-1F08-E511-93BB-00074305CF52.root'),
                            )
sec_file = open('miniaod_files.txt', 'r')
mysecfilelist = []
for line in sec_file:
    # add as many files as you wish this way
    mysecfilelist.append(line.strip())
process.source.fileNames = mysecfilelist
process.TFileService = cms.Service('TFileService', fileName = cms.string('made_ttree.root'))

process.wRtunePMuons = cms.EDProducer("TunePMuonProducer",
                             src = cms.InputTag("slimmedMuons")
                             )

process.ft_nocuts = cms.EDAnalyzer('TTreeMaker',
                            muons_src = cms.InputTag('wRtunePMuons'),
                            jets_src = cms.InputTag('slimmedJets'),
                            genparticles_src = cms.InputTag('prunedGenParticles'),
                            genjets_src = cms.InputTag('slimmedGenJets'),
                            primary_vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                            isolation_dR = cms.double(0.0),
                            leading_lepton_pt_cut = cms.double(0),
                            lepton_eta_cut = cms.double(5),
                            subleading_lepton_pt_cut = cms.double(0),
                            leading_jet_pt_cut = cms.double(0),
                            jet_eta_cut = cms.double(5),
                            subleading_jet_pt_cut = cms.double(0),
                            Mlljj_cut = cms.double(0.0),
                            Mll_cut = cms.double(0.0)
                            )
process.ft_basecuts = process.ft_nocuts.clone(  isolation_dR = cms.double(0.4),
                                                leading_lepton_pt_cut = cms.double(60),
                                                lepton_eta_cut = cms.double(2.4),
                                                subleading_lepton_pt_cut = cms.double(40),
                                                leading_jet_pt_cut = cms.double(40),
                                                jet_eta_cut = cms.double(2.5),
                                                subleading_jet_pt_cut = cms.double(40)
                                                )
process.p = cms.Path(process.wRtunePMuons + process.ft_nocuts)
process.p2 = cms.Path(process.wRtunePMuons + process.ft_basecuts)
#process.p3 = cms.Path(process.wRtunePMuons + process.ft_Mllcut)
#process.p4 = cms.Path(process.wRtunePMuons + process.ft_Mlljjcut)
