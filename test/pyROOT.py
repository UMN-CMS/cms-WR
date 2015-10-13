import ROOT
from ROOT import *
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *
#import print_options
#print_options.set_float_precision(4)
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()
#import EcalDetId
#from DataFormats.EcalDetId import *
#
import sys,os

#file='EXO-Phys14DR-00009.root'
#events = Events (file)
#
#print file
#handleElectrons = Handle('std::vector<pat::Electron>')
#electronTAG = "slimmedElectrons"
#
#
#print "numEvents = ", events.size()
#n=0
#for event in events:
#    print n
#    event.getByLabel(electronTAG, handleElectrons)
#    electrons = handleElectrons.product()
#    print electrons.size()
#    eleIndex=0
#    for electron in electrons:
#        
#        print eleIndex, electron.pt(), electron.eta(), electron.phi()
#        eleIndex+=1
#    
#    n+=1

file='/uscms/home/skalafut/nobackup/WR_starting2015/GEN-SIM_13TeV/00774ADC-8F87-E411-B6EF-00266CF32E78.root'

handleGen = Handle('std::vector<reco::GenParticle>')
genTAG = "genParticles"

events = Events(file)

stopIndex = 100
index=0
for event in events:
	event.getByLabel(genTAG,handleGen)
	genParts = handleGen.product()
	for gen in genParts:
		if(abs(gen.pdgId())!=9900024 and gen.mass() > 1200):
			print gen.pdgId(), gen.mass(), gen.status(), gen.charge(), gen.numberOfDaughters()
		#if(abs(gen.pdgId())==9900024):
		#	print gen.pdgId(), gen.mass(), gen.status(), gen.charge(), gen.numberOfDaughters()
	
	index+=1
	if(index>=stopIndex):
		break


