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

file='EXO-Phys14DR-00009.root'
events = Events (file)

print file
handleElectrons = Handle('std::vector<pat::Electron>')
electronTAG = "slimmedElectrons"

print "numEvents = ", events.size()
n=0
for event in events:
    print n
    event.getByLabel(electronTAG, handleElectrons)
    electrons = handleElectrons.product()
    print electrons.size()
    eleIndex=0
    for electron in electrons:
        
        print eleIndex, electron.pt(), electron.eta(), electron.phi()
        eleIndex+=1
    
    n+=1
    

