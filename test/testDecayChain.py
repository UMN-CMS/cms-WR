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

#from __future__ import print_function

file='/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-6000_MNu-3000_TuneCUETP8M1_pythia8_13TeV.root'
events = Events (file)

print file

handleGenParticles = Handle('std::vector<reco::GenParticle>')
genParticleTAG = 'prunedGenParticles'

print "numEvents = ", events.size()


from collections import *

nWR = Counter()
nWRdaughter0 = Counter()
nNuR = Counter()
nWstarR = Counter()
nNuRdaughter0 = Counter()
nNuRdaughter3 = Counter()

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def genPrint(particle, level):
    for i in range(0,level):
        print " ",
    for i in range(0,level):
        print "-", 
        sys.stdout.write('')
    
    print ">",
    if(particle.status()==1):
        print bcolors.OKBLUE + " ",
    if(-9 < particle.pdgId() < 9):
        print bcolors.OKGREEN + " ",
        
    print particle.pdgId(), particle.mass(), particle.eta(), particle.pt(), particle.status(), particle.numberOfDaughters(), bcolors.ENDC+' '
    

def printDaughters(particle, level):
#    genPrint(particle, level)
    if(level<6 and particle.numberOfDaughters()>0):
        for kDaughter in range(0, particle.numberOfDaughters()):
            daughter = particle.daughter(kDaughter)
            genPrint(daughter, level)
            if((daughter.pdgId() > 9 or daughter.pdgId()<-9) and not daughter.pdgId()==1):
                printDaughters(daughter, level+1)

for event in events:
    event.getByLabel(genParticleTAG, handleGenParticles)
    genParticles = handleGenParticles.product()


    

    for genParticle in genParticles :
        if(genParticle.pdgId()==9900024 or genParticle.pdgId()==-9900024):

            # recursive function printing decay 
            printDaughters(genParticle,1) 

            nuR = genParticle.daughter(1)

            nWR[genParticle.pdgId()] += 1
            nWRdaughter0[genParticle.daughter(0).pdgId()] +=1
            nNuR[nuR.pdgId()] +=1
            nWstarR[nuR.daughter(1).pdgId()] += 1
            nNuRdaughter0[nuR.daughter(0).pdgId()] += 1
            if(nuR.numberOfDaughters()>2):
                nNuRdaughter3[nuR.daughter(2).pdgId()] += 1

            # nNuRdaugh1[nuR.daughter(0).pdgId()] +=1
            # nNuRdaugh2[nuR.daughter(1).pdgId()] +=1
            # nNuRdaugh3[nuR.daughter(2).pdgId()] +=1

            #     print "**   -> ", nuR.daughter(3).pdgId(), nuR.daughter(3).pt(), nuR.daughter(3).eta(), nuR.daughter(3).mass()
            #     nNuRdaugh4[nuR.daughter(3).pdgId()] +=1


            # print " -> ", genParticle.daughter(0).pdgId(), genParticle.daughter(0).pt(), genParticle.daughter(0).eta(), genParticle.daughter(0).mass()
            # print " -> ", genParticle.daughter(1).pdgId(), genParticle.daughter(1).pt(), genParticle.daughter(1).eta(), genParticle.daughter(1).mass()

            # print "   -> ", nuR.daughter(0).pdgId(), nuR.daughter(0).pt(), nuR.daughter(0).eta(), nuR.daughter(0).mass()
            # print "   -> ", nuR.daughter(1).pdgId(), nuR.daughter(1).pt(), nuR.daughter(1).eta(), nuR.daughter(1).mass()
            # print "   -> ", nuR.daughter(2).pdgId(), nuR.daughter(2).pt(), nuR.daughter(2).eta(), nuR.daughter(2).mass()
            # ultraDaughter = nuR.daughter(2)
            # if(ultraDaughter.pdgId() == 9900014 or ultraDaughter.pdgId() == 9900016):
            #     print "      -> ", ultraDaughter.daughter(0).pdgId(), ultraDaughter.daughter(0).pt(), ultraDaughter.daughter(0).eta(), ultraDaughter.daughter(0).mass()
            #     print "      -> ", ultraDaughter.daughter(1).pdgId(), ultraDaughter.daughter(1).pt(), ultraDaughter.daughter(1).eta(), ultraDaughter.daughter(1).mass()
            #     print "      -> ", ultraDaughter.daughter(2).pdgId(), ultraDaughter.daughter(2).pt(), ultraDaughter.daughter(2).eta(), ultraDaughter.daughter(2).mass()



                
print "=============================="
print nWR
print nWRdaughter0
print nNuR
print nNuRdaughter0
print nWstarR
print nNuRdaughter3
