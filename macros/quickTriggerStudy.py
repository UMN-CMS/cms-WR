# use this python program to see which HLT paths are fired most often in WR->ENu->EEJJ and WR->MuNu->MuMuJJ evts

# sys.argv is a list of command line arguments passed to this program
#import sys
#oldargv = sys.argv[:]
#sys.argv = ['-b-']
import ROOT
import array
import math
#import dict
#ROOT.gROOT.SetBatch(True)  #run root without graphics
#sys.argv = oldargv

# load FWLite C++ libs
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# load FWLite python libs
from DataFormats.FWLite import Handle, Events

# define additional fxns

#hltPathName is a string
def skipThisPath(hltPathName):
	doSkip = False
	stringsToAvoid = ["Jet","hoton","MET","HT","AlCa","Flag","step","BTag","Rsq","ReducedIterative","Physics","L1","ZeroBias","Onia","Phi","Random","Tau","Ecal","Final","Hcal","Upsilon","Jpsi","Displaced","NoFilters","L2","SameSign","TripleMu"]
	for i in xrange(int(len(stringsToAvoid))):
		if(hltPathName.find(stringsToAvoid[i]) != -1):
			#if one of the strings which should be avoided is found in the hlt path name, then break out of this
			#for loop and return True
			doSkip = True
			break
	return doSkip

#return uncertainty on efficiency
#effNumerator and effDenominator are floats
def calcEffUnc(effNumerator, effDenominator):
	uncertainty = (1/effDenominator)*math.sqrt(effNumerator*(1-(effNumerator/effDenominator)))
	return uncertainty
#end calcEffUnc

#check that the correct WR decay products are found at GEN lvl
#genHandle is an edm::Handle object linked to the prunedGenParticles collection
#doMuon is True or False to indicate the final state of interest (True --> mumuqq, False --> eeqq)
def foundCorrectWrDecayProducts(genHandle, doMuon):
	#print 'in foundCorrectWrDecayProducts fxn'
	foundAllDecayProducts = False
	maxEtaQrk, maxEtaLept = 99, 99 
	leptPdg, hvyNuPdgId, wRpdgId = 99, 100, 9900024
	if(doMuon == True): leptPdg, hvyNuPdgId, maxEtaLept = 13, 9900014, 2.4
	else: leptPdg, hvyNuPdgId = 11, 9900012
	#now the lepton and hvyNu pdgId are set to the desired values
	
	nWrLeptonDaus, nNuLeptonDaus, nNuQuarkDaus = 0,0,0
	#loop over all particles in genHandle and
	#look for at least one lepton from WR, one lepton from Nu, and two quarks from Nu
	genList = genHandle.product() 	#genList is a vector (treat like a list)
	#print 'there are this many gen particles in the event: ', genList.size()
	for i in xrange(int(genList.size())):
		if( abs(genList[i].pdgId()) == leptPdg and abs(genList[i].mother(0).pdgId()) == wRpdgId and abs(genList[i].eta())<=maxEtaLept): nWrLeptonDaus+=1
		if( abs(genList[i].pdgId()) == leptPdg and abs(genList[i].mother(0).pdgId()) == hvyNuPdgId and abs(genList[i].eta())<=maxEtaLept): nNuLeptonDaus+=1
		if( abs(genList[i].pdgId()) <= 6 and abs(genList[i].mother(0).pdgId()) == hvyNuPdgId and abs(genList[i].eta())<=maxEtaQrk): nNuQuarkDaus+=1

	#print 'nWrLeptonDaus = ', nWrLeptonDaus, ' , nNuLeptonDaus = ', nNuLeptonDaus,' , nNuQuarkDaus = ', nNuQuarkDaus
	if(nWrLeptonDaus > 0 and nNuLeptonDaus > 0 and nNuQuarkDaus > 1): foundAllDecayProducts = True

	return foundAllDecayProducts
#end foundCorrectWrDecayProducts()


# define handles and labels for TriggerResults, TriggerObjectStandAlone, and GenParticle collections
trigResultsHandl, trigResultsLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
trigObjsHandl, trigObjsLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
genParticleHandl, genParticleLabel = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"

# open one or more miniAOD .root files, and create an iterator to loop over the Event objects in the file(s)
allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")
#allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")
doMuonChannel = False


#vars and containers which will be filled in loop over events
trigNamesAndNumPassing = dict()
totalNumEvts = 0
totalNumEvtsAfterGenCuts = 0

# loop over Event objects in input file(s)
# enumerate returns two-tuple objects (an immutable list where each element in the list is a pair), and thus two different
# iterator variables are needed in the for loop declaration
# the first value of each pair in the two tuple is a simple number index starting at 0, and the second value in each pair
# in the two tuple is the object of interest
for evNum, oneEvent in enumerate(allEvents):
	#evNum points to a simple number, while oneEvent points to an edm::Event object
	#if(evNum > 5): break
	totalNumEvts+=1.0
	
	#link the reco::GenParticle collection with an edm::Handle 
	oneEvent.getByLabel(genParticleLabel, genParticleHandl)
	#check that the correct WR decay products are found at GEN level
	if(foundCorrectWrDecayProducts(genParticleHandl, doMuonChannel) == False): continue
	totalNumEvtsAfterGenCuts+=1.0
	oneEvent.getByLabel(trigObjsLabel, trigObjsHandl)
	oneEvent.getByLabel(trigResultsLabel, trigResultsHandl)

	if((evNum%1000)==0): print "evt number: ", evNum
	
	allPathNames = oneEvent.object().triggerNames(trigResultsHandl.product())
	for i in xrange(trigResultsHandl.product().size()):
		#loop over all HLT path names and keep track of the paths which are fired, how many times they fire, and the total number of evts
		#which have been run over
		#print "Trigger ", allPathNames.triggerName(i), (" PASS" if trigResultsHandl.product().accept(i) else " fail (or not run)")
		if(trigResultsHandl.product().accept(i) and skipThisPath(allPathNames.triggerName(i)) == False ):
			if(allPathNames.triggerName(i) in trigNamesAndNumPassing):
				#the trigger path name already exists in the map trigNamesAndNumPassing
				#increment the numerical value stored in this map entry by 1
				oldNumPassing = trigNamesAndNumPassing[allPathNames.triggerName(i)]
				newNumPassing = oldNumPassing+1.0
				trigNamesAndNumPassing[allPathNames.triggerName(i)] = newNumPassing

			else:
				#the trigger path name does not exist in the map trigNamesAndNumPassing
				#add the path name to the map as a key, and set the map value to 1
				trigNamesAndNumPassing[allPathNames.triggerName(i)] = 1.0

	#end loop over HLT path names
# end loop over evts in input file(s)

print "total num evts before any cuts = ", totalNumEvts
print 'total num evts passing gen cuts = ', totalNumEvtsAfterGenCuts
for keyIter in iter(trigNamesAndNumPassing):
	if( (trigNamesAndNumPassing[keyIter]/totalNumEvtsAfterGenCuts) > 0.60):
		print "HLT path named: ", keyIter," has trigger efficiency = ", 100*(trigNamesAndNumPassing[keyIter]/totalNumEvtsAfterGenCuts) ,"% +/- ", 100*(calcEffUnc(trigNamesAndNumPassing[keyIter],totalNumEvtsAfterGenCuts)), "% in evts where gen requirements are passed"
	
#end loop over keys in trigNamesAndNumPassing map

print 'gen requirements are that the four final state particles at GEN lvl from the WR decay are found, and fall within detector eta acceptance'





