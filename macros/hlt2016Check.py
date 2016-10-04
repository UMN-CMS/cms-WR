# use this python program to check the efficiency of triggers in the 2016 HLT menu on miniAOD or miniAODSIM evts
# without any offline selection requirements

import ROOT
import array
import math

# load FWLite C++ libs
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# load FWLite python libs
from DataFormats.FWLite import Handle, Events

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


#switch this to True to study muon channel triggers
doMuonChannel = False 

# define handles and labels for TriggerResults, TriggerObjectStandAlone, and GenParticle collections
trigResultsHandl, trigResultsLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
trigObjsHandl, trigObjsLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
recoJetsHandl, recoJetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"

# open one or more miniAOD .root files, and create an iterator to loop over the Event objects in the file(s)
#absolute path name to input miniAOD or miniAODSIM file
allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")


#absolute file path name used for electron channel
if(doMuonChannel == False):
	allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")
	recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"



#vars and containers which will be filled in loop over events
trigNamesAndNumPassing = dict()
totalNumEvts = 0

# loop over Event objects in input file(s)
# enumerate returns two-tuple objects (an immutable list where each element in the list is a pair), and thus two different
# iterator variables are needed in the for loop declaration
# the first value of each pair in the two tuple is a simple number index starting at 0, and the second value in each pair
# in the two tuple is the object of interest
for evNum, oneEvent in enumerate(allEvents):
	#evNum points to a simple number, while oneEvent points to an edm::Event object
	#if(evNum > 2): break
	totalNumEvts+=1.0
	
	#link pat::Electron or pat::Muon collections with a edm::Handle objs
	oneEvent.getByLabel(recoLeptonLabel, recoLeptonHandl)
	oneEvent.getByLabel(recoJetsLabel, recoJetsHandl)
	
	oneEvent.getByLabel(trigObjsLabel, trigObjsHandl)
	oneEvent.getByLabel(trigResultsLabel, trigResultsHandl)

	if((evNum%1000)==0): print "evt number: ", evNum
	
	allPathNames = oneEvent.object().triggerNames(trigResultsHandl.product())
	for i in xrange(trigResultsHandl.product().size()):
		#loop over all HLT path names and keep track of the paths which are fired, how many times they fire, and the total number of evts
		#which have been run over
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

#save the tree with kinematics info to a file

#triggers which fire with efficiency below minimumEfficiency will not be printed
minimumEfficiency = 0.0001
print "total num evts processed = ", totalNumEvts
for keyIter in iter(trigNamesAndNumPassing):
	if( (trigNamesAndNumPassing[keyIter]/totalNumEvts) > minimumEfficiency):
		print "HLT path named: ", keyIter," has trigger efficiency = ", 100*(trigNamesAndNumPassing[keyIter]/totalNumEvts) ,"% +/- ", 100*(calcEffUnc(trigNamesAndNumPassing[keyIter],totalNumEvts)), "% in evts without any selection requirements"
	
#end loop over keys in trigNamesAndNumPassing map


