# use this python program to see which HLT paths are fired most often in TTBar->Leptons+Jets evts

import ROOT
import array
import math

# load FWLite C++ libs
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# load FWLite python libs
from DataFormats.FWLite import Handle, Events

# define additional fxns

#save a TTree filled with the kinematic attributes (pt, eta, etc) of the four final state leptons and quarks/jets
#to a file 


#hltPathName is a string
def useThisPath(hltPathName):
	dontSkip = False
	stringsToUse = ["Mu23_TrkIsoVVL_Ele12","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1","HLT_Mu40_TkMu11_v1","HLT_Mu50_v1","HLT_Ele105"]
	for i in xrange(int(len(stringsToUse))):
		if(hltPathName.find(stringsToUse[i]) != -1):
			#if one of the strings which should be used is found in the hlt path name, then break out of this
			#for loop and return True
			dontSkip = True
			break
	return dontSkip

# define handles and labels for TriggerResults and TriggerObjectStandAlone collections
trigResultsHandl, trigResultsLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
trigObjsHandl, trigObjsLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"

# open one or more miniAOD .root files, and create an iterator to loop over the Event objects in the file(s)
doMuonChannel = True 
doTtBar = False 
doDyPlusJets = False  

#by default, use TTBar events
allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTJets_TuneCUETP8M1_13TeV_pythia8_1.root")
if(doTtBar == False and doDyPlusJets == True):#use DY+Jets events
	allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/DYJetsToLL_M-50_TuneCUETP8M1_FlatPU_10_to_50_13TeV_pythia8_1.root")

if(doTtBar == False and doDyPlusJets == False and doMuonChannel == True):#use WR to MuMuJJ events
	allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")


if(doTtBar == False and doDyPlusJets == False and doMuonChannel == False):#use WR to EEJJ events
	allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")


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
	#if(evNum > 4500): break
	totalNumEvts+=1.0
	
	oneEvent.getByLabel(trigObjsLabel, trigObjsHandl)
	oneEvent.getByLabel(trigResultsLabel, trigResultsHandl)
	
	if((evNum%1000)==0): print "evt number: ", evNum
	
	allPathNames = oneEvent.object().triggerNames(trigResultsHandl.product())
	for i in xrange(trigResultsHandl.product().size()):
		#loop over all HLT path names and keep track of the paths which are fired, how many times they fire, and the total number of evts
		#which have been run over
		if(trigResultsHandl.product().accept(i) and useThisPath(allPathNames.triggerName(i)) == True ):
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

datasetId = "ttbar"
if(doTtBar==False and doDyPlusJets==True): datasetId = "dyPlusJets"
if(doTtBar==False and doDyPlusJets==False and doMuonChannel==True): datasetId = "Wr_M800_To_Nu_M400_To_MuMuJJ"
if(doTtBar==False and doDyPlusJets==False and doMuonChannel==False): datasetId = "Wr_M800_To_Nu_M400_To_EEJJ"

#map of cross sxns of datasets considered here (ttbar, dy+jets, WR signal)
#keyed by datasetId, values obtained from cmsWR/doc/datasets.dat
#cross sxn values in femtobarns
crossSxns = dict()
crossSxns["ttbar"]=275900.0
crossSxns["dyPlusJets"]=6104000.0
crossSxns["Wr_M800_To_Nu_M400_To_EEJJ"]=3650.0
crossSxns["Wr_M800_To_Nu_M400_To_MuMuJJ"]=4060.0

#integrated lumi in 1/femtobarns
intgLumi = 1.0

print "total num evts which were scanned over = ", totalNumEvts
for keyIter in iter(trigNamesAndNumPassing):
	for stringIter in iter(crossSxns):
		if(stringIter == datasetId): print "HLT path named: ", keyIter," yields ", ((trigNamesAndNumPassing[keyIter]/totalNumEvts)*crossSxns[stringIter]*intgLumi*intgLumi)," ",stringIter ," evts per femtobarn of data"
	
#end loop over keys in trigNamesAndNumPassing map


