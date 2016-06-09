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
	#stringsToUse = ["Mu23_TrkIsoVVL_Ele12","Mu8_TrkIsoVVL_Ele23","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1","HLT_Ele105_CaloIdVT","HLT_Mu40_TkMu11_v1","HLT_Mu50_v1"]
	#stringsToUse = ["Mu23_TrkIsoVVL_Ele12","Mu8_TrkIsoVVL_Ele23","HLT_DoubleEle33","HLT_Ele105_CaloIdVT","HLT_Mu40_TkMu11_v1","HLT_Mu50_v1","Mu30_Ele30"]
	stringsToUse = ["WP60_SC4_Mass55","WP60_Ele8_Mass55"]
	#stringsToUse = ["tagAndProbe"]
	
	for i in xrange(int(len(stringsToUse))):
		if(hltPathName.find(stringsToUse[i]) != -1):
			#if one of the strings which should be used is found in the hlt path name, then break out of this
			#for loop and return True
			dontSkip = True
			break
	return dontSkip

#return uncertainty on efficiency
#effNumerator and effDenominator are floats
def calcEffUnc(effNumerator, effDenominator):
	uncertainty = (1/effDenominator)*math.sqrt(effNumerator*(1-(effNumerator/effDenominator)))
	return uncertainty
#end calcEffUnc

# define handles and labels for TriggerResults and TriggerObjectStandAlone collections
trigResultsHandl, trigResultsLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
trigObjsHandl, trigObjsLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
wRskimPathHandl, skimPathLabel = Handle("edm::TriggerResults"), ("TriggerResults","","SKIM")
#the reco object handles (jets, eles, muons) are only needed if I want to apply cuts to the reco objects
#before counting the number of evts in which each trigger is fired
recoJetsHandl, recoJetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
recoMuHandl, recoMuLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
recoElectronHandl, recoElectronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"


#use this for same flavor lepton final states
#otherwise the default will be mixed flavor lepton (e+mu) final states
doMixedFlavor = False
doMuonChannel = False
recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
if(doMuonChannel == False): recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"

# open one or more miniAOD .root files, and create an iterator to loop over the Event objects in the file(s)
doTtBar = False 
allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTJets_TuneCUETP8M1_13TeV_pythia8_1.root")
if(doTtBar == False):
	allEvents = Events("root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DoubleEG_RunD_v4_SHv12/DoubleEG_RunD_v4_37_2_vnl.root")
	#allEvents = Events("root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_4_1_Alq.root")
	#allEvents = Events("/eos/uscms/store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronLowMassRegionSkim_1.root")
	#allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/DYJetsToLL_M-50_TuneCUETP8M1_FlatPU_10_to_50_13TeV_pythia8_1.root")



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
	if(evNum > 1000): break
	if(evNum > 80000): break
	
	oneEvent.getByLabel(trigObjsLabel, trigObjsHandl)
	oneEvent.getByLabel(trigResultsLabel, trigResultsHandl)
	oneEvent.getByLabel(wRskimPathHandl, skimPathLabel)
	if((evNum%10000)==0): print "evt number: ", evNum

	#require the tagAndProbe path is passed
	wrPathNames = oneEvent.object().triggerNames(wRskimPathHandl.product())
	for r in xrange(wrPathNames.size()):
		if(wrPathNames[r] == "tagAndProbe"):
			wrPathIndex = r
	if(wRskimPathHandl.product().accept(wrPathIndex) == False): continue

	#find and save the indices of the two HLT paths used for tagandprobe
	totalNumEvts+=1.0
	
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

#save the tree with kinematics info to a file
#fout.cd()
#outputTuple.Write()
#fout.Close()

datasetId = "ttbar"
#if(doTtBar==False): datasetId = "dyPlusJets"
if(doTtBar==False): datasetId = "DoubleEGdata"
print "total num evts before any cuts, from the tagAndProbe path = ", totalNumEvts
#print "total num evts passing GEN cuts = ", totalNumEvtsAfterGenCuts
#print 'total num evts passing gen and reco cuts = ', totalNumEvtsAfterGenRecoMatchingAndRecoCuts
for keyIter in iter(trigNamesAndNumPassing):
	if( (trigNamesAndNumPassing[keyIter]/totalNumEvts) > 0.0):
		print "HLT path named: ", keyIter," has trigger efficiency = ", 100*(trigNamesAndNumPassing[keyIter]/totalNumEvts) ,"% +/- ", 100*(calcEffUnc(trigNamesAndNumPassing[keyIter],totalNumEvts)), "% on ", datasetId, ' events'
	
#end loop over keys in trigNamesAndNumPassing map

#print 'gen and reco requirements are that the four final state particles at GEN lvl from the WR decay are found, and fall within detector eta acceptance, four reco objects are matched to these four GEN particles, and that the reco lepton matched to the Nu daughter lepton has pT > 25 GeV (not applied in muon channel evts)'




