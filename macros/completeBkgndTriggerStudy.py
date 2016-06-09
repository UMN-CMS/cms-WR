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

#check that the correct WR decay products are found at GEN lvl
#genHandle is an edm::Handle object linked to the prunedGenParticles collection
#doMuon is True or False to indicate the final state of interest (True --> mumuqq, False --> eeqq)
def foundCorrectWrDecayProducts(genHandle, doMuon):
	#print 'in foundCorrectWrDecayProducts fxn'
	foundAllDecayProducts = False
	maxEtaQrk, maxEtaLept = 2.5, 2.5 
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

#use this fxn to determine if the evt contains reco objects which are matched to their gen counterparts, within
#a distance dRlept (for leptons) or dRjets (for jets/quarks)
#doMuon and doRecoCuts are either True or False
def foundMatchingRecoObjects(genHndl, recoLeptHndl, recoJetsHndl, doMuon, dRlept, dRjets, doRecoCuts):
	objsAndFoundMatchesList = []
	foundRecoMatchesForAllObjects = False
	genList, recoLepts, recoJets = genHndl.product(), recoLeptHndl.product(), recoJetsHndl.product()
	#if(recoLepts.size() == 0 or recoJets.size() == 0 or genList.size() == 0): return foundRecoMatchesForAllObjects
	if(recoLepts.size() == 0 or recoJets.size() == 0 or genList.size() == 0): return objsAndFoundMatchesList 
	recoLeptonMatchedToWrDau, recoLeptonMatchedToNuDau = recoLepts[recoLepts.size()-1], recoLepts[recoLepts.size()-1]
	recoJetMatchedToLeadNuDauQrk, recoJetMatchedToSubleadNuDauQrk = recoJets[recoJets.size()-1], recoJets[recoJets.size()-1]

	leptPdg, hvyNuPdgId, wRpdgId = 99, 100, 9900024
	if(doMuon == True): leptPdg, hvyNuPdgId = 13, 9900014
	else: leptPdg, hvyNuPdgId = 11, 9900012

	nWrLeptonMatchedRecos, nNuLeptonMatchedRecos, nNuQuarkMatchedRecos = 0,0,0
	#print 'there are ', genList.size(), ' gen particles in the evt'
	for i in xrange(int(genList.size())):
		#loop over all gen particles and find if there are any matching reco objects
		#if doRecoCuts == True, then add the reco objects to appropriate lists
		if( abs(genList[i].pdgId()) == leptPdg and abs(genList[i].mother(0).pdgId()) == wRpdgId):
			for l in xrange(int(recoLepts.size())):
				llDr = deltR(genList[i].eta(),genList[i].phi(),recoLepts[l].eta(),recoLepts[l].phi())
				#print 'dR btwn reco lepton and gen WR dau lepton = ', llDr
				if(llDr <= dRlept ):
					#need >= in case the initialization of recoLeptonMatchedToWrDau happens to be the best match
					if(recoLepts[l].pt() >= recoLeptonMatchedToWrDau.pt()):
						recoLeptonMatchedToWrDau = recoLepts[l]
						nWrLeptonMatchedRecos+=1
						#print 'pt of reco lepton matched to Wr dau = ', recoLeptonMatchedToWrDau.pt()
						#print 'incremented nWrLeptonMatchedRecos by 1'
						#print ' '
		if( abs(genList[i].pdgId()) == leptPdg and abs(genList[i].mother(0).pdgId()) == hvyNuPdgId):
			for n in xrange(int(recoLepts.size())):
				slDr = deltR(genList[i].eta(),genList[i].phi(),recoLepts[n].eta(),recoLepts[n].phi())
				#print 'dR btwn reco lepton and gen Nu dau lepton = ', slDr
				if(slDr <= dRlept ):
					if(recoLepts[n].pt() >= recoLeptonMatchedToNuDau.pt() and recoLepts[n].pt() != recoLeptonMatchedToWrDau.pt() and (doRecoCuts==False or recoLepts[n].pt() >= 25) ):
						#print 'reco cuts are applied: ', doRecoCuts
						recoLeptonMatchedToNuDau = recoLepts[n]
						#print 'pt of reco lepton matched to Nu dau = ', recoLeptonMatchedToNuDau.pt()
						nNuLeptonMatchedRecos+=1
						#print 'incremented nNuLeptonMatchedRecos by 1'
		if( abs(genList[i].pdgId()) <= 6 and abs(genList[i].mother(0).pdgId()) == hvyNuPdgId):
			#loop over the jets and see if there are any which are geometrically close (eta, phi) to the gen quarks
			for j in xrange(int(recoJets.size())):
				jetDr = deltR(genList[i].eta(),genList[i].phi(),recoJets[j].eta(),recoJets[j].phi())
				#print 'dR btwn reco jet and gen Nu dau quark = ', jetDr
				if(jetDr <= dRjets):
					ptLeadJet, ptSubleadJet = recoJetMatchedToLeadNuDauQrk.pt(), recoJetMatchedToSubleadNuDauQrk.pt()
					#need >= for the same reason as is listed above for the lepton dau coming directly from the WR decay
					if(recoJets[j].pt() >= ptLeadJet):
						#update the reco jets matched to the leading and subleading gen quarks
						recoJetMatchedToSubleadNuDauQrk = recoJetMatchedToLeadNuDauQrk
						recoJetMatchedToLeadNuDauQrk = recoJets[j]
						nNuQuarkMatchedRecos+=1
						#print 'incremented nNuQuarkMatchedRecos by 1'
					if(recoJets[j].pt() < ptLeadJet and nNuQuarkMatchedRecos == 0):
						#update the subleading jet object
						recoJetMatchedToLeadNuDauQrk = recoJets[j]
						nNuQuarkMatchedRecos+=1
						#print 'incremented nNuQuarkMatchedRecos by 1'
					if(recoJets[j].pt() < ptLeadJet and recoJets[j].pt() >= ptSubleadJet):
						#update the subleading jet object
						recoJetMatchedToLeadNuDauQrk = recoJets[j]
						nNuQuarkMatchedRecos+=1
						#print 'incremented nNuQuarkMatchedRecos by 1'
	#end loop over gen particles in event
	#print 'after looping over all gen particles'
	#print 'there were ', nWrLeptonMatchedRecos,' reco leptons matched to the WR gen lepton daughter'
	#print 'there were ', nNuLeptonMatchedRecos,' reco leptons matched to the Nu gen lepton daughter'
	#print 'there were ', nNuQuarkMatchedRecos,' reco jets matched to the two Nu gen quark daughters'


	#if the number of matched reco objects is sufficient, then update foundRecoMatchesForAllObjects to True
	#and add the four reco objects to the array which is returned
	if(nWrLeptonMatchedRecos > 0 and nNuLeptonMatchedRecos > 0 and nNuQuarkMatchedRecos > 1):
		foundRecoMatchesForAllObjects = True
		objsAndFoundMatchesList.append(1)
		objsAndFoundMatchesList.append(recoLeptonMatchedToWrDau)
		objsAndFoundMatchesList.append(recoLeptonMatchedToNuDau)
		objsAndFoundMatchesList.append(recoJetMatchedToLeadNuDauQrk)
		objsAndFoundMatchesList.append(recoJetMatchedToSubleadNuDauQrk)
	
	if(foundRecoMatchesForAllObjects == False):
		objsAndFoundMatchesList.append(0)
	
	return objsAndFoundMatchesList
#end foundMatchingRecoObjects()


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
#totalNumEvtsAfterGenCuts = 0
#totalNumEvtsAfterGenRecoMatchingAndRecoCuts = 0

#tree with kinematics of final state reco objects, and variables tied to branches
#foutName = ""
#if(doTtBar == False): foutName = "recoElectronKinematicsTree.root"
#if(doTtBar == True): foutName = "recoMuonKinematicsTree.root"
#fout = ROOT.TFile(foutName,"RECREATE")
#fout.cd()
#
##names of the variables in output tuple
#varNames = ["WrDauLeptPt","WrDauLeptEta","WrDauLeptPhi","NuDauLeptPt","NuDauLeptEta","NuDauLeptPhi","leadJetPt","leadJetEta","leadJetPhi","subleadJetPt","subleadJetEta","subleadJetPhi","evNum"]
#
##create the ntuple
#outputTuple = ROOT.TNtuple("recoKinematics","tuple",":".join(varNames))



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
	
	#oneEvent.getByLabel(recoJetsLabel, recoJetsHandl)
	##split workflow between mixed flavor and same flavor lepton final states
	#if(doMixedFlavor==True):
	#	#interested in emu final state
	#	oneEvent.getByLabel(recoMuLabel, recoMuHandl)
	#	oneEvent.getByLabel(recoElectronLabel, recoElectronHandl)
	#	asd


	#if(doMixedFlavor==False):
	#	#interested in ee or mumu final states
	#	oneEvent.getByLabel(recoLeptonLabel, recoLeptonHandl)
	
	#check that the correct WR decay products are found at GEN level
	#if(foundCorrectWrDecayProducts(genParticleHandl, doTtBar) == False): continue
	#totalNumEvtsAfterGenCuts+=1.0
	##check that the event has reco objs which are closely matched to the GEN lvl WR decay products
	##objsAndFoundMatchesArr = []
	#objsAndFoundMatchesArr = foundMatchingRecoObjects(genParticleHandl, recoLeptonHandl, recoJetsHandl, doTtBar,dRforLeptons,dRforJets,doRECOCuts)
	#if(int(len(objsAndFoundMatchesArr)) < 2): continue
	##if the execution reaches this point then there are many entries in objsAndFoundMatchesArr, corresponding to
	##the pt, eta, and phi of each final state particle (lepton, jet)
	#totalNumEvtsAfterGenRecoMatchingAndRecoCuts+=1.0

	#outputVals = []
	##calculate the output values and put them into a list
	##there are only four objects of interest, and we need the pt, eta, and phi for each one
	#for i in xrange(4):
	#	outputVals.append(objsAndFoundMatchesArr[i+1].pt())
	#	outputVals.append(objsAndFoundMatchesArr[i+1].eta())
	#	outputVals.append(objsAndFoundMatchesArr[i+1].phi())
	##end loop over reco objects

	##add an event number to the tuple
	#outputVals.append(evNum)

	##add the values in the list to the TNtuple object named outputTuple
	#outputTuple.Fill(array.array("f",outputVals))
	
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




