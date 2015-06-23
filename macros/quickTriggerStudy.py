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

#use this fxn to compute deltaR with +/-pi safe deltaPhi fxn from ROOT.TMath
#etaOne, phiOne, etaTwo, and phiTwo are all floats
def deltR(etaOne, phiOne, etaTwo, phiTwo):
	dRval = -1
	dEta = etaOne - etaTwo
	#print 'dEta = ', dEta
	dPhi = ROOT.TVector2.Phi_mpi_pi(phiOne - phiTwo)
	#print 'dPhi = ', dPhi
	dRval = math.sqrt(math.pow(dEta,2) + math.pow(dPhi,2))
	return dRval
#end deltR()


#use this fxn to determine if the evt contains reco objects which are matched to their gen counterparts, within
#a distance dRlept (for leptons) or dRjets (for jets/quarks)
#doMuon and doRecoCuts are either True or False
def foundMatchingRecoObjects(genHndl, recoLeptHndl, recoJetsHndl, doMuon, dRlept, dRjets, doRecoCuts):
	foundRecoMatchesForAllObjects = False
	genList, recoLepts, recoJets = genHndl.product(), recoLeptHndl.product(), recoJetsHndl.product()
	if(recoLepts.size() == 0 or recoJets.size() == 0 or genList.size() == 0): return foundRecoMatchesForAllObjects
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
						#print 'incremented nWrLeptonMatchedRecos by 1'
		if( abs(genList[i].pdgId()) == leptPdg and abs(genList[i].mother(0).pdgId()) == hvyNuPdgId):
			for n in xrange(int(recoLepts.size())):
				slDr = deltR(genList[i].eta(),genList[i].phi(),recoLepts[n].eta(),recoLepts[n].phi())
				#print 'dR btwn reco lepton and gen Nu dau lepton = ', slDr
				if(slDr <= dRlept ):
					if(recoLepts[n].pt() >= recoLeptonMatchedToNuDau.pt() and recoLepts[n].pt() != recoLeptonMatchedToWrDau.pt() and (doRecoCuts==False or recoLepts[n].pt() >= 25) ):
						recoLeptonMatchedToNuDa = recoLepts[n]
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
	if(nWrLeptonMatchedRecos > 0 and nNuLeptonMatchedRecos > 0 and nNuQuarkMatchedRecos > 1): foundRecoMatchesForAllObjects = True

	return foundRecoMatchesForAllObjects
#end foundMatchingRecoObjects()


# define handles and labels for TriggerResults, TriggerObjectStandAlone, and GenParticle collections
trigResultsHandl, trigResultsLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
trigObjsHandl, trigObjsLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
genParticleHandl, genParticleLabel = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"
recoJetsHandl, recoJetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"

# open one or more miniAOD .root files, and create an iterator to loop over the Event objects in the file(s)
doMuonChannel = False 
doRECOCuts = True	#apply cuts on reco objects which are matched to gen particles from WR decay chain
allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")

#dR matching thresholds for leptons and jets
dRforJets = 0.4
dRforLeptons = 0.2

if(doMuonChannel == False):
	allEvents = Events("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root")
	recoLeptonHandl, recoLeptonLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"



#vars and containers which will be filled in loop over events
trigNamesAndNumPassing = dict()
totalNumEvts = 0
totalNumEvtsAfterGenCuts = 0
totalNumEvtsAfterGenRecoMatchingAndRecoCuts = 0

# loop over Event objects in input file(s)
# enumerate returns two-tuple objects (an immutable list where each element in the list is a pair), and thus two different
# iterator variables are needed in the for loop declaration
# the first value of each pair in the two tuple is a simple number index starting at 0, and the second value in each pair
# in the two tuple is the object of interest
for evNum, oneEvent in enumerate(allEvents):
	#evNum points to a simple number, while oneEvent points to an edm::Event object
	#if(evNum > 2): break
	totalNumEvts+=1.0
	
	#link reco::GenParticle, and pat::Electron or pat::Muon collections with a edm::Handle objs
	oneEvent.getByLabel(genParticleLabel, genParticleHandl)
	oneEvent.getByLabel(recoLeptonLabel, recoLeptonHandl)
	oneEvent.getByLabel(recoJetsLabel, recoJetsHandl)
	#check that the correct WR decay products are found at GEN level
	if(foundCorrectWrDecayProducts(genParticleHandl, doMuonChannel) == False): continue
	totalNumEvtsAfterGenCuts+=1.0
	#check that the event has reco objs which are closely matched to the GEN lvl WR decay products
	if(foundMatchingRecoObjects(genParticleHandl, recoLeptonHandl, recoJetsHandl, doMuonChannel,dRforLeptons,dRforJets,doRECOCuts) == False): continue
	totalNumEvtsAfterGenRecoMatchingAndRecoCuts+=1.0
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

print "total num evts before any cuts = ", totalNumEvts
print "total num evts passing GEN cuts = ", totalNumEvtsAfterGenCuts
print 'total num evts passing gen and reco cuts = ', totalNumEvtsAfterGenRecoMatchingAndRecoCuts
for keyIter in iter(trigNamesAndNumPassing):
	if( (trigNamesAndNumPassing[keyIter]/totalNumEvtsAfterGenRecoMatchingAndRecoCuts) > 0.70):
		print "HLT path named: ", keyIter," has trigger efficiency = ", 100*(trigNamesAndNumPassing[keyIter]/totalNumEvtsAfterGenRecoMatchingAndRecoCuts) ,"% +/- ", 100*(calcEffUnc(trigNamesAndNumPassing[keyIter],totalNumEvtsAfterGenRecoMatchingAndRecoCuts)), "% in evts where gen and reco requirements are met"
	
#end loop over keys in trigNamesAndNumPassing map

print 'gen and reco requirements are that the four final state particles at GEN lvl from the WR decay are found, and fall within detector eta acceptance, four reco objects are matched to these four GEN particles, and that the reco lepton matched to the Nu daughter lepton has pT > 25 GeV (not applied in muon channel evts)'




