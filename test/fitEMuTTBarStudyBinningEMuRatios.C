#include "../interface/rooFitFxns.h"


/*
 * use this fxn to loop over events in a TChain, pull the event weights from the TChain,
 * combine the weights with the cross section*luminosity (or some other normalization
 * factor), and store them in a RooDataSet.  This fxn returns a pointer to the RooDataSet.
 * the RooRealVar objects named WR_mass and weight must be in the RooArgSet named dataSetVars
 */
RooDataSet applyNormalization(TChain * chain, TString dataSetName, Float_t normalization, RooArgSet & dataSetVars, RooRealVar & fourObjMass, RooRealVar & evWgt, Bool_t doFxnWgt, Float_t intercept, Float_t slope){
	///declare local vars to read info from the TChain
	Float_t massWR_, weight_;

	///link the local vars to branches in the TChain
	chain->SetBranchAddress("WR_mass",&massWR_);
	chain->SetBranchAddress("weight",&weight_);

	///make an empty RooDataset without identifying the event weight branch, and save a pointer to it
	RooDataSet * tempDataSetPtr = new RooDataSet("temp","temp",dataSetVars);

	Long64_t evts = chain->GetEntries();	///< get number of entries in the TChain

#ifdef DEBUG
	std::cout<< evts <<" events in the TChain"<<std::endl;
#endif
	
	Float_t emuWgt=1.;
	///loop over all entries in the TChain
	for(Long64_t i=0; i<evts; i++){
		chain->GetEntry(i);
		if(massWR_<600.) continue;

		///set the values of the RooRealVar objects using local vars tied to TChain branches
		fourObjMass.setVal(massWR_);
		if(!doFxnWgt) evWgt.setVal(weight_*normalization);
		else{
			emuWgt = intercept + (slope*massWR_);
			evWgt.setVal(weight_*normalization*emuWgt);
		}

		///add the values of the RooRealVar objects to the temporary RooDataSet
		tempDataSetPtr->add(dataSetVars);
	}///end loop over entries in TChain

	tempDataSetPtr->Print();
	chain->ResetBranchAddresses();
	
	///make the RooDataSet from the RooDataSet pointer, and identify the RooRealVar with weights
	RooDataSet dataSetWithWeights(dataSetName, dataSetName, tempDataSetPtr, dataSetVars,"", evWgt.GetName());
	return dataSetWithWeights;

}///end applyNormalization()


void fitEMuTTBarStudyBinningEMuRatios(){
	///////////////////////////////////////////////////////////////////////////////////////////
	//adjust these parameters to determine the number of toys thrown, the number of events generated, and the number of bins in the M_LLJJ distributions which are generated
	int nToys=375;
	//Int_t evtsArr[]={982,2000,4000};	///<there are 982 events in real data from MuonEG in 2015 which pass emu flavor sideband selection
	Int_t evtsArr[]={30000};	///<there are 982 events in real data from MuonEG in 2015 which pass emu flavor sideband selection
	//int nbins[]={10, 20, 30, 40, 50, 70, 90, 110};
	int nbins[]={30};
	Float_t lljjMin = 600.;
	//Float_t lljjMaxArr[]={2000.,2350.,2700.,3000.};
	Float_t lljjMaxArr[]={2000.};
	///////////////////////////////////////////////////////////////////////////////////////////
	vector<Float_t> lljjMaxVect(lljjMaxArr,lljjMaxArr + sizeof(lljjMaxArr)/sizeof(Float_t));
	vector<Int_t> evtsVect(evtsArr,evtsArr + sizeof(evtsArr)/sizeof(Int_t));
	Int_t maxLLJJSize = lljjMaxVect.size();
	Int_t evtsSize = evtsVect.size();
	gStyle->SetOptStat("oRMe");
	gStyle->SetOptFit(1111);
	
	TChain * ttbaremuchain = new TChain("Tree_Iter0");
	ttbaremuchain->Add("/afs/cern.ch/work/s/skalafut/public/WR_starting2015/forPeterRootFiles/selected_tree_TT_flavoursidebandEMu.root");

	//make a RooDataSet from ttbaremu MLLJJ, fit a RooExp to the distribution
	RooRealVar massWR("fourObjectMass","fourObjectMass", 600, 6500);
	RooRealVar evtWeight("evWeightSign", "evWeightSign", -2000000,200000);
	RooArgSet vars(massWR,evtWeight);

	RooDataSet ttbaremuDataset = applyNormalization(ttbaremuchain, "ttbarEMuDataSet",1, vars, massWR, evtWeight, false, 0.0, 0.0);
	Fits::expPower.setVal(-0.004);
	RooFitResult * fitResult = Fits::expPdfRooAbsPdf->fitTo(ttbaremuDataset, RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
	vector<int> nbinsVector(nbins,nbins + sizeof(nbins)/sizeof(int));
	int diffBins = nbinsVector.size();
	RooArgSet testVars(massWR);

	for(Int_t r=0; r<evtsSize; r++){
		for(Int_t q=0; q<maxLLJJSize; q++){
			Float_t lljjMax = lljjMaxVect[q];
			Int_t eventsToGenerate = evtsVect[r];

			TTree * linearFitParamTree = new TTree("t","");
			Float_t intercept=0, slope=0, interceptUncertainty=0, slopeUncertainty=0;
			Int_t numbins=0;
			linearFitParamTree->Branch("intercept",&intercept);
			linearFitParamTree->Branch("interceptUncertainty",&interceptUncertainty);
			linearFitParamTree->Branch("slope",&slope);
			linearFitParamTree->Branch("slopeUncertainty",&slopeUncertainty);
			linearFitParamTree->Branch("numbins",&numbins);
			TFile outputFile(("outputFileWithTTBarTree"+to_string(nToys)+"toys"+"_"+to_string((int) lljjMax)+"maxFitRange_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".root").c_str(), "recreate");

			for(int i=0; i<nToys; i++){

				for(int j=0; j<diffBins; j++){
					intercept=-1, slope=-1, interceptUncertainty=-1, slopeUncertainty=-1, numbins=-1;

					///generate two RooDataSets with a fixed number of events using the exponential fitted to the EMuJJ distribution as a template
					RooDataSet * rooDataSetOne = Fits::expPdfRooAbsPdf->generate(testVars, RooFit::NumEvents(eventsToGenerate));
					RooDataSet * rooDataSetTwo = Fits::expPdfRooAbsPdf->generate(testVars, RooFit::NumEvents(eventsToGenerate));

					///create two histograms from the two RooDataHists
					TH1* binnedHistOne = rooDataSetOne->createHistogram("massLLJJOne", massWR, RooFit::Binning(nbinsVector[j], lljjMin, lljjMax));
					binnedHistOne->Sumw2(true);
					TH1* binnedHistTwo = rooDataSetTwo->createHistogram("massLLJJTwo", massWR, RooFit::Binning(nbinsVector[j], lljjMin, lljjMax));
					binnedHistTwo->Sumw2(true);
		
					TH1* ratioHistOneTwo = (TH1*) binnedHistOne->Clone();
					ratioHistOneTwo->Sumw2(true);
					ratioHistOneTwo->Divide(binnedHistTwo);
					TF1* ratioFit = new TF1("ratioFit","[0]*x+[1]", lljjMin, lljjMax);
					ratioHistOneTwo->Fit("ratioFit","Q");

					if((i==0 || i==50 || i==100 || i==150 || i==200 || i==250 || i==300) && (j==0 || j==2 || j==4) ){
						TCanvas* canv = new TCanvas("canv","canv",800,800);
						canv->cd();

						binnedHistOne->Draw();
						canv->SaveAs(("binnedHistOne_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".pdf").c_str(),"recreate");
						canv->SaveAs(("binnedHistOne_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".png").c_str(),"recreate");
						canv->Clear();
						binnedHistTwo->Draw();
						canv->SaveAs(("binnedHistTwo_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".pdf").c_str(),"recreate");
						canv->SaveAs(("binnedHistTwo_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".png").c_str(),"recreate");

						canv->Clear();
						ratioHistOneTwo->Draw();
						ratioFit->SetLineColor(kRed);
						ratioFit->Draw("same");
						canv->SaveAs(("ratioHistoWithFit_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".pdf").c_str(), "recreate");
						canv->SaveAs(("ratioHistoWithFit_nbins_"+to_string(nbinsVector[j])+"_toy_"+to_string(i)+"_maxMLLJJ_"+to_string((int) lljjMax)+"_"+to_string((int) eventsToGenerate)+"eventsGenerated"+".png").c_str(), "recreate");
						delete canv;
					}

					numbins=nbinsVector[j];
					intercept=ratioFit->GetParameter(1), slope=ratioFit->GetParameter(0), interceptUncertainty=ratioFit->GetParError(1), slopeUncertainty=ratioFit->GetParError(0);
					linearFitParamTree->Fill();

					///delete the two RooDataHists
					delete ratioHistOneTwo;
					delete ratioFit;
					delete binnedHistOne;
					delete binnedHistTwo;
					delete rooDataSetOne;
					delete rooDataSetTwo;
				}///end loop over different bin sizes

			}///end loop over toys

			linearFitParamTree->Write();
		}///end loop over different fit ranges
	}///end loop over different event counts to be generated


}///end fitEMuTTBarStudyBinningEMuRatios()

