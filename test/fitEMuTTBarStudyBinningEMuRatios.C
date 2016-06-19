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
	TChain * ttbaremuchain = new TChain("Tree_Iter0");
	ttbaremuchain->Add("/afs/cern.ch/work/s/skalafut/public/WR_starting2015/forPeterRootFiles/selected_tree_TT_flavoursidebandEMu.root");

	//make a RooDataSet from ttbaremu MLLJJ, fit a RooExp to the distribution
	RooRealVar massWR("fourObjectMass","fourObjectMass", 600, 6500);
	RooRealVar evtWeight("evWeightSign", "evWeightSign", -2000000,200000);
	RooArgSet vars(massWR,evtWeight);
	
	RooDataSet ttbaremuDataset = applyNormalization(ttbaremuchain, "ttbarEMuDataSet",1, vars, massWR, evtWeight, false, 0.0, 0.0);
	Fits::expPower.setVal(-0.004);
	RooFitResult * fitResult = Fits::expPdfRooAbsPdf->fitTo(ttbaremuDataset, RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
	int nToys=1;
	Int_t eventsToGenerate=4000;
	int nbins[]={10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	vector<int> nbinsVector(nbins,nbins + sizeof(nbins)/sizeof(int));
	int diffBins = nbinsVector.size();
	Float_t lljjMin = 600., lljjMax = 2000.;
	RooArgSet testVars(massWR);
	map<int,Double_t> mapNBinsToIntercept;
	map<int,Double_t> mapNBinsToSlope;
	Double_t tempIntercepts=0, tempSlopes=0;
	for(int j=0; j<diffBins; j++){
		///generate two RooDataHists with a fixed number of events using the exponential fitted to the EMuJJ distribution as a template
		RooDataHist * rooHistOne = Fits::expPdfRooAbsPdf->generateBinned(testVars, RooFit::NumEvents(eventsToGenerate));
		RooDataHist * rooHistTwo = Fits::expPdfRooAbsPdf->generateBinned(testVars, RooFit::NumEvents(eventsToGenerate));

		//throw many toys for each value of nbinsVector
		for(int i=0; i<nToys; i++){
			///create two histograms from the two RooDataHists
			TH1* binnedHistOne = rooHistOne->createHistogram("massLLJJOne", massWR, RooFit::Binning(nbinsVector[j], lljjMin, lljjMax));
			TH1* binnedHistTwo = rooHistTwo->createHistogram("massLLJJTwo", massWR, RooFit::Binning(nbinsVector[j], lljjMin, lljjMax));
			TH1* ratioHistOneTwo = (TH1*) binnedHistOne->Clone();
			ratioHistOneTwo->Divide(binnedHistTwo);
			TF1* ratioFit = new TF1("ratioFit","[0]*x+[1]", lljjMin, lljjMax);
			ratioHistOneTwo->Fit("ratioFit");
			
			if(i==0 && (j==0 || j==2) ){
				TCanvas* canv = new TCanvas("canv","canv",800,800);
				canv->cd();
				ratioHistOneTwo->Draw();
				ratioFit->SetLineColor(kRed);
				ratioFit->Draw("same");
				canv->SaveAs(("ratioHistoWithFit_nbins_"+to_string(j)+".pdf").c_str(), "recreate");
				canv->SaveAs(("ratioHistoWithFit_nbins_"+to_string(j)+".png").c_str(), "recreate");
				delete canv;
			}
			tempIntercepts += ratioFit->GetParameter(1);
			tempSlopes += ratioFit->GetParameter(0);

			///delete the two RooDataHists
			delete ratioHistOneTwo;
			delete ratioFit;
			delete binnedHistOne;
			delete binnedHistTwo;
		}///end loop over toys
		mapNBinsToIntercept[nbinsVector[j]] = (Double_t) (tempIntercepts/nToys);
		mapNBinsToSlope[nbinsVector[j]] = (Double_t) (tempSlopes/nToys);
		
		tempIntercepts=0, tempSlopes=0;
		delete rooHistOne;
		delete rooHistTwo;
	}///end loop over different bin sizes
	for(int k=0; k<diffBins; k++){
		std::cout<<"when there are\t"<< nbinsVector[k] <<"\tbins, the average intercept and slope of the ratio plot are\t"<< mapNBinsToIntercept[nbinsVector[k]] <<"\tand\t"<< mapNBinsToSlope[nbinsVector[k]] << std::endl;
	}//end loop over elements in maps
	
}///end fitEMuTTBarStudyBinningEMuRatios()

