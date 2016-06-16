//#include "../interface/rooFitFxns.h"
#define DOTTBAR	///comment this line out to run over DY datasets

/*
//if a file containing normalized pdf integral values (like from 2.0 to 3.0 TeV) already exists, then delete its contents
void clearExistingPdfIntegralTxtFile(std::string pdfIntgrlFileName){
	std::fstream checkFileExistence(.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);
	if(){

	}

}
*/

void multiStepFitTo(RooDataSet dataSetToFit, RooAbsPdf * pdfToFit, RooFitResult * fitRsltPtr, RooRealVar observable){
	RooAbsReal * nll = pdfToFit->createNLL(dataSetToFit);
	RooMinuit m(*nll);
	m.setVerbose(kTRUE);
	m.hesse();	///<run HESSE to estimate step size
	m.migrad();	///<run migrad to do minimization
	pdfToFit->getParameters(observable)->Print("s");
	std::cout<<" "<<std::endl;
	std::cout<<"/////////////// HESSE errors on fit parameters ///////////////"<<std::endl;
	std::cout<<"/////////////// HESSE errors on fit parameters ///////////////"<<std::endl;
	std::cout<<"/////////////// HESSE errors on fit parameters ///////////////"<<std::endl;
	m.hesse();	///<run HESSE to calculate errors of fit params
	pdfToFit->getParameters(observable)->Print("s");
	
	std::cout<<" "<<std::endl;
	std::cout<<" "<<std::endl;
	std::cout<<"/////////////// MINOS errors on fit parameters ///////////////"<<std::endl;
	std::cout<<"/////////////// MINOS errors on fit parameters ///////////////"<<std::endl;
	std::cout<<"/////////////// MINOS errors on fit parameters ///////////////"<<std::endl;
	m.minos();	///<run MINOS to calculate asymmetric errors of fit params
	pdfToFit->getParameters(observable)->Print("s");
	fitRsltPtr = m.save();


}///end multiStepFitTo()

void saveIntegralToFile(std::string fitParamTxtFile, Float_t intgrl, std::string pdfDetails){
	std::ofstream fitParamStream(fitParamTxtFile.c_str(), std::ofstream::app);
	fitParamStream << pdfDetails << "\t" << intgrl << std::endl;
	fitParamStream.close();
}///end saveIntegralToFile()


///add a std::pair which contains a RooAbsPdf pointer, and a RooFitResult pointer to a RooWorkspace
void saveArgSetToRooWorkspace(RooArgSet & rooArgSet, RooWorkspace & rooWS, Bool_t writeFile, TString wsFileName){
	rooWS.import(rooArgSet);
	if(writeFile) rooWS.writeToFile(wsFileName);
}

void saveDataSetToRooWorkspace(RooDataSet & rooDST, RooWorkspace & rooWS, Bool_t writeFile, TString wsFileName){
	rooWS.import(rooDST);
	if(writeFile) rooWS.writeToFile(wsFileName);
}

void saveRooPlot(RooPlot * rFrame, TString cnm, TString plotFileName){
	TCanvas * con = new TCanvas(cnm, cnm, 800,800);
	con->cd();
	rFrame->Draw();
	con->Update();
	con->SaveAs(plotFileName+".png","recreate");
	con->SaveAs(plotFileName+".pdf","recreate");
	con->Close();
	con->Delete();
}///end saveRooPlot()

/*
 * use this fxn to loop over events in a TChain, pull the event weights from the TChain,
 * combine the weights with the cross section*luminosity (or some other normalization
 * factor), and store them in a RooDataSet.  This fxn returns a pointer to the RooDataSet.
 * the RooRealVar objects named fourObjMass and evWgt must be in the RooArgSet named dataSetVars
 */
RooDataSet applyNormalization(TChain * chain, TString dataSetName, Float_t normalization, RooArgSet & dataSetVars, RooRealVar & fourObjMass, RooRealVar & evWgt, Bool_t doFxnWgt, Float_t intercept, Float_t slope){
	///declare local vars to read info from the TChain
	Float_t massWR_, weight_;

	///link the local vars to branches in the TChain
	chain->SetBranchAddress("fourObjectMass",&massWR_);
	chain->SetBranchAddress("evWeightSign",&weight_);

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

void fitPdfReturnPdfAndResults(RooFitResult * ftRslt, RooDataSet & dataSetInForPlotting, RooDataSet & dataSetInForFitting, RooAbsPdf * pdf, Bool_t savePlots, TString plotTitle, TString outputFileName, TString cvname, TString xAxisTitle, RooRealVar & independentPlotVar, Int_t nVisualizationBins, Double_t & intgrlChk, Float_t normFactor, RooArgSet & pdfArgSet, RooArgSet & pdfIntegralArgSet, Bool_t setPdfNorm, std::string integralFileName, std::string fitPdfDetails){
	///fit the pdf to the appropriate dataset
	ftRslt = pdf->fitTo(dataSetInForFitting, RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
	//multiStepFitTo(dataSetInForFitting, pdf, ftRslt, independentPlotVar);
	
	//make a RooPlot to show the RooDataSet and fit pdf
	RooPlot * frm = independentPlotVar.frame(nVisualizationBins);
	frm->GetXaxis()->SetTitle(xAxisTitle);
	frm->SetTitle(plotTitle);
	//plot the pdf and RooDataSet on the RooPlot, then make the residual and pull distributions
	//the RooDataSet for plotting will usually be different from the RooDataSet used for fitting in TTBar estimations
	dataSetInForPlotting.plotOn(frm);
	if(setPdfNorm) pdf->plotOn(frm,RooFit::LineColor(kRed), RooFit::Normalization(normFactor,-1));
	else pdf->plotOn(frm,RooFit::LineColor(kRed));
	pdfArgSet.add(*pdf);
	
	RooArgSet *vset = pdf->getVariables();	
	RooRealVar *var_pdf = (RooRealVar*)vset->find("fourObjectMass");
	var_pdf->setRange("range_f",2000,3000);

	RooAbsReal * integral = pdf->createIntegral(RooArgSet(*var_pdf),RooFit::Range("range_f"));
	Float_t intVl = (Float_t) integral->getVal();
	pdfIntegralArgSet.add(*integral);
	intgrlChk = integral->getVal();
	cout<<"pdf named "<< pdf->GetName() << " has INTEGRAL=\t"<<integral->getVal()<<endl;

	saveIntegralToFile(integralFileName, intVl, fitPdfDetails);

	if(savePlots){
		RooHist * resid = frm->residHist();
		RooHist * pull = frm->pullHist();
		
		RooPlot * residFrm = independentPlotVar.frame(nVisualizationBins);
		residFrm->GetXaxis()->SetTitle(xAxisTitle);
		residFrm->SetTitle("Residual "+plotTitle);
		RooPlot * pullFrm = independentPlotVar.frame(nVisualizationBins);
		pullFrm->GetXaxis()->SetTitle(xAxisTitle);
		pullFrm->SetTitle("Residual/Error "+plotTitle);
		residFrm->addPlotable(resid,"P");
		pullFrm->addPlotable(pull,"P");
	
		saveRooPlot(frm, cvname, outputFileName);
		saveRooPlot(residFrm,cvname+"1", outputFileName+"_resid");
		saveRooPlot(pullFrm,cvname+"2", outputFileName+"_pull");

		residFrm->Delete();
		pullFrm->Delete();
	}
	
	frm->Delete();
	
}///end fitPdfReturnPdfAndResults()

void commonRooFitBkgndEstimation(){

	///declare TChains to real data and MC
	float xsecs[] = {19.32,2.731,0.241,0.01678,0.00139,0.00008948,0.000004135,4.56e-7,2.066e-8};
	TString DY_MuMu_names[] = {"DY_MuMuDataSet_120to200","DY_MuMuDataSet_200to400","DY_MuMuDataSet_400to800","DY_MuMuDataSet_800to1400","DY_MuMuDataSet_1400to2300","DY_MuMuDataSet_2300to3500","DY_MuMuDataSet_3500to4500","DY_MuMuDataSet_4500to6000","DY_MuMuDataSet_6000toInf"};
	std::vector<TChain*> DY_MuMuChain;
#ifndef DOTTBAR
	TChain * DY_120to200Chain = new TChain("120to200");
	DY_120to200Chain->Add("DY_Selected_MuMu120to200.root");
	DY_MuMuChain.push_back(DY_120to200Chain);
	TChain * DY_200to400Chain = new TChain("200to400");
	DY_200to400Chain->Add("DY_Selected_MuMu200to400.root");
	DY_MuMuChain.push_back(DY_200to400Chain);
	TChain * DY_400to800Chain = new TChain("400to800");
	DY_400to800Chain->Add("DY_Selected_MuMu400to800.root");
	DY_MuMuChain.push_back(DY_400to800Chain);
	TChain * DY_800to1400Chain = new TChain("800to140");
	DY_800to1400Chain->Add("DY_Selected_MuMu800to1400.root");
	DY_MuMuChain.push_back(DY_800to1400Chain);
	TChain * DY_1400to2300Chain = new TChain("1400to23");
	DY_1400to2300Chain->Add("DY_Selected_MuMu1400to2300.root");
	DY_MuMuChain.push_back(DY_1400to2300Chain);
	TChain * DY_2300to3500Chain = new TChain("2300to35");
	DY_2300to3500Chain->Add("DY_Selected_MuMu2300to3500.root");
	DY_MuMuChain.push_back(DY_2300to3500Chain);
	TChain * DY_3500to4500Chain = new TChain("3500to45");
	DY_3500to4500Chain->Add("DY_Selected_MuMu3500to4500.root");
	DY_MuMuChain.push_back(DY_3500to4500Chain);
	TChain * DY_4500to6000Chain = new TChain("4500to60");
	DY_4500to6000Chain->Add("DY_Selected_MuMu4500to6000.root");
	DY_MuMuChain.push_back(DY_4500to6000Chain);
	TChain * DY_6000toInfChain = new TChain("6000toIn");
	DY_6000toInfChain->Add("DY_Selected_MuMu6000toInf.root");
	DY_MuMuChain.push_back(DY_6000toInfChain);
#endif

#ifdef DOTTBAR
	TChain * eMuMCChain = new TChain("recoAnalyzerOne/recoObjectsAllCuts","");
	eMuMCChain->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_emujj_signal_region/analyzed_tree_TTJetsDiLeptPartTwo_Madgraph_ToEMuJJSignalRegion_HLTMu30Ele30_emujj.root");
	TChain * eEMCChain = new TChain("unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts","");
	eEMCChain->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/analyzed_TTJetsDiLept_Madgraph_25ns_eejj_signal_region.root");
	TChain * muMuMCChain = new TChain("unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts","");
	muMuMCChain->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_mumujj_signal_region/analyzed_TTJetsDiLept_Madgraph_25ns_mumujj_signal_region.root");
	TChain * dataChain = new TChain("recoAnalyzerOne/recoObjectsAllCuts","");
	dataChain->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_mumujj_signal_region/analyzed_tree_MuonEGRun2015CandDEMuJJSignalRegion_HLTMu30Ele30_allSilver.root");
#endif

	RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,6500);
	RooRealVar evtWeight("evWeightSign", "evWeightSign", -2,2);
	RooArgSet vars(massWR,evtWeight);
	
#ifdef DOTTBAR
	RooDataSet eMuMCDataSet = applyNormalization(eMuMCChain,"eMuMCDataSet",intLumi*ttBarXsxn/nTTBarEvts, vars, massWR, evtWeight, true, emuInter, 0.);
	RooDataSet eMuRealDataSet = applyNormalization(dataChain,"eMuRealDataSet",1, vars, massWR, evtWeight, false, 0., 0.);
	RooDataSet eEMCDataSet = applyNormalization(eEMCChain,"eEMCDataSet",intLumi*ttBarXsxn/nTTBarEvts, vars, massWR, evtWeight, true, emuInter, 0.);
	RooDataSet muMuMCDataSet = applyNormalization(muMuMCChain,"muMuMCDataSet",intLumi*ttBarXsxn/nTTBarEvts, vars, massWR, evtWeight, true, emuInter, 0.);
#endif

	/// Drell-Yan ///	
	std::vector<RooDataSet> DY_MuMuDataSet;
#ifndef DOTTBAR	
	for(int i=0;i<9;i++){
	  RooDataSet DY = applyNormalization(DY_MuMuChain[i],DY_MuMu_names[i],intLumi*xsecs[i]/100000., vars, massWR, evtWeight, false, 0, 0.);
	  DY_MuMuDataSet.push_back(DY);
	}

	for(int i=1;i<9;i++){
	  DY_MuMuDataSet[0].append(DY_MuMuDataSet[i]);
	}
#endif
	
	RooFitResult * tempFitRslt = NULL;
	Double_t integralCheck = 0.;
	RooArgSet elePdfs("elePdfsArgSet");
	RooArgSet eleIntegralChecks("eleIntegralChecksArgSet");
	
	///set reasonable initial values for fit params
	pwrOne.setVal(25.), pwrTwo.setVal(2.), modExpPow.setVal(0.9), expPower.setVal(-0.003), multExpPowFactor.setVal(0.00001), expPowerForQuadExp.setVal(-0.003), expPowerForModExp.setVal(-0.1);
	

#ifdef DOTTBAR
	TString fileNameAppend = "_multiStepFit";
	//electron channel
	std::string ttBarElePdfIntegralFile = "ttBarElePdfIntegrals.txt";
	fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, expPdfRooAbsPdf, false, "TTBar EEJJ Estimate with Exponential fit to MuonEG data","MEEJJ_ttBar_data_driven_estimate_exponentialFit"+fileNameAppend,"c1","M_EEJJ [GeV]", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true, ttBarElePdfIntegralFile,"ttBarEleExponentialFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, pwrLawRooAbsPdf, false, "TTBar EEJJ Estimate with Power Law fit to MuonEG data","MEEJJ_ttBar_data_driven_estimate_powerLawFit"+fileNameAppend,"c2","M_EEJJ [GeV]", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true, ttBarElePdfIntegralFile,"ttBarElePowerLawFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, modExpRooAbsPdf, false, "TTBar EEJJ Estimate with Mod Exponential fit to MuonEG data","MEEJJ_ttBar_data_driven_estimate_modExponentialFit"+fileNameAppend,"c3","M_EEJJ [GeV]", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true, ttBarElePdfIntegralFile,"ttBarEleModExponentialFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, quadExpRooAbsPdf, false, "TTBar EEJJ Estimate with Quadratic Exponential fit to MuonEG data","MEEJJ_ttBar_data_driven_estimate_quadExponentialFit"+fileNameAppend,"c4","M_EEJJ [GeV]", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true, ttBarElePdfIntegralFile,"ttBarEleQuadExponentialFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, linearExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true);

	//fitPdfReturnPdfAndResults(tempFitRslt, eEMCDataSet, eMuRealDataSet, sumTwoExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcEEtoEMuRatioIntercept*nEMuDataEvts, elePdfs, eleIntegralChecks, true);

	
	//make, fill, and save the RooWorkspace
	RooWorkspace electronEstimate("ttBarElectronEstimate");
	saveArgSetToRooWorkspace(elePdfs, electronEstimate, false, "");
	std::cout<<" "<<std::endl;
	std::cout<<"/////////// making a RooArgSet in the RooWorkspace by calling defineSet() ////////"<<std::endl;
	//TString pdfNameList[] = {"rescaledExpPdf","rescaledPowerLawSystematicPdf","modExp","quadraticExp"};
	electronEstimate.defineSet("elePdfsRooArgSet", electronEstimate.allPdfs(), kTRUE);
	std::cout<<"/////////// made a RooArgSet in the RooWorkspace by calling defineSet() ////////"<<std::endl;
	std::cout<<" "<<std::endl;
	//saveArgSetToRooWorkspace(eleIntegralChecks, electronEstimate, false, "");
	saveDataSetToRooWorkspace(eMuRealDataSet, electronEstimate, false,"");
	saveDataSetToRooWorkspace(eEMCDataSet, electronEstimate, true,"ttBarBkgndEleEstimate.root");


	/*
	//redo the same steps, but use the mumu/emu interpolation factor and muMuMCDataSet
	//muon channel
	RooArgSet muonPdfs("muonPdfsArgSet");
	RooArgSet muonIntegralChecks("muonIntegralChecksArgSet");
	std::string ttBarMuonPdfIntegralFile = "ttBarMuonPdfIntegrals.txt";
	
	fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, expPdfRooAbsPdf, true, "TTBar MuMuJJ Estimate with Exponential fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_exponentialFit"+fileNameAppend,"m1","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true, ttBarMuonPdfIntegralFile,"ttBarMuonExponentialFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, pwrLawRooAbsPdf, true, "TTBar MuMuJJ Estimate with Power Law fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_powerLawFit"+fileNameAppend,"m2","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true, ttBarMuonPdfIntegralFile,"ttBarMuonPowerLawFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, modExpRooAbsPdf, true, "TTBar MuMuJJ Estimate with Mod Exponential fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_modExponentialFit"+fileNameAppend,"m3","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true, ttBarMuonPdfIntegralFile,"ttBarMuonModExponentialFit integral from 2.0 to 3.0 TeV = ");

	fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, quadExpRooAbsPdf, true, "TTBar MuMuJJ Estimate with Quadratic Exponential fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_quadExponentialFit"+fileNameAppend,"m4","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true, ttBarMuonPdfIntegralFile,"ttBarMuonQuadExponentialFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, linearExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true);

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, sumTwoExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, true);

	
	RooWorkspace muonEstimate("ttBarMuonEstimate");
	saveArgSetToRooWorkspace(muonPdfs, muonEstimate, false, "");
	//saveArgSetToRooWorkspace(muonIntegralChecks, muonEstimate, false, "");
	saveDataSetToRooWorkspace(eMuRealDataSet, muonEstimate, false,"");
	saveDataSetToRooWorkspace(muMuMCDataSet, muonEstimate, true,"ttBarBkgndMuonEstimate.root");
	*/

#endif

#ifndef DOTTBAR
	//DY background fits
	//muon channel
	RooArgSet DYmuonPdfs("DYmuonPdfsArgSet");
	RooArgSet DYmuonIntegralChecks("DYmuonIntegralChecksArgSet");
	std::string DYMuonPdfIntegralFile = "DYMuonPdfIntegrals.txt";
	
	fitPdfReturnPdfAndResults(tempFitRslt, DY_MuMuDataSet[0], DY_MuMuDataSet[0], expPdfRooAbsPdf, true, "DY MuMuJJ Exponential fit","DY_MMuMuJJ_exponentialFit","m1","M_MuMuJJ [GeV]", massWR, 100, integralCheck, 1, DYmuonPdfs, DYmuonIntegralChecks, false, DYMuonPdfIntegralFile, "DYMuonExponentialFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, pwrLawRooAbsPdf, true, "TTBar MuMuJJ Estimate with Power Law fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_powerLawFit","m2","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, false, DYMuonPdfIntegralFile, "DYMuonPowerLawFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, modExpRooAbsPdf, true, "TTBar MuMuJJ Estimate with Mod Exponential fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_modExponentialFit","m3","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, false, DYMuonPdfIntegralFile, "DYMuonModExponentialFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, quadExpRooAbsPdf, true, "TTBar MuMuJJ Estimate with Quadratic Exponential fit to MuonEG data","MMuMuJJ_ttBar_data_driven_estimate_quadExponentialFit","m4","M_MuMuJJ [GeV]", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, false, DYMuonPdfIntegralFile, "DYMuonQuadExponentialFit integral from 2.0 to 3.0 TeV = ");

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, linearExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, false);

	//fitPdfReturnPdfAndResults(tempFitRslt, muMuMCDataSet, eMuRealDataSet, sumTwoExpRooAbsPdf, false, "","","","", massWR, 100, integralCheck, ttBarMcMuMutoEMuRatioIntercept*nEMuDataEvts, muonPdfs, muonIntegralChecks, false);

	
	RooWorkspace DYmuonFits("DYMuonFits");
	saveArgSetToRooWorkspace(DYmuonPdfs, DYmuonFits, false, "");
	//saveArgSetToRooWorkspace(DYmuonIntegralChecks, DYmuonFits, false, "");	///<import complains about this
	saveDataSetToRooWorkspace(DY_MuMuDataSet[0], DYmuonFits, true,"DYMuonFits.root");

#endif

}//commonRooFitBkgndEstimation()

