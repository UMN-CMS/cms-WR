//#define DEBUG

/*
 * use this fxn to loop over events in a TChain, pull the event weights from the TChain,
 * combine the weights with the cross section*luminosity (or some other normalization
 * factor), and store them in a RooDataSet.  This fxn returns a pointer to the RooDataSet.
 * the RooRealVar objects named fourObjMass and evWgt must be in the RooArgSet named dataSetVars
 */
RooDataSet applyNormalization(TChain * chain, TString dataSetName, Float_t normalization, RooArgSet & dataSetVars, RooRealVar & fourObjMass, RooRealVar & evWgt){
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
	
	///loop over all entries in the TChain
	for(Long64_t i=0; i<evts; i++){
		chain->GetEntry(i);

		///set the values of the RooRealVar objects using local vars tied to TChain branches
		fourObjMass.setVal(massWR_);
		evWgt.setVal(weight_*normalization);

		///add the values of the RooRealVar objects to the temporary RooDataSet
		tempDataSetPtr->add(dataSetVars);
	}///end loop over entries in TChain

	tempDataSetPtr->Print();
	chain->ResetBranchAddresses();
	
	///make the RooDataSet from the RooDataSet pointer, and identify the RooRealVar with weights
	RooDataSet dataSetWithWeights(dataSetName, dataSetName, tempDataSetPtr, dataSetVars,"", evWgt.GetName());
	return dataSetWithWeights;

}///end applyNormalization()

void macro(){

	Float_t maxMassWR = 13000;	///< use this for the RooRealVar massWR, and the fit range near the end
	RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,maxMassWR);
	RooRealVar genEvtWeights("evWeightSign", "evWeightSign", -2,2);

	
	RooArgSet vars(massWR,genEvtWeights);

	///declare TChains used to construct RooDataSet objects for each bkgnd source
	TString treeName = "unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts";
	TString dirName = "root://cmsxrootd.fnal.gov//store/user/skalafut/analyzed_25ns_eejj_signal_region/";

	
	TChain * WRToEEJJTree = new TChain(treeName,"");
	WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_2600_MNu_1300_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain
	//WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_6000_MNu_3000_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain

	///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
	Float_t intLumi = 1000.;	///< integrated lumi
	//RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",(0.015*intLumi/50000), vars, massWR, genEvtWeights);
	RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",(0.0000144*intLumi/16897), vars, massWR, genEvtWeights);


	/*
	WR.append(ttBar);
	WR.append(singleTopW);
	//WR.append(wJets);
	WR.append(wz);
	WR.append(zz);
	WR.append(dyJets);
	*/

	///alpha, power, meanMassWR, and sigmaMassWR are needed for crystal ball fit
	//RooRealVar alpha("alpha","",0,6);
	//RooRealVar power("power","",1,20);
	//RooRealVar meanMassWR("meanMassWR","",2400,2800);
	//RooRealVar sigmaMassWR("sigmaMassWR","",10,300);
	//RooCBShape CBLowToHighMass("CBLowToHighMass","",massWR,meanMassWR,sigmaMassWR,alpha,power);
	//RooCBShape signalPDF("wrFit","",massWR,meanMassWR,sigmaMassWR,alpha,power);
	
	//RooRealVar logNormMean("logNormMean","",2300,2850);
	//RooRealVar logNormShape("logNormShape","",1,50);
	//RooLognormal logNorm("logNorm","",massWR,logNormMean,logNormShape);
	
	//RooRealVar gaussMean("gaussMean","",1700,3000);
	//RooRealVar gaussWidth("gaussWidth","",30,300);
	//RooGaussian gaussHighMass("gaussHighMass","",massWR,gaussMean,gaussWidth);

	//define coefs for fxns
	//RooRealVar coefOne("coefOne","",0,1);
	//RooRealVar coefTwo("coefTwo","",0,1);
	//RooAddPdf signalPDF("wrFit","",RooArgList(logNorm,CBLowToHighMass),RooArgList(coefOne,coefTwo));



	/**/
	//RooRealVar mass("measWR", "", 600, 13000); 
	RooRealVar meanPeak("meanPeak", "", 1000, 600, maxMassWR);
	RooRealVar meanShiftSecondGauss("meanShiftSecondGauss","", 0.01, -0.3, 0.3);
	RooRealVar meanShiftThirdGauss("meanShiftThirdGauss","", 0.01, -0.5, 0.5);

	RooFormulaVar meanSecondGauss("meanSecondGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftSecondGauss));
	RooFormulaVar meanThirdGauss("meanThirdGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftThirdGauss));
   
	RooRealVar sigmaRelPeak("sigmaRelPeak", "", 50./2600, 10./2600, 300./2600);
	
	RooRealVar sigmaThirdGauss("sigmaThirdGauss", "sigmaThirdGauss", 600, 300, 1200);
	RooRealVar sigmaShiftSecondGaus("sigmaShiftSecondGaus", "", 1, 0.9, 3);

	RooFormulaVar sigmaPeak("sigmaPeak", "@0*@1", RooArgSet(meanPeak, sigmaRelPeak));

	RooFormulaVar sigmaSecondGauss("sigmaSecondGauss", "@0*@1", RooArgSet(sigmaPeak, sigmaShiftSecondGaus));
 
	RooRealVar meanRight("meanRight","",600,maxMassWR);		//right side of mass peak
	RooRealVar sigmaRight("sigmaRight","",10,300);
	RooGaussian gaussRight("gaussRight","",massWR,meanRight,sigmaRight);
	RooRealVar meanLeft("meanLeft","",600,maxMassWR);	///left side of mass peak
	RooRealVar sigmaLeft("sigmaLeft","",50,1000);
	RooGaussian gaussLeft("gaussLeft","",massWR,meanLeft,sigmaLeft);
	RooRealVar meanLowTail("meanLowTail","",600,maxMassWR);	///low mass tail
	RooRealVar sigmaLowTail("sigmaLowTail","",300,1000);
	RooGaussian gaussLowTail("gaussLowTail","",massWR,meanLowTail,sigmaLowTail);
	RooRealVar meanHighTail("meanHighTail","",600,maxMassWR);	///high mass tail
	RooRealVar sigmaHighTail("sigmaHighTail","",50,400);
	RooGaussian gaussHighTail("gaussHighTail","",massWR,meanHighTail,sigmaHighTail);
	
	RooGaussian gausPeak("gausPeak","", massWR, meanPeak, sigmaPeak);
	RooGaussian gausSecondGauss("gausSecondGauss", "", massWR, meanSecondGauss, sigmaSecondGauss);
	RooGaussian gausThirdGauss("gausThirdGauss", "", massWR, meanThirdGauss, sigmaThirdGauss);
	
	//define coefficients for the gaussian fxns
	RooRealVar coefRight("coefRight","",1,0,3.);
	RooRealVar coefLeft("coefLeft","",1,0,3.);
	RooRealVar coefLowTail("coefLowTail","",1,0,3.);
	RooRealVar coefHighTail("coefHighTail","",1,0,2.);
//	RooAddPdf signalPDF("wrFit","",RooArgList(gaussRight,gaussLeft,gaussLowTail,gaussHighTail),RooArgList(coefRight,coefLeft,coefLowTail,coefHighTail));
	//RooAddPdf signalPDF("wrFit","",RooArgList(gaussRight,gaussLeft,gaussLowTail),RooArgList(coefRight,coefLeft,coefLowTail));
	RooAddPdf signalPDF("wrFit","",RooArgList(gausPeak,gausSecondGauss,gausThirdGauss),RooArgList(coefRight,coefLeft, coefLowTail));
	//RooAddPdf signalPDF("wrFit","",RooArgList(logNorm,gaussLeft,gaussLowTail),RooArgList(coefRight,coefLeft,coefLowTail));


	/**/

	//Exponential
	//RooRealVar decay("decay","",-0.1,-1,0.);
	//RooExponential bkgndPDF("bkgndFit","",massWR,decay);

	RooPlot *frame = massWR.frame();
	frame->GetXaxis()->SetRangeUser(600, 9000);
	frame->GetXaxis()->SetTitle("EEJJ Mass [GeV]");
	frame->SetTitle("EEJJ Mass for WR MC M_{WR} = 2600 GeV  #intlumi = 1000/pb");
	//ttBar.plotOn(frame);
	WR.plotOn(frame);


	//bkgndPDF.fitTo(ttBar, RooFit::Range(650,1750,kTRUE));
	//bkgndPDF.plotOn(frame);
//	signalPDF.fitTo(WR, RooFit::Range(600,maxMassWR,kTRUE), RooFit::Verbose(kTRUE));
	signalPDF.fitTo(WR,  RooFit::Verbose(kTRUE));
	signalPDF.plotOn(frame, RooFit::LineColor(kRed));
	signalPDF.plotOn(frame, RooFit::Components(gausSecondGauss), RooFit::LineColor(kMagenta), RooFit::LineStyle(2));
	signalPDF.plotOn(frame, RooFit::Components(gausPeak), RooFit::LineColor(kGreen), RooFit::LineStyle(2));
	signalPDF.plotOn(frame, RooFit::Components(gausThirdGauss), RooFit::LineColor(kBlue), RooFit::LineStyle(2));
	signalPDF.plotOn(frame, RooFit::Components(gaussHighTail), RooFit::LineColor(kBlue), RooFit::LineStyle(1));

	TCanvas * c1 = new TCanvas("c1","c1",600,600);
	c1->cd();
	frame->Draw();
	TString plotDir = "tempPlots/RooDataSetMC25ns/";
	TString plotFileName = "eejj_mass_WR_signal_MWR_2600_MNu_1300_with_quadruple_gauss_and_four_floating_coefficients_fit.png";
	//TString plotFileName = "eejj_mass_WR_signal_MWR_2600_MNu_1300_with_quadruple_gauss_and_no_floating_coefficients_fit.png";
	
	c1->SaveAs(plotDir+plotFileName,"recreate");

}///end macro()

