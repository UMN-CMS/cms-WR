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
	TString dirName = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";

	//backgrounds
	/*
	TChain * dyJetsTree = new TChain(treeName,"");
	dyJetsTree->Add(dirName+"analyzed_DYJets_Madgraph_M_50_25ns_eejj_signal_region.root");
	TChain * ttBarTree = new TChain(treeName,"");
	ttBarTree->Add(dirName+"analyzed_TTOnly_PowhegPythia_25ns_eejj_signal_region_reMiniAOD.root");
	TChain * singleTopWTree = new TChain(treeName,"");
	singleTopWTree->Add(dirName+"analyzed_SingleTopPlusW_25ns_eejj_signal_region_reMiniAOD.root");
	*/

	/* only two events pass the WJets signal selection cuts!
	TChain * wJetsTree = new TChain(treeName,"");
	wJetsTree->Add(dirName+"analyzed_WJets_Madgraph_25ns_eejj_signal_region.root");
	//wJetsTree->Add(dirName+"analyzed_WJets_25ns_eejj_signal_region.root");	///< wJets aMCNLO
	*/

	/*
	TChain * wzTree = new TChain(treeName,"");
	wzTree->Add(dirName+"analyzed_WZ_25ns_eejj_signal_region.root");
	//wzTree->Add(dirName+"analyzed_WZPlusJets_25ns_eejj_signal_region.root");	///< wzPlusJets aMCNLO
	TChain * zzTree = new TChain(treeName,"");
	zzTree->Add(dirName+"analyzed_ZZ_25ns_eejj_signal_region.root");
	*/
	
	TChain * WRToEEJJTree = new TChain(treeName,"");
	TString genWrMass = "2600";	///<gen WR mass shown in plot titles and names of saved images
	Float_t maxRange = 3200;	///<max value of M_LLJJ shown in RooPlot (RooPlot pointer named frame)
	WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_2600_MNu_1300_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain
	//WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_6000_MNu_3000_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain

	///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
	Float_t intLumi = 1000.;	///< integrated lumi
	RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",(0.015*intLumi/50000), vars, massWR, genEvtWeights);
	//RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",(0.0000144*intLumi/16897), vars, massWR, genEvtWeights);


	/*
	RooDataSet dyJets = applyNormalization(dyJetsTree, "dyJets",(6025.2*intLumi/9052671), vars, massWR, genEvtWeights);
	RooDataSet ttBar = applyNormalization(ttBarTree, "ttBar",(831.76*intLumi/19899500), vars, massWR, genEvtWeights);
	RooDataSet singleTopW = applyNormalization(singleTopWTree, "singleTopW",(35.6*intLumi/995600), vars, massWR, genEvtWeights);
	//RooDataSet wJets = applyNormalization(wJetsTree, "wJets",(61500*intLumi/72121586), vars, massWR, genEvtWeights);
	RooDataSet wz = applyNormalization(wzTree, "wz",(66.1*intLumi/991232), vars, massWR, genEvtWeights);
	RooDataSet zz = applyNormalization(zzTree, "zz",(15.4*intLumi/996168), vars, massWR, genEvtWeights);

	//ttBar.append(dyJets);
	ttBar.append(singleTopW);
	////ttBar.append(wJets);
	ttBar.append(wz);
	ttBar.append(zz);
	//ttBar.append(WR);
	*/


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

	RooRealVar meanPeak("meanPeak", "", 1000, 600, maxMassWR);
	//mean of the two other gaussians must be less than the mean of the peak gaussian
	RooRealVar meanShiftSecondGauss("meanShiftSecondGauss","", -0.01, -0.3, 0.);
	RooRealVar meanShiftThirdGauss("meanShiftThirdGauss","", -0.01, -0.5, 0.);

	RooFormulaVar meanSecondGauss("meanSecondGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftSecondGauss));
	RooFormulaVar meanThirdGauss("meanThirdGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftThirdGauss));

	RooRealVar sigmaRelPeak("sigmaRelPeak", "", 50./2600, 10./2600, 300./2600);
	RooRealVar sigmaShiftSecondGaus("sigmaShiftSecondGaus", "", 1, 0.9, 3);
	RooRealVar sigmaShiftThirdGaus("sigmaShiftThirdGaus", "", 1, 0.1, 3);

	RooFormulaVar sigmaPeak("sigmaPeak", "@0*@1", RooArgSet(meanPeak, sigmaRelPeak));
	RooFormulaVar sigmaSecondGauss("sigmaSecondGauss", "@0*@1", RooArgSet(sigmaPeak, sigmaShiftSecondGaus));
	RooFormulaVar sigmaThirdGauss("sigmaThirdGauss", "300+@0*@1", RooArgSet(sigmaPeak, sigmaShiftThirdGaus));
	//RooRealVar sigmaThirdGauss("sigmaThirdGauss", "sigmaThirdGauss", 600, 300, 1200);
	

	/*
	 * these are no longer useful.  they were tied to a specific mass point, and were not flexible with respect to expected M_WR
	RooRealVar meanRight("meanRight","",2450,maxMassWR);		//right side of mass peak
	RooRealVar sigmaRight("sigmaRight","",30,1000);
	RooGaussian gaussRight("gaussRight","",massWR,meanRight,sigmaRight);
	RooRealVar meanLeft("meanLeft","",2300,maxMassWR);	///left side of mass peak
	RooRealVar sigmaLeft("sigmaLeft","",50,1000);
	RooGaussian gaussLeft("gaussLeft","",massWR,meanLeft,sigmaLeft);
	RooRealVar meanLowTail("meanLowTail","",1100,maxMassWR);	///low mass tail
	RooRealVar sigmaLowTail("sigmaLowTail","",150,1000);
	RooGaussian gaussLowTail("gaussLowTail","",massWR,meanLowTail,sigmaLowTail);
	RooRealVar meanHighTail("meanHighTail","",2500,maxMassWR);	///high mass tail
	RooRealVar sigmaHighTail("sigmaHighTail","",50,1000);
	RooGaussian gaussHighTail("gaussHighTail","",massWR,meanHighTail,sigmaHighTail);
	*/

	RooGaussian gausPeak("gausPeak","", massWR, meanPeak, sigmaPeak);
	RooGaussian gausSecondGauss("gausSecondGauss", "", massWR, meanSecondGauss, sigmaSecondGauss);
	RooGaussian gausThirdGauss("gausThirdGauss", "", massWR, meanThirdGauss, sigmaThirdGauss);


	//define coefficients for the gaussian fxns
	RooRealVar coefRight("coefRight","",1,0,3.);
	RooRealVar coefLeft("coefLeft","",1,0,3.);
	RooRealVar coefLowTail("coefLowTail","",1,0,3.);
	RooRealVar coefHighTail("coefHighTail","",1,0,2.);
	//RooAddPdf signalPDF("wrFit","",RooArgList(gaussRight,gaussLeft,gaussLowTail,gaussHighTail),RooArgList(coefRight,coefLeft,coefLowTail,coefHighTail));
	RooAddPdf signalPDF("wrFit","",RooArgList(gausPeak,gausSecondGauss,gausThirdGauss),RooArgList(coefRight,coefLeft,coefLowTail));
	
	//RooAddPdf signalPDF("wrFit","",RooArgList(logNorm,gaussLeft,gaussLowTail),RooArgList(coefRight,coefLeft,coefLowTail));


	/**/

	//Exponential
	//RooRealVar decay("decay","",-0.1,-1,0.);
	//RooExponential bkgndPDF("bkgndFit","",massWR,decay);

	RooPlot *frame = massWR.frame();
	frame->GetXaxis()->SetRangeUser(600,maxRange);
	frame->GetXaxis()->SetTitle("EEJJ Mass [GeV]");
	frame->SetTitle("EEJJ Mass for WR MC M_{WR} = " + genWrMass + " GeV  #intlumi = 1000/pb");
	//ttBar.plotOn(frame);
	WR.plotOn(frame);


	//bkgndPDF.fitTo(ttBar);
	//bkgndPDF.plotOn(frame);
	signalPDF.fitTo(WR, RooFit::Verbose(kTRUE));
	signalPDF.plotOn(frame, RooFit::LineColor(kRed));
	signalPDF.plotOn(frame, RooFit::Components(gausSecondGauss), RooFit::LineColor(kMagenta), RooFit::LineStyle(2));
	signalPDF.plotOn(frame, RooFit::Components(gausPeak), RooFit::LineColor(kGreen), RooFit::LineStyle(2));
	signalPDF.plotOn(frame, RooFit::Components(gausThirdGauss), RooFit::LineColor(kBlue), RooFit::LineStyle(2));

	TCanvas * c1 = new TCanvas("c1","c1",600,600);
	c1->cd();
	frame->Draw();
	TString plotDir = "tempPlots/RooDataSetMC25ns/";
	TString plotFileName = "eejj_mass_WR_signal_MWR_" + genWrMass + "_MNu_halfMWR_with_triple_gauss_and_three_floating_coefficients_fit.png";
	//TString plotFileName = "eejj_mass_WR_signal_MWR_" + genWrMass + "_MNu_halfMWR_with_triple_gauss_and_no_floating_coefficients_fit.png";
	
	c1->SaveAs(plotDir+plotFileName,"recreate");

}///end macro()

