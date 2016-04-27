//#define DEBUG
#define DOFIT

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

void bkgndFitRooMacro(){

	///declare TChains used to construct RooDataSet objects for each bkgnd source
	TString treeName = "recoAnalyzerOne/recoObjectsAllCuts";
	TString dirName = "/eos/uscms/store/user/skalafut/analyzed_25ns_emujj_signal_region/";

	//backgrounds
	TChain * dyJetsTree = new TChain(treeName,"");
	dyJetsTree->Add(dirName+"analyzed_DYJets_Madgraph_25ns_emujj_signal_region_reMiniAOD.root");
	TChain * ttBarTree = new TChain(treeName,"");
	ttBarTree->Add(dirName+"analyzed_TTOnly_PowhegPythia_25ns_emujj_signal_region_reMiniAOD.root");
	//TChain * wJetsTree = new TChain(treeName,"");
	//wJetsTree->Add(dirName+"analyzed_WJets_Madgraph_25ns_emujj_signal_region_reMiniAOD.root");
	TChain * wzTree = new TChain(treeName,"");
	wzTree->Add(dirName+"analyzed_WZ_25ns_emujj_signal_region_reMiniAOD.root");
	TChain * zzTree = new TChain(treeName,"");
	zzTree->Add(dirName+"analyzed_ZZ_25ns_emujj_signal_region_reMiniAOD.root");
	std::vector<Float_t> xSxnBkgnds = {6025.2,831.76,66.1,15.4};	///< elements are in this order: DYJets, TTBar, WZ, ZZ
	std::vector<Float_t> numEvtsBkgnds = {9042031,static_cast<float>(96834559),978512,996944};

	///make a RooDataset for each bkgnd source using applyNormalization(), then append all of the RooDataset objects to one RooDataset
	Float_t maxMassWR = 2100;	///< use this for the RooRealVar massWR, and the fit range near the end
	Float_t intLumi = 1570.674;
	RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,maxMassWR);
	RooRealVar evtWeight("evWeightSign", "evWeightSign", -2,2);
	RooArgSet vars(massWR,evtWeight);

	RooDataSet ttBar = applyNormalization(ttBarTree,"ttBar",intLumi*xSxnBkgnds[1]/numEvtsBkgnds[1], vars, massWR, evtWeight);
	RooDataSet dyJets = applyNormalization(dyJetsTree,"dyJets",intLumi*xSxnBkgnds[0]/numEvtsBkgnds[0], vars, massWR, evtWeight);
	RooDataSet wz = applyNormalization(wzTree,"wz",intLumi*xSxnBkgnds[2]/numEvtsBkgnds[2], vars, massWR, evtWeight);
	RooDataSet zz = applyNormalization(zzTree,"zz",intLumi*xSxnBkgnds[3]/numEvtsBkgnds[3], vars, massWR, evtWeight);

	ttBar.append(dyJets);
	ttBar.append(wz);
	ttBar.append(zz);

	///now define fit parameters and fit a function (exponential, quadratic, etc) to the massWR distribution
	///and plot the distribution with the fit
	RooPlot *frame = massWR.frame(300);
	frame->GetXaxis()->SetRangeUser(600,maxMassWR);
	frame->GetXaxis()->SetTitle("RECO EMuJJ Mass [GeV]");
	frame->SetTitle("RECO EMuJJ Mass #SigmaAll Backgrounds  #intlumi = 1570.7/pb");
	ttBar.plotOn(frame);

#ifdef DOFIT
	RooRealVar expPower("expPower","",-0.05,-0.2,0.);
	RooExponential expPdf("expPdf","", massWR, expPower);
	RooRealVar coef("coef","",100,0,600);
	RooAddPdf * scaledExpPdf = new RooAddPdf("fourObjMassFit","",RooArgList(expPdf),RooArgList(coef));
	
	//scaledExpPdf->fitTo(ttBar, RooFit::Verbose(kTRUE));
	scaledExpPdf->fitTo(ttBar);
	scaledExpPdf->plotOn(frame, RooFit::LineColor(kRed));
	//RooDataSet * scaledExpPdfDataSet = scaledExpPdf->generate(RooArgSet(massWR));
	/*
	Int_t numBins = 1000;
	TH1 * fitHisto = scaledExpPdf->createHistogram("fourObjectMass",numBins);
	TH1 * dataPtsHisto = ttBar.createHistogram("fourObjectMass",numBins);
	///loop over the bin contents of fitHisto and reset them to (fit - data)/data
	for(Int_t i=1; i<=fitHisto->GetNbinsX(); i++){
		Double_t dataBinContent = dataPtsHisto->GetBinContent(i);
		Double_t fitBinContent = fitHisto->GetBinContent(i);
		if(dataBinContent != 0){
			Double_t residualFrxn = (fitBinContent-dataBinContent)/dataBinContent;
			fitHisto->SetBinContent(i,residualFrxn);
			if(i==10||i==30||i==50) std::cout<<"residualFrxn =\t"<< residualFrxn << std::endl;
		}///end if(dataPtsHisto bin has nonzero entries)
		else{
			///set the bin content to 1 if there are no entries in dataPtsHisto, but there are entries in fitHisto
			if(fitBinContent == 0) fitHisto->SetBinContent(i,0);
			else fitHisto->SetBinContent(i,1);
		}

	}///end loop over fitHisto and dataPtsHisto bin contents
	fitHisto->SetTitle("RECO M_{EMuJJ} #Delta(Exp Fit, Data) / Data");
	fitHisto->GetXaxis()->SetTitle("RECO M_{EMuJJ} [GeV]");
	fitHisto->GetYaxis()->SetTitle("#Delta(Fit, Data) / Data");
	TCanvas * c2 = new TCanvas("c2","c2",800,800);
	c2->cd();
	fitHisto->Draw("");
	c2->SaveAs("tempPlots/rooFitBkgndsEMuSignalRegion/residualExpFitAndBkgndMCsEMuJJBkgndsSignalRegion.png","recreate");
	c2->SaveAs("tempPlots/rooFitBkgndsEMuSignalRegion/residualExpFitAndBkgndMCsEMuJJBkgndsSignalRegion.pdf","recreate");
	*/

#endif

	TCanvas * c1 = new TCanvas("c1","c1",800,800);
	c1->cd();
	frame->Draw();
	c1->Update();
	TString plotDir = "tempPlots/rooFitBkgndsEMuSignalRegion/";
	TString fitType = "noFit";
#ifdef DOFIT
	fitType = "two_DOF_scaledExpFit";
#endif
	TString plotFileName = "emujj_mass_all_bkgnds_signal_region_" + fitType;
	TString pngFile = ".png";
	TString pdfFile = ".pdf";

	c1->SaveAs(plotDir+plotFileName+pngFile,"recreate");
	c1->SaveAs(plotDir+plotFileName+pdfFile,"recreate");

	/*
	for(Int_t i=0; i<max; i++){
		///fit a function to each M_EEJJ distribution, save the image, then move on to a different input file
		RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,maxMassWR);
		RooRealVar evtWeight("evWeightSign", "evWeightSign", -2,2);

		RooArgSet vars(massWR,evtWeight);

		TChain * WRToEEJJTree = new TChain(treeName,"");
		WRToEEJJTree->Add(dirName+"wr" + genWrMass + "nu" + genNuMass + "Tree.root");	///< WR signal TChain

		///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
		RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",intLumi*xSxn[i]/50000, vars, massWR, evtWeight);	///<overall normalization won't affect shape of M_EEJJ distribution

		RooRealVar meanPeak("meanPeak", "", 1000, 600, maxMassWR);
		//mean of the two other gaussians must be less than the mean of the peak gaussian
		RooRealVar meanShiftSecondGauss("meanShiftSecondGauss","", -0.05, -0.3, 0.);
		RooRealVar meanShiftThirdGauss("meanShiftThirdGauss","", -0.06, -0.5, 0.);

		RooFormulaVar meanSecondGauss("meanSecondGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftSecondGauss));
		RooFormulaVar meanThirdGauss("meanThirdGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftThirdGauss));

		RooRealVar sigmaRelPeak("sigmaRelPeak", "", 50./2600, 10./2600, 300./2600);
		RooRealVar sigmaShiftSecondGaus("sigmaShiftSecondGaus", "", 1, 0.9, 3);
		RooRealVar sigmaShiftThirdGaus("sigmaShiftThirdGaus", "", 2, 1, 5);

		RooFormulaVar sigmaPeak("sigmaPeak", "@0*@1", RooArgSet(meanPeak, sigmaRelPeak));
		RooFormulaVar sigmaSecondGauss("sigmaSecondGauss", "@0*@1", RooArgSet(sigmaPeak, sigmaShiftSecondGaus));
		RooFormulaVar sigmaThirdGauss("sigmaThirdGauss", "300+@0*@1", RooArgSet(sigmaPeak, sigmaShiftThirdGaus));
		//RooRealVar sigmaThirdGauss("sigmaThirdGauss", "sigmaThirdGauss", 600, 300, 1200);

		RooGaussian gausPeak("gausPeak","", massWR, meanPeak, sigmaPeak);
		RooGaussian gausSecondGauss("gausSecondGauss", "", massWR, meanSecondGauss, sigmaSecondGauss);
		RooGaussian gausThirdGauss("gausThirdGauss", "", massWR, meanThirdGauss, sigmaThirdGauss);


		//define coefficients for the gaussian fxns
		RooRealVar coefRight("coefRight","",1,0,3.);
		RooRealVar coefLeft("coefLeft","",1,0,3.);
		RooRealVar coefLowTail("coefLowTail","",1,0,3.);
		RooRealVar coefHighTail("coefHighTail","",1,0,2.);
		//RooAddPdf signalPDF("wrFit","",RooArgList(gaussRight,gaussLeft,gaussLowTail,gaussHighTail),RooArgList(coefRight,coefLeft,coefLowTail,coefHighTail));
		RooAddPdf * signalPDF = new RooAddPdf("wrFit","",RooArgList(gausPeak,gausSecondGauss,gausThirdGauss),RooArgList(coefRight,coefLeft,coefLowTail));

		//RooAddPdf signalPDF("wrFit","",RooArgList(logNorm,gaussLeft,gaussLowTail),RooArgList(coefRight,coefLeft,coefLowTail));

		RooPlot *frame = massWR.frame(300);
		frame->GetXaxis()->SetRangeUser(600,maxMassWR);
		frame->GetXaxis()->SetTitle("EEJJ Mass [GeV]");
		frame->SetTitle("EEJJ Mass for WR MC M_{WR} = " + genWrMass + " GeV  #intlumi = 1000/pb");
		//frame->SetTitle("EEJJ Mass for WR MC M_{WR} = " + genWrMass + " GeV   arbitrary normalization");
		WR.plotOn(frame);


		//signalPDF.fitTo(WR, RooFit::Verbose(kTRUE));
		signalPDF->fitTo(WR);
		signalPDF->plotOn(frame, RooFit::LineColor(kRed));
		signalPDF->plotOn(frame, RooFit::Components(gausSecondGauss), RooFit::LineColor(kMagenta), RooFit::LineStyle(2));
		signalPDF->plotOn(frame, RooFit::Components(gausPeak), RooFit::LineColor(kGreen), RooFit::LineStyle(2));
		signalPDF->plotOn(frame, RooFit::Components(gausThirdGauss), RooFit::LineColor(kBlue), RooFit::LineStyle(2));

		TCanvas * c1 = new TCanvas("c1","c1",600,600);
		c1->cd();
		frame->Draw();
		c1->Update();
		TString plotDir = "tempPlots/RooDataSetMC25ns/";
		TString plotFileName = "eejj_mass_WR_signal_MWR_" + genWrMass + "_MNu_halfMWR_with_triple_gauss_and_three_floating_coefficients_fit.png";
		//TString plotFileName = "eejj_mass_WR_signal_MWR_" + genWrMass + "_MNu_halfMWR_with_triple_gauss_and_no_floating_coefficients_fit.png";

		c1->SaveAs(plotDir+plotFileName,"recreate");
		signalPDF->Delete();

	}///end loop over all input files

	*/
	

	/*
	Float_t maxRange = 1.6*(genWrMass.Atof());	///<max value of M_LLJJ shown in RooPlot (RooPlot pointer named frame)
	//WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_2600_MNu_1300_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain
	//WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_6000_MNu_3000_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain

	//TString fileName = "wr"+genWrMass+"nu"+genNuMass+"Tree.root";
	//std::cout<<"opening file named:\t"<< fileName <<std::endl;

	WRToEEJJTree->Add(dirName+"wr" + genWrMass + "nu" + genNuMass + "Tree.root");	///< WR signal TChain

	///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
	Float_t intLumi = 1000.;	///< integrated lumi
	RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",0.5, vars, massWR, evtWeight);	///<overall normalization won't affect shape of M_EEJJ distribution



	//RooDataSet dyJets = applyNormalization(dyJetsTree, "dyJets",(6025.2*intLumi/9052671), vars, massWR, evtWeight);
	//RooDataSet ttBar = applyNormalization(ttBarTree, "ttBar",(831.76*intLumi/19899500), vars, massWR, evtWeight);
	//RooDataSet singleTopW = applyNormalization(singleTopWTree, "singleTopW",(35.6*intLumi/995600), vars, massWR, evtWeight);
	////RooDataSet wJets = applyNormalization(wJetsTree, "wJets",(61500*intLumi/72121586), vars, massWR, evtWeight);
	//RooDataSet wz = applyNormalization(wzTree, "wz",(66.1*intLumi/991232), vars, massWR, evtWeight);
	//RooDataSet zz = applyNormalization(zzTree, "zz",(15.4*intLumi/996168), vars, massWR, evtWeight);

	//ttBar.append(dyJets);
	//ttBar.append(singleTopW);
	//////ttBar.append(wJets);
	//ttBar.append(wz);
	//ttBar.append(zz);
	//ttBar.append(WR);


	//WR.append(ttBar);
	//WR.append(singleTopW);
	////WR.append(wJets);
	//WR.append(wz);
	//WR.append(zz);
	//WR.append(dyJets);

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
	RooRealVar meanShiftSecondGauss("meanShiftSecondGauss","", -0.05, -0.3, 0.);
	RooRealVar meanShiftThirdGauss("meanShiftThirdGauss","", -0.2, -3.5, 0.);

	RooFormulaVar meanSecondGauss("meanSecondGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftSecondGauss));
	RooFormulaVar meanThirdGauss("meanThirdGauss", "", "@0+@1*@0", RooArgSet(meanPeak, meanShiftThirdGauss));

	RooRealVar sigmaRelPeak("sigmaRelPeak", "", 50./2600, 10./2600, 300./2600);
	RooRealVar sigmaShiftSecondGaus("sigmaShiftSecondGaus", "", 1, 0.9, 6);
	RooRealVar sigmaShiftThirdGaus("sigmaShiftThirdGaus", "", 2, 1, 30);

	RooFormulaVar sigmaPeak("sigmaPeak", "@0*@1", RooArgSet(meanPeak, sigmaRelPeak));
	RooFormulaVar sigmaSecondGauss("sigmaSecondGauss", "@0*@1", RooArgSet(sigmaPeak, sigmaShiftSecondGaus));
	RooFormulaVar sigmaThirdGauss("sigmaThirdGauss", "220+@0*@1", RooArgSet(sigmaPeak, sigmaShiftThirdGaus));
	//RooRealVar sigmaThirdGauss("sigmaThirdGauss", "sigmaThirdGauss", 600, 300, 1200);
	
	RooGaussian gausPeak("gausPeak","", massWR, meanPeak, sigmaPeak);
	RooGaussian gausSecondGauss("gausSecondGauss", "", massWR, meanSecondGauss, sigmaSecondGauss);
	RooGaussian gausThirdGauss("gausThirdGauss", "", massWR, meanThirdGauss, sigmaThirdGauss);


	//define coefficients for the gaussian fxns
	RooRealVar coefRight("coefRight","",1,0,10.);
	RooRealVar coefLeft("coefLeft","",1,0,10.);
	RooRealVar coefLowTail("coefLowTail","",1,0,10.);
	RooRealVar coefHighTail("coefHighTail","",1,0,2.);
	//RooAddPdf signalPDF("wrFit","",RooArgList(gaussRight,gaussLeft,gaussLowTail,gaussHighTail),RooArgList(coefRight,coefLeft,coefLowTail,coefHighTail));
	RooAddPdf signalPDF("wrFit","",RooArgList(gausPeak,gausSecondGauss,gausThirdGauss),RooArgList(coefRight,coefLeft,coefLowTail));
	
	//RooAddPdf signalPDF("wrFit","",RooArgList(logNorm,gaussLeft,gaussLowTail),RooArgList(coefRight,coefLeft,coefLowTail));



	//Exponential
	//RooRealVar decay("decay","",-0.1,-1,0.);
	//RooExponential bkgndPDF("bkgndFit","",massWR,decay);

	RooPlot *frame = massWR.frame();
	frame->GetXaxis()->SetRangeUser(600,maxRange);
	frame->GetXaxis()->SetTitle("EEJJ Mass [GeV]");
	//frame->SetTitle("EEJJ Mass for WR MC M_{WR} = " + genWrMass + " GeV  #intlumi = 1000/pb");
	frame->SetTitle("EEJJ Mass for WR MC M_{WR} = " + genWrMass + " GeV   arbitrary normalization");
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

	*/

}///end BkgndFitRooMacro()

