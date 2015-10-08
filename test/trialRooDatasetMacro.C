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

	Float_t maxMassWR = 3600;	///< use this for the RooRealVar massWR, and the fit range near the end
	RooRealVar massWR("fourObjectMass", "fourObjectMass", 500,maxMassWR);
	RooRealVar genEvtWeights("evWeightSign", "evWeightSign", -2,2);

	
	RooArgSet vars(massWR,genEvtWeights);

	///declare TChains used to construct RooDataSet objects for each bkgnd source
	TString treeName = "unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts";
	TString dirName = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	TChain * dyJetsTree = new TChain(treeName,"");
	dyJetsTree->Add(dirName+"analyzed_DYJets_Madgraph_25ns_eejj_signal_region.root");
	TChain * ttBarTree = new TChain(treeName,"");
	ttBarTree->Add(dirName+"analyzed_TTOnly_PowhegPythia_25ns_eejj_signal_region.root");
	TChain * singleTopWTree = new TChain(treeName,"");
	singleTopWTree->Add(dirName+"analyzed_SingleTopPlusW_25ns_eejj_signal_region.root");
	
	/* only two events pass the WJets signal selection cuts!
	TChain * wJetsTree = new TChain(treeName,"");
	wJetsTree->Add(dirName+"analyzed_WJets_Madgraph_25ns_eejj_signal_region.root");
	//wJetsTree->Add(dirName+"analyzed_WJets_25ns_eejj_signal_region.root");	///< wJets aMCNLO
	*/

	TChain * wzTree = new TChain(treeName,"");
	wzTree->Add(dirName+"analyzed_WZ_25ns_eejj_signal_region.root");
	//wzTree->Add(dirName+"analyzed_WZPlusJets_25ns_eejj_signal_region.root");	///< wzPlusJets aMCNLO
	TChain * zzTree = new TChain(treeName,"");
	zzTree->Add(dirName+"analyzed_ZZ_25ns_eejj_signal_region.root");
	
	TChain * WRToEEJJTree = new TChain(treeName,"");
	WRToEEJJTree->Add(dirName+"analyzed_WRToENuToEEJJ_MWR_2600_MNu_1300_25ns_eejj_signal_region_no_gen_matching.root");	///< WR signal TChain

	///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
	Float_t intLumi = 1000.;	///< integrated lumi
	RooDataSet WR = applyNormalization(WRToEEJJTree, "WR",(0.015*intLumi/50000), vars, massWR, genEvtWeights);

	RooDataSet dyJets = applyNormalization(dyJetsTree, "dyJets",(6025.2*intLumi/9052671), vars, massWR, genEvtWeights);
	RooDataSet ttBar = applyNormalization(ttBarTree, "ttBar",(831.76*intLumi/19899500), vars, massWR, genEvtWeights);
	RooDataSet singleTopW = applyNormalization(singleTopWTree, "singleTopW",(35.6*intLumi/995600), vars, massWR, genEvtWeights);
	//RooDataSet wJets = applyNormalization(wJetsTree, "wJets",(61500*intLumi/72121586), vars, massWR, genEvtWeights);
	RooDataSet wz = applyNormalization(wzTree, "wz",(66.1*intLumi/991232), vars, massWR, genEvtWeights);
	RooDataSet zz = applyNormalization(zzTree, "zz",(15.4*intLumi/996168), vars, massWR, genEvtWeights);

	dyJets.append(ttBar);
	dyJets.append(singleTopW);
	//dyJets.append(wJets);
	dyJets.append(wz);
	dyJets.append(zz);
	dyJets.append(WR);
	

	/*
	WR.append(ttBar);
	WR.append(singleTopW);
	//WR.append(wJets);
	WR.append(wz);
	WR.append(zz);
	WR.append(dyJets);
	*/

	//RooGuassian somePDF();
	
	///alpha, power, meanMassWR, and sigmaMassWR are needed for crystal ball fit
	//RooRealVar alpha("alpha","",0,5);
	//RooRealVar power("power","",0,5);
	//RooRealVar meanMassWR("meanMassWR","",2500,2600);
	//RooRealVar sigmaMassWR("sigmaMassWR","",10,200);
	//RooCBShape signalPDF("wrFit","",massWR,meanMassWR,sigmaMassWR,alpha,power);
	
	RooRealVar bkgndAlpha("bkgndAlpha","",-2,2);
	RooRealVar bkgndPower("bkgndPower","",0,3);
	RooRealVar bkgndMeanMassWR("bkgndMeanMassWR","",600,700);
	RooRealVar bkgndSigmaMassWR("bkgndSigmaMassWR","",40,200);
	RooCBShape bkgndPDF("bkgndFit","",massWR,bkgndMeanMassWR,bkgndSigmaMassWR,bkgndAlpha,bkgndPower);


	RooPlot *frame = massWR.frame();
	frame->GetXaxis()->SetTitle("EEJJ Mass [GeV]");
	frame->SetTitle("EEJJ Mass for #SigmaBackground MCs and WR signal M_{WR} = 2600 GeV");
	dyJets.plotOn(frame);

	bkgndPDF.fitTo(dyJets);
	bkgndPDF.plotOn(frame);
	
	frame->Draw();
	TString plotDir = "tempPlots/RooDataSetMC25ns/";
	TString plotFileName = "eejj_mass_sum_of_bkgnds_and_WR_signal_M_2600_with_CB_fit.png";
	frame->SaveAs(plotDir+plotFileName,"recreate");

}///end macro()

