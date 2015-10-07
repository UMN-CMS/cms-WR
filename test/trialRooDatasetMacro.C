{

	RooRealVar massWR("fourObjectMass", "fourObjectMass", 0,3600);
	RooRealVar massNuR("subleadingLeptonThreeObjMass", "subleadingLeptonThreeObjMass", 0,2000);
	RooRealVar massDilepton("dileptonMass", "dileptonMass", 0,2500);
	RooRealVar genEvtWeights("evWeightSign", "evWeightSign", -2,2);
	
	RooArgSet vars(massWR,massNuR,massDilepton,genEvtWeights);

	///declare TChains used to construct RooDataSet objects for each bkgnd source
	TString treeName = "unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts";
	TString dirName = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	TChain * ttBarTree = new TChain(treeName,"");
	ttBarTree->Add(dirName+"analyzed_TTOnly_PowhegPythia_25ns_eejj_signal_region.root");
	/*
	TChain * dyJetsTree = new TChain(treeName,"");
	dyJetsTree->Add(dirName+"analyzed_DYJets_Madgraph_25ns_eejj_signal_region.root");
	TChain * singleTopWTree = new TChain(treeName,"");
	singleTopWTree->Add(dirName+"analyzed_SingleTopPlusW_25ns_eejj_signal_region.root");
	TChain * wJetsTree = new TChain(treeName,"");
	wJetsTree->Add(dirName+"analyzed_WJets_Madgraph_25ns_eejj_signal_region.root");
	//wJetsTree->Add(dirName+"analyzed_WJets_25ns_eejj_signal_region.root");	///< wJets aMCNLO
	TChain * wzTree = new TChain(treeName,"");
	wzTree->Add(dirName+"analyzed_WZ_25ns_eejj_signal_region.root");
	//wzTree->Add(dirName+"analyzed_WZPlusJets_25ns_eejj_signal_region.root");	///< wzPlusJets aMCNLO
	TChain * zzTree = new TChain(treeName,"");
	zzTree->Add(dirName+"analyzed_ZZ_25ns_eejj_signal_region.root");
	*/
	
	///declare RooDataSet objects, and add all of them into a single RooDataSet using append()
	//RooDataSet ttBar("ttBar","ttBar",ttBarTree, vars,"(evWeightSign != 0 ? evWeightSign = evWeightSign*(831.76*1000/19899500) : evWeightSign = evWeightSign*(831.76*1000/19899500) )","evWeightSign");
	RooDataSet ttBar("ttBar","ttBar",ttBarTree, vars,"","evWeightSign");
	
	/*
	RooDataSet dyJets("dyJets","dyJets",dyJetsTree, vars,"","evWeightSign");
	RooDataSet singleTopW("singleTopW","singleTopW",singleTopWTree, vars,"","evWeightSign");
	RooDataSet wJets("wJets","wJets",wJetsTree, vars,"","evWeightSign");
	RooDataSet wz("wz","wz",wzTree, vars,"","evWeightSign");
	RooDataSet zz("zz","zz",zzTree, vars,"","evWeightSign");

	ttBar.append(dyJets);
	ttBar.append(singleTopW);
	ttBar.append(wJets);
	ttBar.append(wz);
	ttBar.append(zz);
	*/

	//RooGuassian signalPDF();

	RooPlot *frame = massWR.frame();
	ttBar.plotOn(frame);


	//ttBar.Fit(signalPDF);
	//signalPDF.plotOn(frame);
	
	frame->Draw();

}

