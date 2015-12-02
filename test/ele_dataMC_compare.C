#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TColor.h"
#include "TLatex.h"
#include "TPaveText.h"
#include <vector>
#include <iostream>
#include <string>

//#define DEBUG
//#define PrintIntegral
//#define CHKTTBAR

TString DetermineYaxisName(TH1F * ptrDataHist, TString xLabel);
void Fill_Histo(std::vector<TH1F*> h1, TTree* tree, std::vector<float> PUW, bool pileup_reweight, bool is_data);
std::vector<TH1F*> MakeNHistos(TString hname, int n, int bins, float x_min, float x_max);
std::vector<float> PileUpWeights(TTree* tree,TTree* tree_data);

void ele_dataMC_compare(){

#ifdef DEBUG
	std::cout<<"in ele_dataMC_compare()"<<std::endl;
#endif
	
	TString directory = "/eos/uscms/store/user/skalafut/analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj/";
	TFile * hfile0 = new TFile(directory+"analyzed_DYJets_Madgraph_25ns_skim_low_mass_region_eejj_reMiniAOD.root");//dyjets
	//TFile * hfile1 = new TFile(directory+"analyzed_TTJets_25ns_skim_low_mass_region_eejj.root");//ttbar
	TFile * hfile1 = new TFile(directory+"analyzed_TTOnly_PowhegPythia_25ns_skim_low_mass_region_eejj_reMiniAOD.root");//ttbar
	TFile * hfile3 = new TFile(directory+"analyzed_WZ_25ns_skim_low_mass_region_eejj.root");//wz
	TFile * hfile4 = new TFile(directory+"analyzed_ZZ_25ns_skim_low_mass_region_eejj.root");//zz
	TFile * hfile5 = new TFile(directory+"analyzed_WJets_Madgraph_25ns_skim_low_mass_region_eejj_reMiniAOD.root");//wjets

	TFile * hfile_data = new TFile(directory+"analyzed_DoubleEG_25ns_skim_low_mass_region_eejj_Run2015D_all_silver.root");//data

	int numProcesses = 6;

#ifdef DEBUG
	std::cout<<"declared pointers to input files"<<std::endl;
#endif

	TString treeName = "recoAnalyzerTwo/recoObjectsWithPtEtaCuts";
	TTree *tree0 = (TTree*)hfile0->Get(treeName);  
	TTree *tree1 = (TTree*)hfile1->Get(treeName);  
	//TTree *tree2 = (TTree*)hfile2->Get(treeName);  
	TTree *tree3 = (TTree*)hfile3->Get(treeName);  
	TTree *tree4 = (TTree*)hfile4->Get(treeName);  
	TTree *tree5 = (TTree*)hfile5->Get(treeName);  

	TTree *tree_data = (TTree*)hfile_data->Get(treeName);  

#ifdef DEBUG
	std::cout<<"declared pointers to trees in input files"<<std::endl;
#endif

	///set histo names, number of bins, and axis limits here
	vector<TH1F*> h_Mlljj = MakeNHistos("h_Mlljj",numProcesses,20,0,650);
	vector<TH1F*> h_Mll = MakeNHistos("h_Mll",numProcesses,22,0,220);
	vector<TH1F*> h_l1pt = MakeNHistos("h_l1pt",numProcesses,22,0,220);
	vector<TH1F*> h_l2pt = MakeNHistos("h_l2pt",numProcesses,14,0,140);
	vector<TH1F*> h_j1pt = MakeNHistos("h_j1pt",numProcesses,22,0,220);
	vector<TH1F*> h_j2pt = MakeNHistos("h_j2pt",numProcesses,14,0,140);
	vector<TH1F*> h_l1eta = MakeNHistos("h_l1eta", numProcesses,20,-3.,3.);
	vector<TH1F*> h_l2eta = MakeNHistos("h_l2eta", numProcesses,20,-3.,3.);
	vector<TH1F*> h_j1eta = MakeNHistos("h_j1eta", numProcesses,20,-3.,3.);
	vector<TH1F*> h_j2eta = MakeNHistos("h_j2eta", numProcesses,20,-3.,3.);
	vector<TH1F*> h_l1phi = MakeNHistos("h_l1phi", numProcesses,20,-3.15,3.15);
	vector<TH1F*> h_l2phi = MakeNHistos("h_l2phi", numProcesses,20,-3.15,3.15);
	vector<TH1F*> h_j1phi = MakeNHistos("h_j1phi", numProcesses,20,-3.15,3.15);
	vector<TH1F*> h_j2phi = MakeNHistos("h_j2phi", numProcesses,20,-3.15,3.15);
	vector<TH1F*> h_nleptons = MakeNHistos("h_nleptons", numProcesses,4,0,4);
	vector<TH1F*> h_njets = MakeNHistos("h_njets", numProcesses,6,0,6);
	vector<TH1F*> h_nvertices = MakeNHistos("h_nvertices", numProcesses,33,0,33);
	vector<TH1F*> h_dR_l1l2 = MakeNHistos("h_dR_l1l2", numProcesses,20,0,5);
	vector<TH1F*> h_dR_j1j2 = MakeNHistos("h_dR_j1j2", numProcesses,20,0,5);
	vector<TH1F*> h_dR_l1j1 = MakeNHistos("h_dR_l1j1", numProcesses,20,0,5);
	vector<TH1F*> h_dR_l1j2 = MakeNHistos("h_dR_l1j2", numProcesses,20,0,5);
	vector<TH1F*> h_dR_l2j1 = MakeNHistos("h_dR_l2j1", numProcesses,20,0,5);
	vector<TH1F*> h_dR_l2j2 = MakeNHistos("h_dR_l2j2", numProcesses,20,0,5);
	vector<TH1F*> h_Mjj = MakeNHistos("h_Mjj",numProcesses,40,0,400);
	vector<TH1F*> h_MET = MakeNHistos("h_MET",numProcesses,20,0,400);
	
#ifdef DEBUG
	std::cout<<"made vectors of TH1F pointers to histos"<<std::endl;
#endif

	int nhistos = 25;	//max is 25
	std::vector<THStack*> ths; // Stacks;
	for(int istacks=0; istacks<nhistos; istacks++)
	{
		TString full_name = Form("th%d",istacks);      
		THStack *th = new THStack(full_name.Data(),"");
		ths.push_back(th);
	}
	
#ifdef DEBUG
	std::cout<<"made vector of THStack pointers"<<std::endl;
#endif

	///declare one vector of TH1F pointers for each MC process and real data sample
	///the number of elements in each vector should be equal to the number of unique vector<TH1F*> containers which were declared above 
	std::vector<TH1F*> histos0(nhistos); // DY
	std::vector<TH1F*> histos1(nhistos); // TTbar
	//std::vector<TH1F*> histos2(nhistos); // WR signal
	std::vector<TH1F*> histos3(nhistos); // WZ
	std::vector<TH1F*> histos4(nhistos); // ZZ
	std::vector<TH1F*> histos5(nhistos); // WJets
	std::vector<TH1F*> histos_data(nhistos); // data

	///load the vectors tied to all MC and real data samples into one vector<vector<TH1F*> > object
	vector<vector<TH1F*> > histos(numProcesses);
	histos[0] = histos0;
	histos[1] = histos1;
	histos[2] = histos3;
	histos[3] = histos4;
	histos[4] = histos5;
	histos[numProcesses-1] = histos_data;

	///link the histos made by calls to MakeNHistos() above to specific elements of vector<vector<TH1F*> > container 
	for(int j=0; j<numProcesses; j++){
		histos[j][0] = h_Mlljj[j];
		histos[j][1] = h_Mll[j];
		histos[j][2] = h_l1pt[j];
		histos[j][3] = h_l2pt[j];
		histos[j][4] = h_j1pt[j];
		histos[j][5] = h_j2pt[j];
		histos[j][6] = h_l1eta[j];
		histos[j][7] = h_l2eta[j];
		histos[j][8] = h_j1eta[j];
		histos[j][9] = h_j2eta[j];
		histos[j][10] = h_l1phi[j];
		histos[j][11] = h_l2phi[j];
		histos[j][12] = h_j1phi[j];
		histos[j][13] = h_j2phi[j];
		histos[j][14] = h_nleptons[j];
		histos[j][15] = h_njets[j];
		histos[j][16] = h_nvertices[j];
		histos[j][17] = h_dR_l1l2[j];
		histos[j][18] = h_dR_j1j2[j];
		histos[j][19] = h_dR_l1j1[j];
		histos[j][20] = h_dR_l1j2[j];
		histos[j][21] = h_dR_l2j1[j];
		histos[j][22] = h_dR_l2j2[j];
		histos[j][23] = h_Mjj[j];
		histos[j][24] = h_MET[j];
	
	}
	
#ifdef DEBUG
	std::cout<<"linked MC processes and real data to histos of specific quantities"<<std::endl;
#endif

	///calculate PU weights and store them in persistent vectors of floats
	std::vector<float> PUW0 = PileUpWeights(tree0,tree_data);
	std::vector<float> PUW1 = PileUpWeights(tree1,tree_data);
	std::vector<float> PUW3 = PileUpWeights(tree3,tree_data);
	std::vector<float> PUW4 = PileUpWeights(tree4,tree_data);
	std::vector<float> PUW5 = PileUpWeights(tree5,tree_data);
	std::vector<float> PUW_data = PileUpWeights(tree_data,tree_data);

	///fill the histos with content from TTrees
	Fill_Histo(histos[0],tree0,PUW0,true,false); // DY
	Fill_Histo(histos[1],tree1,PUW1,true,false); // TTbar
	Fill_Histo(histos[2],tree3,PUW3,true,false); // WZ
	Fill_Histo(histos[3],tree4,PUW4,true,false); // ZZ
	Fill_Histo(histos[4],tree5,PUW5,true,false); // WJets

	Fill_Histo(histos[numProcesses-1],tree_data,PUW_data,false,true);	///real data

	//Float_t intLumi = 64.11;	///50ns from Run2015B and 2015C
	//Float_t intLumi = 1570.674;	///25ns from Run2015C and 2015D golden
	Float_t intLumi = 351.885;	///25ns from Run2015D silver
	// Scale = xsection*luminosity/events
	for(std::vector<TH1F*>::size_type i = 0; i != nhistos; i++){
#ifdef DEBUG
		std::cout<<"setting line and fill colors, and adding histos to THStack objects"<<std::endl;
#endif
		Double_t bkgndIntegral = 0;	///< integral of all bkgnd MC histos
		
		//histos[0][i]->Scale(6025.2*(intLumi)/19925500);	///aMC@NLO DYJets 50ns
		//histos[0][i]->Scale(6025.2*(intLumi)/9051899);	///madgraph DYJets 50ns
		histos[0][i]->Scale(6025.2*(intLumi)/9042031);	///madgraph DYJets 25ns
		histos[0][i]->SetFillColor(5);
		bkgndIntegral += histos[0][i]->Integral();
		
		//histos[1][i]->Scale(815.96*(intLumi)/4994250);	///aMCatNLO ttBar to all 50ns
		//histos[1][i]->Scale(831.76*(intLumi)/19665194);	///powheg-pythia ttBar to all 50ns
		histos[1][i]->Scale(831.76*(intLumi)/96834559);	///powheg-pythia ttBar to all 25ns
		//histos[1][i]->Scale(57.35*(intLumi)/24512786);	///ttBar to dilepton 25ns
		histos[1][i]->SetFillColor(3);
		bkgndIntegral += histos[1][i]->Integral();
#ifdef CHKTTBAR
		std::cout<<"ee chnl ttBar integral =\t"<< histos[1][i]->Integral() <<std::endl;
#endif
	
		//histos[2][i]->Scale(66.1*(intLumi)/996920);	///WZ to all 50ns
		//histos[2][i]->Scale(5.52*(intLumi)/31054519);	///WZ to 2L2Q 25ns
		histos[2][i]->Scale(66.1*(intLumi)/978512);	///WZ to all 25ns
		
		histos[2][i]->SetFillColor(4);
		bkgndIntegral += histos[2][i]->Integral();
	
	
		histos[3][i]->Scale(15.4*(intLumi)/996944);		///ZZ to all 25ns
		//histos[3][i]->Scale(15.4*(intLumi)/998848);		///ZZ to all 50ns
		//histos[3][i]->Scale(3.38*(intLumi)/18898680);		///ZZ to 2L2Q 25ns
		histos[3][i]->SetFillColor(7);
		bkgndIntegral += histos[3][i]->Integral();
		
		//histos[4][i]->Scale(6.15e4*(intLumi)/24089991);		///WJetsToLNu 50ns
		histos[4][i]->Scale(6.15e4*(intLumi)/72207128);		///WJetsToLNu 25ns
		histos[4][i]->SetFillColor(6);
		bkgndIntegral += histos[4][i]->Integral();
	
		histos[numProcesses-1][i]->SetMarkerStyle(20);
#ifdef PrintIntegral
		std::cout<<"integral of bkgnd histos =\t"<< bkgndIntegral <<std::endl;
		std::cout<<"integral of real data histo =\t"<< histos[numProcesses-1][i]->Integral() <<std::endl;
#endif

		ths[i]->Add(histos[3][i]);
		ths[i]->Add(histos[2][i]);
		ths[i]->Add(histos[4][i]);
		ths[i]->Add(histos[1][i]);
		ths[i]->Add(histos[0][i]);
	
		/*
		///rescale the max y value to make room for the legend
		Double_t oldMax = ( (ths[i]->GetMaximum() > histos[5][i]->GetMaximum() ) ? ths[i]->GetMaximum() : histos[5][i]->GetMaximum() );
		ths[i]->SetMaximum(20*oldMax);
		histos[5][i]->SetMaximum(20*oldMax);
		//ths[i]->SetMaximum(1.5*oldMax);
		//histos[5][i]->SetMaximum(1.5*oldMax);
		ths[i]->SetMinimum(0.05);
		histos[5][i]->SetMinimum(0.05);
		*/
	}


	///make a legend with appropriate labels for MC processes
	TLegend *leg = new TLegend( 0.62, 0.64, 0.89, 0.87 ) ;
	leg->SetNColumns(2);
	leg->AddEntry( histos[numProcesses-1][0], "Data" ,"ep") ;
	leg->AddEntry( histos[0][0], "DY" ) ; 
	leg->AddEntry( histos[1][0], "ttbar" ) ;
	leg->AddEntry( histos[4][0], "WJets" ) ;  
	leg->AddEntry( histos[2][0], "WZ" ) ; 
	leg->AddEntry( histos[3][0], "ZZ" ) ;
	leg->SetFillColor( kWhite ) ;

	TString xtitles[] = {"M_{EEJJ} [GeV]","M_{EE} [GeV]","leading electron p_{T} [GeV]","subleading electron p_{T} [GeV]","leading jet p_{T} [GeV]","subleading jet p_{T} [GeV]","leading electron #eta","subleading electron #eta","leading jet #eta","subleading jet #eta","leading electron #phi","subleading electron #phi","leading jet #phi","subleading jet #phi","number of electrons","number of jets","number of vertices","#DeltaR lead ele sublead ele","#DeltaR lead jet sublead jet","#DeltaR lead ele lead jet","#DeltaR lead ele sublead jet","#DeltaR sublead ele lead jet","#DeltaR sublead ele sublead jet","M_{JJ} [GeV]","Missing E_{T} [GeV]"};
	
	TString titles[] = {"CMS Preliminary M_{EEJJ}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary DiElectron Mass  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Lead Electron p_{T}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Sublead Electron p_{T}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Lead Jet p_{T}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Sublead Jet p_{T}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Lead Electron #eta  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Sublead Electron #eta  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Lead Jet #eta  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Subleading jet #eta  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary leading electron #phi  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Subleading electron #phi  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary leading jet #phi  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Subleading jet #phi  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary number of electrons  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary number of jets  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary number of vertices  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR lead ele Sublead ele  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR lead jet Sublead jet  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR lead ele lead jet  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR lead ele Sublead jet  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR Sublead ele lead jet  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary #DeltaR Sublead ele Sublead jet  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Dijet Mass  #surds = 13 TeV 25ns  #intlumi = 351.9/pb","CMS Preliminary Missing E_{T}  #surds = 13 TeV 25ns  #intlumi = 351.9/pb"};

	TString fnames[] = {"MEEJJ","MEE","l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","nleptons","njets","nvertices","dR_l1l2","dR_j1j2","dR_l1j1","dR_l1j2","dR_l2j1","dR_l2j2","MJJ","MET"};


	for(int icanvas=0; icanvas<nhistos; icanvas++){
#ifdef DEBUG
		std::cout<<"drawing histos on canvas"<<std::endl;
#endif
		TString name = "";
		name += icanvas;
		TCanvas* mycanvas = new TCanvas( name, name, 800, 800 ) ;
#ifdef DEBUG
		std::cout<<"made TCanvas and initialized a pointer to it"<<std::endl;
#endif
		mycanvas->cd();

		///rescale vertical axis for linear scale plot
		Double_t oldMax = ( (ths[icanvas]->GetMaximum() > histos[numProcesses-1][icanvas]->GetMaximum() ) ? ths[icanvas]->GetMaximum() : histos[numProcesses-1][icanvas]->GetMaximum() );
		ths[icanvas]->SetMaximum(1.5*oldMax);
		histos[numProcesses-1][icanvas]->SetMaximum(1.5*oldMax);
		
		ths[icanvas]->Draw("hist");
		histos[numProcesses-1][icanvas]->Draw("epsame");
		mycanvas->Update();

		///set the x axis title
		ths[icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
		histos[numProcesses-1][icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
	
		///set the y axis title using the bin size from the data histo, and the x axis title
		TString yLabel = DetermineYaxisName(histos[numProcesses-1][icanvas], xtitles[icanvas]);
		ths[icanvas]->GetYaxis()->SetTitle(yLabel );
		histos[numProcesses-1][icanvas]->GetYaxis()->SetTitle(yLabel );

		///set the histogram title to show a name
		ths[icanvas]->SetTitle(titles[icanvas]);
		histos[numProcesses-1][icanvas]->SetTitle(titles[icanvas]);

		///draw the legend, and a text box above the legend with the COM energy and integrated lumi
		leg->Draw();
		
		mycanvas->Update();
		TString tag = "_25ns";
		
		TString fn = "tempPlots/electrons25ns/";
		TString fn_pdf = fn + fnames[icanvas].Data() + tag + ".pdf";
		TString fn_png = fn + fnames[icanvas].Data() + tag + ".png";
		mycanvas->Print(fn_pdf.Data());
		mycanvas->Print(fn_png.Data());
		mycanvas->SetLogy();

		///rescale vertical axis for log scale plot
		ths[icanvas]->SetMaximum(27*oldMax);
		histos[numProcesses-1][icanvas]->SetMaximum(27*oldMax);
		ths[icanvas]->SetMinimum(0.2);
		histos[numProcesses-1][icanvas]->SetMinimum(0.2);

		leg->Draw();
		mycanvas->Update();
		TString fn_log_pdf = fn + fnames[icanvas].Data() + tag + "_log.pdf";
		TString fn_log_png = fn + fnames[icanvas].Data() + tag + "_log.png";
		mycanvas->Print(fn_log_pdf.Data());
		mycanvas->Print(fn_log_png.Data());
		mycanvas->Close();

	}
}


vector<TH1F*> MakeNHistos(TString hname, int n, int bins, float x_min, float x_max){
#ifdef DEBUG
	std::cout<<"in MakeNHistos"<<std::endl;
#endif

	vector<TH1F*> NHistos;
	for(int i=0; i<n; i++)
	{
		TString full_name = Form("%s_%d",hname.Data(),i);      
		TH1F *h = new TH1F(full_name.Data(),"",bins,x_min,x_max);
		NHistos.push_back(h);
	}

#ifdef DEBUG
	std::cout<<"leaving MakeNHistos"<<std::endl;
#endif

	return NHistos;
}

/**
 *
 * determine the name for the y axis using a pointer to a TH1F histogram made with real data
 * and the label applied to the x axis 
 */
TString DetermineYaxisName(TH1F * ptrDataHist, TString xLabel){
	char tempTitle[120];
	Float_t binWidth = ptrDataHist->GetXaxis()->GetBinWidth(1);
	if(xLabel.Contains("GeV")){
		if(binWidth >= 0.1) sprintf(tempTitle,"Events / %.1f / GeV", binWidth);
		else sprintf(tempTitle,"Events / %.3f / GeV", binWidth);
	}///end if xLabel contains GeV
	else{
		///xLabel does not contain GeV
		if(binWidth >= 0.1) sprintf(tempTitle,"Events / %.1f", binWidth);
		else sprintf(tempTitle,"Events / %.3f", binWidth);
	}///end else
	TString yTitle(tempTitle);
	return yTitle;

}///end DetermineYaxisName()

void Fill_Histo(std::vector<TH1F*> h1, TTree* tree, std::vector<float> PUW, bool pileup_reweight, bool is_data){  
#ifdef DEBUG
	std::cout<<"in Fill_Histo"<<std::endl;
#endif

	int nentries = tree->GetEntries();
	Float_t fourObjectMass;
	Float_t ptEle[2];
	Float_t etaEle[2];
	Float_t phiEle[2];
	Float_t ptJet[2];
	Float_t etaJet[2];
	Float_t phiJet[2];
	Float_t dileptonMass;
	Float_t dijetMass;
	Float_t missingET;
	
	Float_t dR_l1l2;
	Float_t dR_j1j2;
	Float_t dR_l1j1;
	Float_t dR_l1j2;
	Float_t dR_l2j1;
	Float_t dR_l2j2;
	Int_t nLeptons;
	Int_t nJets;
	Int_t nVertices;
	Float_t evWeightSign;

	tree->SetBranchAddress("fourObjectMass",&fourObjectMass);
	tree->SetBranchAddress("ptEle",ptEle);
	tree->SetBranchAddress("etaEle",etaEle);
	tree->SetBranchAddress("phiEle",phiEle);

	tree->SetBranchAddress("ptJet",ptJet);
	tree->SetBranchAddress("etaJet",etaJet);
	tree->SetBranchAddress("phiJet",phiJet);

	tree->SetBranchAddress("dileptonMass",&dileptonMass);
	tree->SetBranchAddress("dijetMass",&dijetMass);
	tree->SetBranchAddress("missingET",&missingET);
	
	tree->SetBranchAddress("dR_leadingLeptonSubleadingLepton",&dR_l1l2);
	tree->SetBranchAddress("dR_leadingJetSubleadingJet",&dR_j1j2);
	tree->SetBranchAddress("dR_leadingLeptonLeadingJet",&dR_l1j1);
	tree->SetBranchAddress("dR_leadingLeptonSubleadingJet",&dR_l1j2);
	tree->SetBranchAddress("dR_subleadingLeptonLeadingJet",&dR_l2j1);
	tree->SetBranchAddress("dR_subleadingLeptonSubleadingJet",&dR_l2j2);

	tree->SetBranchAddress("nLeptons",&nLeptons);
	tree->SetBranchAddress("nJets",&nJets);
	tree->SetBranchAddress("nVertices",&nVertices);

	if(!is_data)
		tree->SetBranchAddress("evWeightSign",&evWeightSign);

	for (Int_t ev = 0; ev < nentries; ev++) {
		Float_t reweight = 1;
		tree->GetEntry(ev); 

		if(dileptonMass > 50)//l1_pt>60 && l2_pt>50 && j1_pt>40 && j2_pt>40 && dileptonMass<200)// && (dR_l1l2 < 2.5 || dR_l1l2 > 3.5))// && dR_l1j1 > 0.4 && dR_l1j2 > 0.4 && dR_l2j1 > 0.4 && dR_l2j2 > 0.4)
		{
			if(!is_data){
				if(pileup_reweight)
					reweight = evWeightSign*PUW[nVertices];
				else
					reweight = evWeightSign;

				h1[0]->Fill(fourObjectMass,reweight);
				h1[1]->Fill(dileptonMass,reweight);
				h1[2]->Fill(ptEle[0],reweight);
				h1[3]->Fill(ptEle[1],reweight);
				h1[4]->Fill(ptJet[0],reweight);
				h1[5]->Fill(ptJet[1],reweight);
				h1[6]->Fill(etaEle[0],reweight);
				h1[7]->Fill(etaEle[1],reweight);
				h1[8]->Fill(etaJet[0],reweight);
				h1[9]->Fill(etaJet[1],reweight);
				h1[10]->Fill(phiEle[0],reweight);
				h1[11]->Fill(phiEle[1],reweight);
				h1[12]->Fill(phiJet[0],reweight);
				h1[13]->Fill(phiJet[1],reweight);	
				h1[14]->Fill(nLeptons,reweight);	
				h1[15]->Fill(nJets,reweight);	
				h1[16]->Fill(nVertices,reweight);	
				h1[17]->Fill(dR_l1l2,reweight);
				h1[18]->Fill(dR_j1j2,reweight);
				h1[19]->Fill(dR_l1j1,reweight);
				h1[20]->Fill(dR_l1j2,reweight);
				h1[21]->Fill(dR_l2j1,reweight);
				h1[22]->Fill(dR_l2j2,reweight);
				h1[23]->Fill(dijetMass,reweight);
				h1[24]->Fill(missingET,reweight);
			}///end if(!is_data)
			else {
				h1[0]->Fill(fourObjectMass);
				h1[1]->Fill(dileptonMass);
				h1[2]->Fill(ptEle[0]);
				h1[3]->Fill(ptEle[1]);
				h1[4]->Fill(ptJet[0]);
				h1[5]->Fill(ptJet[1]);
				h1[6]->Fill(etaEle[0]);
				h1[7]->Fill(etaEle[1]);
				h1[8]->Fill(etaJet[0]);
				h1[9]->Fill(etaJet[1]);
				h1[10]->Fill(phiEle[0]);
				h1[11]->Fill(phiEle[1]);
				h1[12]->Fill(phiJet[0]);
				h1[13]->Fill(phiJet[1]);	
				h1[14]->Fill(nLeptons);	
				h1[15]->Fill(nJets);	
				h1[16]->Fill(nVertices);	
				h1[17]->Fill(dR_l1l2);
				h1[18]->Fill(dR_j1j2);
				h1[19]->Fill(dR_l1j1);
				h1[20]->Fill(dR_l1j2);
				h1[21]->Fill(dR_l2j1);
				h1[22]->Fill(dR_l2j2);
				h1[23]->Fill(dijetMass);
				h1[24]->Fill(missingET);
			}///end else (for real data)
		}///end if(true)
	}///end loop over events ev in tree

#ifdef DEBUG
	std::cout<<"leaving Fill_Histo"<<std::endl;
#endif
}///end Fill_Histo()

std::vector<float> PileUpWeights(TTree* tree,TTree* tree_data){
	TH1F *h_data = new TH1F("h_data","",50,0,50);
	TH1F *h0 = new TH1F("h0","",50,0,50);  

	Int_t nVertices;
	Int_t nVertices_data;
	tree->SetBranchAddress("nVertices",&nVertices);
	tree_data->SetBranchAddress("nVertices",&nVertices_data);

	for (Int_t ev = 0; ev < tree->GetEntries(); ev++) {
		tree->GetEntry(ev);
		h0->Fill(nVertices);
	}
	for (Int_t ev = 0; ev < tree_data->GetEntries(); ev++) {
		tree_data->GetEntry(ev);
		h_data->Fill(nVertices_data);
	}

	h_data->Scale(1/h_data->Integral());
	h0->Scale(1/h0->Integral());

	std::vector<float> PUW(51);

	for(int i = 0;i<51;i++)
		if(h0->GetBinContent(i) != 0)
			PUW[i] = h_data->GetBinContent(i)/h0->GetBinContent(i);
		else
			PUW[i] = 1.0;

	h0->Delete();
	h_data->Delete();

	return PUW;

}
