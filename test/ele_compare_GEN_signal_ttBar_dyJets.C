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
//#define CHKTTBAR
//#define printBkgndIntegrals

TString DetermineYaxisName(TH1F * ptrDataHist, TString xLabel);
void Fill_Histo(std::vector<TH1F*> h1, TTree* tree, std::vector<float> PUW, bool pileup_reweight, bool is_data);
std::vector<TH1F*> MakeNHistos(TString hname, int n, int bins, float x_min, float x_max);
void Fill_TwoDimHisto(std::vector<TH2F*> h1, TTree* tree, std::vector<float> PUW, bool pileup_reweight, bool is_data);
std::vector<TH2F*> MakeTwoNHistos(TString hname, int n, int binsX, float x_min, float x_max, int binsY, float y_min, float y_max);
std::vector<float> PileUpWeights(TTree* tree,TTree* tree_data);

void ele_signal_and_bkgndMC_compare(){

#ifdef DEBUG
	std::cout<<"in ele_dataMC_compare()"<<std::endl;
#endif
	
	TString directory = "/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/";
	TFile * hfile0 = new TFile(directory+"all_genDYJetsAMCNLO_M200to400_DecayKinematics.root");//dyjets MLL 200 to 400
	TFile * hfile1 = new TFile(directory+"complete_genTTBarDecayKinematics.root");//ttbar
	TFile * hfile3 = new TFile(directory+"all_genDYJetsAMCNLO_M400to500_DecayKinematics.root");//dyjets MLL 400 to 500
	TFile * hfile4 = new TFile(directory+"all_genDYJetsAMCNLO_M500to700_DecayKinematics.root");//dyjets MLL 500 to 700
	TFile * hfile5 = new TFile(directory+"all_genDYJetsAMCNLO_M700to800_DecayKinematics.root");//dyjets MLL 700 to 800

	TFile * hfile_signal0 = new TFile(directory+"all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-1400_Nu_M-700.root");
	TFile * hfile_signal1 = new TFile(directory+"all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-2000_Nu_M-1000.root");
	//TFile * hfile_signal2 = new TFile(directory+"wr3200nu1600Tree.root");

	int numDiffProcesses = 7;
#ifdef DEBUG
	std::cout<<"declared pointers to input files"<<std::endl;
#endif

	TString signalTreeName = "genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts";
	TString bkgndTreeName = "genMatchedAnalyzerThree/genLeptonsAndJetsAfterAllCuts";
	TTree *tree0 = (TTree*)hfile0->Get(bkgndTreeName);  
	TTree *tree1 = (TTree*)hfile1->Get(bkgndTreeName);  
	TTree *tree3 = (TTree*)hfile3->Get(bkgndTreeName);  
	TTree *tree4 = (TTree*)hfile4->Get(bkgndTreeName);  
	TTree *tree5 = (TTree*)hfile5->Get(bkgndTreeName);  

	TTree *tree_signal0 = (TTree*)hfile_signal0->Get(signalTreeName);  
	TTree *tree_signal1 = (TTree*)hfile_signal1->Get(signalTreeName);  
	//TTree *tree_signal2 = (TTree*)hfile_signal2->Get(signalTreeName);  

#ifdef DEBUG
	std::cout<<"declared pointers to trees in input files"<<std::endl;
#endif

	///set histo names, number of different physics processes, number of bins, and axis limits here
	vector<TH1F*> h_Mlljj = MakeNHistos("h_Mlljj",numDiffProcesses,20,600,2400);
	vector<TH1F*> h_Mll = MakeNHistos("h_Mll",numDiffProcesses,22,200,2000);
	vector<TH1F*> h_l1pt = MakeNHistos("h_l1pt",numDiffProcesses,22,0,1800);
	vector<TH1F*> h_l2pt = MakeNHistos("h_l2pt",numDiffProcesses,14,0,1400);
	vector<TH1F*> h_j1pt = MakeNHistos("h_j1pt",numDiffProcesses,22,0,1800);
	vector<TH1F*> h_j2pt = MakeNHistos("h_j2pt",numDiffProcesses,14,0,1400);
	vector<TH1F*> h_l1eta = MakeNHistos("h_l1eta", numDiffProcesses,20,-3.,3.);
	vector<TH1F*> h_l2eta = MakeNHistos("h_l2eta", numDiffProcesses,20,-3.,3.);
	vector<TH1F*> h_j1eta = MakeNHistos("h_j1eta", numDiffProcesses,20,-3.,3.);
	vector<TH1F*> h_j2eta = MakeNHistos("h_j2eta", numDiffProcesses,20,-3.,3.);
	vector<TH1F*> h_l1phi = MakeNHistos("h_l1phi", numDiffProcesses,20,-3.15,3.15);
	vector<TH1F*> h_l2phi = MakeNHistos("h_l2phi", numDiffProcesses,20,-3.15,3.15);
	vector<TH1F*> h_j1phi = MakeNHistos("h_j1phi", numDiffProcesses,20,-3.15,3.15);
	vector<TH1F*> h_j2phi = MakeNHistos("h_j2phi", numDiffProcesses,20,-3.15,3.15);
	vector<TH1F*> h_nleptons = MakeNHistos("h_nleptons", numDiffProcesses,4,0,4);
	vector<TH1F*> h_njets = MakeNHistos("h_njets", numDiffProcesses,6,0,6);
	vector<TH1F*> h_nvertices = MakeNHistos("h_nvertices", numDiffProcesses,33,0,33);
	vector<TH1F*> h_dR_l1l2 = MakeNHistos("h_dR_l1l2", numDiffProcesses,20,0,5);
	vector<TH1F*> h_dR_j1j2 = MakeNHistos("h_dR_j1j2", numDiffProcesses,20,0,5);
	vector<TH1F*> h_dR_l1j1 = MakeNHistos("h_dR_l1j1", numDiffProcesses,20,0,5);
	vector<TH1F*> h_dR_l1j2 = MakeNHistos("h_dR_l1j2", numDiffProcesses,20,0,5);
	vector<TH1F*> h_dR_l2j1 = MakeNHistos("h_dR_l2j1", numDiffProcesses,20,0,5);
	vector<TH1F*> h_dR_l2j2 = MakeNHistos("h_dR_l2j2", numDiffProcesses,20,0,5);
	vector<TH1F*> h_Mjj = MakeNHistos("h_Mjj",numDiffProcesses,40,0,1700);
	vector<TH1F*> h_MET = MakeNHistos("h_MET",numDiffProcesses,20,0,400);
	vector<TH1F*> h_Msubleadljj = MakeNHistos("h_Msubleadljj",numDiffProcesses,20,0,1400);
	vector<TH1F*> h_Mleadljj = MakeNHistos("h_Mleadljj",numDiffProcesses,20,0,1800);
	
	
	vector<TH2F*> h_Mljj_Mlljj = MakeTwoNHistos("h_Mljj_Mlljj",numDiffProcesses,15,600,2600,20,50,1600);	///x axis bins and range first, then y axis bins and range
	
#ifdef DEBUG
	std::cout<<"made vectors of TH1F and TH2F pointers to histos"<<std::endl;
#endif

	int nhistos = 1;	//max is 27
	int nTwoDimHistos = 1;	//max is 1
	std::vector<THStack*> ths, std::vector<THStack*> thsTwoDim; // Stacks;
	for(int istacks=0; istacks<nhistos; istacks++)
	{
		TString full_name = Form("th%d",istacks);      
		THStack *th = new THStack(full_name.Data(),"");
		ths.push_back(th);
	}
	for(int istacks=0; istacks<nTwoDimHistos; istacks++)
	{
		TString full_name = Form("th%d",istacks);      
		THStack *thtwo = new THStack(full_name.Data(),"");
		thsTwoDim.push_back(thtwo);
	}
	
#ifdef DEBUG
	std::cout<<"made vector of THStack pointers"<<std::endl;
#endif

	///declare one vector of TH1F pointers for each MC process and real data sample
	///the number of elements in each vector should be equal to the number of unique vector<TH1F*> containers which were declared above 
	std::vector<TH1F*> histos0(nhistos); // DY M200to400
	std::vector<TH1F*> histos1(nhistos); // TTbar
	std::vector<TH1F*> histos3(nhistos); // DY M400to500
	std::vector<TH1F*> histos4(nhistos); // DY M500to700
	std::vector<TH1F*> histos5(nhistos); // DY M700to800
	
	std::vector<TH1F*> histos_signal0(nhistos);
	std::vector<TH1F*> histos_signal1(nhistos);
	//std::vector<TH1F*> histos_signal2(nhistos);

	///declare one vector of TH2F pointers for each process which will be plotted
	std::vector<TH2F*> histosTwoDim0(nTwoDimHistos); // DY M200to400
	std::vector<TH2F*> histosTwoDim1(nTwoDimHistos); // TTbar
	std::vector<TH2F*> histosTwoDim3(nTwoDimHistos); // DY M400to500
	std::vector<TH2F*> histosTwoDim4(nTwoDimHistos); // DY M500to700
	std::vector<TH2F*> histosTwoDim5(nTwoDimHistos); // DY M700to800
	
	std::vector<TH2F*> histosTwoDim_signal0(nTwoDimHistos);
	std::vector<TH2F*> histosTwoDim_signal1(nTwoDimHistos);
	//std::vector<TH2F*> histosTwoDim_signal2(nTwoDimHistos);



	///load the vectors tied to all MC and real data samples into one vector<vector<TH1F*> > object
	vector<vector<TH1F*> > histos(numDiffProcesses);
	histos[0] = histos0;
	histos[1] = histos1;
	histos[2] = histos3;
	histos[3] = histos4;
	histos[4] = histos5;
	histos[5] = histos_signal0;
	histos[6] = histos_signal1;
	//histos[7] = histos_signal2;

	vector<vector<TH2F*> > histosTwoDim(numDiffProcesses);
	histosTwoDim[0] = histosTwoDim0;
	histosTwoDim[1] = histosTwoDim1;
	histosTwoDim[2] = histosTwoDim3;
	histosTwoDim[3] = histosTwoDim4;
	histosTwoDim[4] = histosTwoDim5;
	histosTwoDim[5] = histosTwoDim_signal0;
	histosTwoDim[6] = histosTwoDim_signal1;
	//histosTwoDim[7] = histosTwoDim_signal2;


	
	///link the histos made by calls to MakeNHistos() above to specific elements of vector<vector<TH1F*> > container 
	for(int j=0; j<numDiffProcesses; j++){
		//one dimensional histos
		histos[j][0] = h_Mlljj[j];
		//histos[j][1] = h_Mll[j];
		//histos[j][2] = h_l1pt[j];
		//histos[j][3] = h_l2pt[j];
		//histos[j][4] = h_j1pt[j];
		//histos[j][5] = h_j2pt[j];
		//histos[j][6] = h_l1eta[j];
		//histos[j][7] = h_l2eta[j];
		//histos[j][8] = h_j1eta[j];
		//histos[j][9] = h_j2eta[j];
		//histos[j][10] = h_l1phi[j];
		//histos[j][11] = h_l2phi[j];
		//histos[j][12] = h_j1phi[j];
		//histos[j][13] = h_j2phi[j];
		//histos[j][14] = h_nleptons[j];
		//histos[j][15] = h_njets[j];
		//histos[j][16] = h_nvertices[j];
		//histos[j][17] = h_dR_l1l2[j];
		//histos[j][18] = h_dR_j1j2[j];
		//histos[j][19] = h_dR_l1j1[j];
		//histos[j][20] = h_dR_l1j2[j];
		//histos[j][21] = h_dR_l2j1[j];
		//histos[j][22] = h_dR_l2j2[j];
		//histos[j][23] = h_Mjj[j];
		//histos[j][24] = h_MET[j];
		//histos[j][25] = h_Msubleadljj[j];
		//histos[j][26] = h_Mleadljj[j];

		//two dimensional histos
		histosTwoDim[j][0] = h_Mljj_Mlljj[j];

	}
	
#ifdef DEBUG
	std::cout<<"linked MC processes to histos of specific quantities"<<std::endl;
#endif

	///calculate PU weights and store them in persistent vectors of floats
	std::vector<float> PUW0 = PileUpWeights(tree0,tree_signal0);
	std::vector<float> PUW1 = PileUpWeights(tree1,tree_signal0);
	std::vector<float> PUW3 = PileUpWeights(tree3,tree_signal0);
	std::vector<float> PUW4 = PileUpWeights(tree4,tree_signal0);
	std::vector<float> PUW5 = PileUpWeights(tree5,tree_signal0);
	
	std::vector<float> PUW_signal0 = PileUpWeights(tree_signal0,tree_signal0);
	std::vector<float> PUW_signal1 = PileUpWeights(tree_signal1,tree_signal0);
	//std::vector<float> PUW_signal2 = PileUpWeights(tree_signal2,tree_signal0);

	///fill the 1D and 2D histos with content from TTrees
	Fill_Histo(histos[0],tree0,PUW0,false,false); // DY
	Fill_Histo(histos[1],tree1,PUW1,false,false); // TTbar
	Fill_Histo(histos[2],tree3,PUW3,false,false); // WZ
	Fill_Histo(histos[3],tree4,PUW4,false,false); // ZZ
	Fill_Histo(histos[4],tree5,PUW5,false,false); // WJets

	Fill_Histo(histos[5],tree_signal0,PUW_signal0,false,false);
	Fill_Histo(histos[6],tree_signal1,PUW_signal1,false,false);
	//Fill_Histo(histos[7],tree_signal2,PUW_signal2,false,false);

	Fill_TwoDimHisto(histosTwoDim[0],tree0,PUW0,false,false);
	Fill_TwoDimHisto(histosTwoDim[1],tree1,PUW1,false,false);
	Fill_TwoDimHisto(histosTwoDim[2],tree3,PUW3,false,false);
	Fill_TwoDimHisto(histosTwoDim[3],tree4,PUW4,false,false);
	Fill_TwoDimHisto(histosTwoDim[4],tree5,PUW5,false,false);
	Fill_TwoDimHisto(histosTwoDim[5],tree_signal0,PUW_signal0,false,false);
	Fill_TwoDimHisto(histosTwoDim[6],tree_signal1,PUW_signal1,false,false);
	//Fill_TwoDimHisto(histosTwoDim[7],tree_signal2,PUW_signal2,false,false);
	
	Float_t intLumi = 1500.;
	///rescale the 1D histos
	for(std::vector<TH1F*>::size_type i = 0; i != nhistos; i++){
#ifdef DEBUG
		std::cout<<"setting line and fill colors, and adding histos to THStack objects"<<std::endl;
#endif
		//histos[0][i]->Scale(6025.2*(intLumi)/9042031);	///madgraph MLL 50 DYJetsToLL 25ns
		//histos[0][i]->Scale(6025.2*(intLumi)/28747969);	///amcatnlo MLL 50 DYJetsToLL 25ns
		histos[0][i]->Scale(7.67*(intLumi)/96310);	///amcatnlo MLL 200to400 DYJetsToLL 25ns
		histos[0][i]->SetFillColor(5);
		
#ifdef printBkgndIntegrals
		cout<<"DYJets MLL 200to400 contributes\t"<< histos[0][i]->Integral() <<endl;
#endif
	
		histos[1][i]->Scale(831.76*(intLumi)/96834559);	///powheg-pythia ttBar to all 25ns
		histos[1][i]->SetFillColor(3);
#ifdef printBkgndIntegrals
		cout<<"TTBar contributes\t"<< histos[1][i]->Integral() <<endl;
#endif

#ifdef CHKTTBAR
		std::cout<<"ee chnl ttBar integral =\t"<< histos[1][i]->Integral() <<std::endl;
#endif
	
		//histos[2][i]->Scale(66.1*(intLumi)/991232);	///WZ to all 25ns
		histos[2][i]->Scale(0.423*(intLumi)/98441);	///amcnlo DYJetsToLL MLL 400to500 25ns
		histos[2][i]->SetFillColor(4);
#ifdef printBkgndIntegrals
		cout<<"DYJets MLL 400to500 contributes\t"<< histos[2][i]->Integral() <<endl;
#endif

		//histos[3][i]->Scale(15.4*(intLumi)/996168);		///ZZ to all 25ns
		histos[3][i]->Scale(0.24*(intLumi)/101029);		///amcnlo DYJetsToLL MLL 500to700 25ns
		histos[3][i]->SetFillColor(7);
#ifdef printBkgndIntegrals
		cout<<"DYJets MLL 500to700 contributes\t"<< histos[3][i]->Integral() <<endl;
#endif
	
		//histos[4][i]->Scale(6.15e4*(intLumi)/24151270);		///WJetsToLNu 25ns
		histos[4][i]->Scale(0.035*(intLumi)/96011);		///amcnlo DYJetsToLL MLL 700to800 25ns
		histos[4][i]->SetFillColor(6);
#ifdef printBkgndIntegrals
		cout<<"DYJets MLL 700to800 contributes\t"<< histos[4][i]->Integral() <<endl;
#endif
	
		//histos[5][i]->Scale(1.78*intLumi/50000);
		//histos[5][i]->SetMarkerStyle(20);

		histos[5][i]->Scale(0.389*intLumi/50000);
		histos[5][i]->SetLineStyle(1);
		histos[5][i]->SetLineColor(1);
		histos[5][i]->SetLineWidth(3);

		histos[6][i]->Scale(0.0707*intLumi/50000);
		histos[6][i]->SetLineStyle(1);
		histos[6][i]->SetLineColor(2);
		histos[6][i]->SetLineWidth(3);


		//histos[7][i]->Scale(0.0034*intLumi/50000);
		//histos[7][i]->SetLineStyle(1);
		//histos[7][i]->SetLineColor(4);
		//histos[7][i]->SetLineWidth(3);

		ths[i]->Add(histos[4][i]);
		ths[i]->Add(histos[3][i]);
		ths[i]->Add(histos[2][i]);
		ths[i]->Add(histos[0][i]);
		ths[i]->Add(histos[1][i]);
	
	}//end rescaling of 1D histos

	///rescale the 2D histos
	for(std::vector<TH2F*>::size_type i = 0; i != nTwoDimHistos; i++){
		histosTwoDim[0][i]->Scale(7.67*(intLumi)/96310);	///amcatnlo MLL 200to400 DYJetsToLL 25ns
		histosTwoDim[0][i]->SetFillColor(5);

		histosTwoDim[1][i]->Scale(831.76*(intLumi)/96834559);	///powheg-pythia ttBar to all 25ns
		histosTwoDim[1][i]->SetFillColor(3);

		histosTwoDim[2][i]->Scale(0.423*(intLumi)/98441);	///amcnlo DYJetsToLL MLL 400to500 25ns
		histosTwoDim[2][i]->SetFillColor(4);
		
		histosTwoDim[3][i]->Scale(0.24*(intLumi)/101029);		///amcnlo DYJetsToLL MLL 500to700 25ns
		histosTwoDim[3][i]->SetFillColor(7);

		histosTwoDim[4][i]->Scale(0.035*(intLumi)/96011);		///amcnlo DYJetsToLL MLL 700to800 25ns
		histosTwoDim[4][i]->SetFillColor(6);

		histosTwoDim[5][i]->Scale(0.389*intLumi/50000);
		histosTwoDim[5][i]->SetLineStyle(1);
		histosTwoDim[5][i]->SetLineColor(1);
		histosTwoDim[5][i]->SetLineWidth(3);

		histosTwoDim[6][i]->Scale(0.0707*intLumi/50000);
		histosTwoDim[6][i]->SetLineStyle(1);
		histosTwoDim[6][i]->SetLineColor(2);
		histosTwoDim[6][i]->SetLineWidth(3);

		thsTwoDim[i]->Add(histosTwoDim[4][i]);
		thsTwoDim[i]->Add(histosTwoDim[3][i]);
		thsTwoDim[i]->Add(histosTwoDim[2][i]);
		thsTwoDim[i]->Add(histosTwoDim[0][i]);
		thsTwoDim[i]->Add(histosTwoDim[1][i]);
	
	}//end rescaling of 2D histos

	///make a legend with appropriate labels for MC processes
	TLegend *leg = new TLegend( 0.62, 0.64, 0.89, 0.87 );
	leg->SetNColumns(2);
	//leg->AddEntry( histos[5][0], "M_{WR}=1 M_{NU}=0.5 TeV" ,"ep");
	leg->AddEntry( histos[5][0], "M_{WR}=1.4 TeV","l");
	leg->AddEntry( histos[6][0], "M_{WR}=2 TeV","l");
	leg->AddEntry( histos[0][0], "DY MLL 200to400" ); 
	leg->AddEntry( histos[1][0], "ttbar" );
	leg->AddEntry( histos[2][0], "DY MLL 400to500" );  
	leg->AddEntry( histos[3][0], "DY MLL 500to700" ); 
	leg->AddEntry( histos[4][0], "DY MLL 700to800" );
	leg->SetFillColor( kWhite );

	TString xtitles[] = {"M_{EEJJ} [GeV]","M_{EE} [GeV]","leading electron p_{T} [GeV]","subleading electron p_{T} [GeV]","leading jet p_{T} [GeV]","subleading jet p_{T} [GeV]","leading electron #eta","subleading electron #eta","leading jet #eta","subleading jet #eta","leading electron #phi","subleading electron #phi","leading jet #phi","subleading jet #phi","number of electrons","number of jets","number of vertices","#DeltaR lead ele sublead ele","#DeltaR lead jet sublead jet","#DeltaR lead ele lead jet","#DeltaR lead ele sublead jet","#DeltaR sublead ele lead jet","#DeltaR sublead ele sublead jet","M_{JJ} [GeV]","Missing E_{T} [GeV]"};
	
	TString titles[] = {"CMS Preliminary M_{EEJJ}  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary DiElectron Mass  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Lead Electron p_{T}  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Sublead Electron p_{T}  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Lead Jet p_{T}  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Sublead Jet p_{T}  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Lead Electron #eta  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Sublead Electron #eta  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Lead Jet #eta  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Subleading jet #eta  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary leading electron #phi  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Subleading electron #phi  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary leading jet #phi  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Subleading jet #phi  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary number of electrons  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary number of jets  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary number of vertices  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR lead ele Sublead ele  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR lead jet Sublead jet  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR lead ele lead jet  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR lead ele Sublead jet  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR Sublead ele lead jet  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary #DeltaR Sublead ele Sublead jet  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Dijet Mass  #surds = 13 TeV   #intlumi = 1.5/fb","CMS Preliminary Missing E_{T}  #surds = 13 TeV   #intlumi = 1.5/fb"};

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
		Double_t oldMax = ( (ths[icanvas]->GetMaximum() > histos[5][icanvas]->GetMaximum() ) ? ths[icanvas]->GetMaximum() : histos[5][icanvas]->GetMaximum() );
		ths[icanvas]->SetMaximum(1.5*oldMax);
		histos[5][icanvas]->SetMaximum(1.5*oldMax);
		//histos[6][icanvas]->SetMaximum(1.5*oldMax);
		//histos[7][icanvas]->SetMaximum(1.5*oldMax);
	
		ths[icanvas]->Draw("hist");
		histos[5][icanvas]->Draw("Csame");
		//histos[6][icanvas]->Draw("Csame");
		//histos[7][icanvas]->Draw("lsame");
		mycanvas->Update();

		///set the x axis title
		ths[icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
		histos[5][icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
		//histos[6][icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
		//histos[7][icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
	
		///set the y axis title using the bin size from the data histo, and the x axis title
		TString yLabel = DetermineYaxisName(histos[5][icanvas], xtitles[icanvas]);
		ths[icanvas]->GetYaxis()->SetTitle(yLabel );
		histos[5][icanvas]->GetYaxis()->SetTitle(yLabel );
		//histos[6][icanvas]->GetYaxis()->SetTitle(yLabel );
		//histos[7][icanvas]->GetYaxis()->SetTitle(yLabel );

		///set the histogram title to show a name
		ths[icanvas]->SetTitle(titles[icanvas]);
		histos[5][icanvas]->SetTitle(titles[icanvas]);
		//histos[6][icanvas]->SetTitle(titles[icanvas]);
		//histos[7][icanvas]->SetTitle(titles[icanvas]);

		///draw the legend, and a text box above the legend with the COM energy and integrated lumi
		leg->Draw();
		
		mycanvas->Update();
		TString tag = "_amcnloDYJets_25ns";
		
		TString fn = "tempPlots/wrSignalAndBkgndMCs25ns/";
		TString fn_pdf = fn + fnames[icanvas].Data() + tag + ".pdf";
		TString fn_png = fn + fnames[icanvas].Data() + tag + ".png";
		mycanvas->Print(fn_pdf.Data());
		mycanvas->Print(fn_png.Data());
		mycanvas->SetLogy();

		///rescale vertical axis for log scale plot
		ths[icanvas]->SetMaximum(27*oldMax);
		histos[5][icanvas]->SetMaximum(27*oldMax);
		//histos[6][icanvas]->SetMaximum(27*oldMax);
		//histos[7][icanvas]->SetMaximum(27*oldMax);
		ths[icanvas]->SetMinimum(0.2);
		histos[5][icanvas]->SetMinimum(0.2);
		//histos[6][icanvas]->SetMinimum(0.2);
		//histos[7][icanvas]->SetMinimum(0.2);

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

vector<TH2F*> MakeTwoNHistos(TString hname, int n, int binsX, float x_min, float x_max, int binsY, float y_min, float y_max){
#ifdef DEBUG
	std::cout<<"in MakeTwoNHistos"<<std::endl;
#endif

	vector<TH2F*> NHistos;
	for(int i=0; i<n; i++)
	{
		TString full_name = Form("%s_%d",hname.Data(),i);      
		TH2F *h = new TH2F(full_name.Data(),"",binsX,x_min,x_max,binsY,y_min,y_max);
		NHistos.push_back(h);
	}

#ifdef DEBUG
	std::cout<<"leaving MakeTwoNHistos"<<std::endl;
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
	Float_t subleadingLeptonThreeObjMass;
	Float_t leadLeptonThreeObjMass;
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
	tree->SetBranchAddress("subleadingLeptonThreeObjMass",&subleadingLeptonThreeObjMass);
	tree->SetBranchAddress("leadLeptonThreeObjMass",&leadLeptonThreeObjMass);
	
	tree->SetBranchAddress("missingET",&missingET);
	tree->SetBranchAddress("ptEle",ptEle);
	tree->SetBranchAddress("etaEle",etaEle);
	tree->SetBranchAddress("phiEle",phiEle);

	tree->SetBranchAddress("ptJet",ptJet);
	tree->SetBranchAddress("etaJet",etaJet);
	tree->SetBranchAddress("phiJet",phiJet);

	tree->SetBranchAddress("dileptonMass",&dileptonMass);
	tree->SetBranchAddress("dijetMass",&dijetMass);
	
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
				//h1[1]->Fill(dileptonMass,reweight);
				//h1[2]->Fill(ptEle[0],reweight);
				//h1[3]->Fill(ptEle[1],reweight);
				//h1[4]->Fill(ptJet[0],reweight);
				//h1[5]->Fill(ptJet[1],reweight);
				//h1[6]->Fill(etaEle[0],reweight);
				//h1[7]->Fill(etaEle[1],reweight);
				//h1[8]->Fill(etaJet[0],reweight);
				//h1[9]->Fill(etaJet[1],reweight);
				//h1[10]->Fill(phiEle[0],reweight);
				//h1[11]->Fill(phiEle[1],reweight);
				//h1[12]->Fill(phiJet[0],reweight);
				//h1[13]->Fill(phiJet[1],reweight);	
				//h1[14]->Fill(nLeptons,reweight);	
				//h1[15]->Fill(nJets,reweight);	
				//h1[16]->Fill(nVertices,reweight);	
				//h1[17]->Fill(dR_l1l2,reweight);
				//h1[18]->Fill(dR_j1j2,reweight);
				//h1[19]->Fill(dR_l1j1,reweight);
				//h1[20]->Fill(dR_l1j2,reweight);
				//h1[21]->Fill(dR_l2j1,reweight);
				//h1[22]->Fill(dR_l2j2,reweight);
				//h1[23]->Fill(dijetMass,reweight);
				//h1[24]->Fill(missingET,reweight);
				//h1[25]->Fill(subleadingLeptonThreeObjMass,reweight);
				//h1[26]->Fill(leadLeptonThreeObjMass,reweight);
			}///end if(!is_data)
			else {
				h1[0]->Fill(fourObjectMass);
				//h1[1]->Fill(dileptonMass);
				//h1[2]->Fill(ptEle[0]);
				//h1[3]->Fill(ptEle[1]);
				//h1[4]->Fill(ptJet[0]);
				//h1[5]->Fill(ptJet[1]);
				//h1[6]->Fill(etaEle[0]);
				//h1[7]->Fill(etaEle[1]);
				//h1[8]->Fill(etaJet[0]);
				//h1[9]->Fill(etaJet[1]);
				//h1[10]->Fill(phiEle[0]);
				//h1[11]->Fill(phiEle[1]);
				//h1[12]->Fill(phiJet[0]);
				//h1[13]->Fill(phiJet[1]);	
				//h1[14]->Fill(nLeptons);	
				//h1[15]->Fill(nJets);	
				//h1[16]->Fill(nVertices);	
				//h1[17]->Fill(dR_l1l2);
				//h1[18]->Fill(dR_j1j2);
				//h1[19]->Fill(dR_l1j1);
				//h1[20]->Fill(dR_l1j2);
				//h1[21]->Fill(dR_l2j1);
				//h1[22]->Fill(dR_l2j2);
				//h1[23]->Fill(dijetMass);
				//h1[24]->Fill(missingET);
				//h1[25]->Fill(subleadingLeptonThreeObjMass);
				//h1[26]->Fill(leadLeptonThreeObjMass);
			}///end else (for real data)
		}///end if(true)
	}///end loop over events ev in tree

#ifdef DEBUG
	std::cout<<"leaving Fill_Histo"<<std::endl;
#endif
}///end Fill_Histo()

/**
 * this method is very similar to Fill_Histo(), but is designed for two dimensional histos
 */
void Fill_TwoDimHisto(std::vector<TH2F*> h1, TTree* tree, std::vector<float> PUW, bool pileup_reweight, bool is_data){
#ifdef DEBUG
	std::cout<<"in Fill_Histo"<<std::endl;
#endif

	int nentries = tree->GetEntries();
	Float_t fourObjectMass;
	Float_t subleadingLeptonThreeObjMass;
	Float_t leadLeptonThreeObjMass;
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
	tree->SetBranchAddress("subleadingLeptonThreeObjMass",&subleadingLeptonThreeObjMass);
	tree->SetBranchAddress("leadLeptonThreeObjMass",&leadLeptonThreeObjMass);
	
	tree->SetBranchAddress("missingET",&missingET);
	tree->SetBranchAddress("ptEle",ptEle);
	tree->SetBranchAddress("etaEle",etaEle);
	tree->SetBranchAddress("phiEle",phiEle);

	tree->SetBranchAddress("ptJet",ptJet);
	tree->SetBranchAddress("etaJet",etaJet);
	tree->SetBranchAddress("phiJet",phiJet);

	tree->SetBranchAddress("dileptonMass",&dileptonMass);
	tree->SetBranchAddress("dijetMass",&dijetMass);
	
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

		if(dileptonMass > 50)
		{
			if(!is_data){
				if(pileup_reweight)
					reweight = evWeightSign*PUW[nVertices];
				else
					reweight = evWeightSign;

				h1[0]->Fill(fourObjectMass,subleadingLeptonThreeObjMass,reweight);
			}///end if(!is_data)
			else {
				h1[0]->Fill(fourObjectMass,subleadingLeptonThreeObjMass);
			}///end else (for real data)
		}///end if(true)
	}///end loop over events ev in tree

#ifdef DEBUG
	std::cout<<"leaving Fill_TwoDimHisto"<<std::endl;
#endif
}///end Fill_TwoDimHisto()

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
