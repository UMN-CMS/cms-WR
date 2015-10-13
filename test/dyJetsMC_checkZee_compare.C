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

void dyJetsMC_checkZee_compare(){

#ifdef DEBUG
	std::cout<<"in ele_dataMC_compare()"<<std::endl;
#endif

	///swap btwn 50ns and 25ns in directory
	TString directory = "/eos/uscms/store/user/skalafut/analyzed_50ns_compare_DYJets_using_Zee_with_two_HEEP/";

	//TString fileTag = "";
	
	TFile * hfile0 = new TFile(directory+"analyzed_DYJets_AMCNLO_50ns_gen_and_reco_Zee.root");//dyjets
	//TFile * hfile1 = new TFile(directory+"analyzed_TTOnly_PowhegPythia_"+fileTag);//ttbar
	//TFile * hfile3 = new TFile(directory+"analyzed_WZ_"+fileTag);//wz
	//TFile * hfile4 = new TFile(directory+"analyzed_ZZ_"+fileTag);//zz
	//TFile * hfile5 = new TFile(directory+"analyzed_WJets_"+fileTag);//wjets

	TFile * hfile_data = new TFile(directory+"analyzed_DYJets_Madgraph_50ns_gen_and_reco_Zee.root");//data

#ifdef DEBUG
	std::cout<<"declared pointers to input files"<<std::endl;
#endif

	TString treeName = "genAnalyzerOne/genZedEleEle";
	//TString treeName = "recoAnalyzerOne/recoZedEleEle";
	TTree *tree0 = (TTree*)hfile0->Get(treeName);  
	//TTree *tree1 = (TTree*)hfile1->Get(treeName);  
	//TTree *tree3 = (TTree*)hfile3->Get(treeName);  
	//TTree *tree4 = (TTree*)hfile4->Get(treeName);  
	//TTree *tree5 = (TTree*)hfile5->Get(treeName);  

	TTree *tree_data = (TTree*)hfile_data->Get(treeName);  

#ifdef DEBUG
	std::cout<<"declared pointers to trees in input files"<<std::endl;
#endif

	///set histo names, number of bins, and axis limits here
	vector<TH1F*> h_Mll = MakeNHistos("h_Mll",2,20,50,250);
	vector<TH1F*> h_l1pt = MakeNHistos("h_l1pt",2,22,0,220);
	vector<TH1F*> h_l2pt = MakeNHistos("h_l2pt",2,13,0,130);
	vector<TH1F*> h_l1eta = MakeNHistos("h_l1eta", 2,20,-3.,3.);
	vector<TH1F*> h_l2eta = MakeNHistos("h_l2eta", 2,20,-3.,3.);
	vector<TH1F*> h_l1phi = MakeNHistos("h_l1phi", 2,20,-3.15,3.15);
	vector<TH1F*> h_l2phi = MakeNHistos("h_l2phi", 2,20,-3.15,3.15);
	vector<TH1F*> h_nleptons = MakeNHistos("h_nleptons", 2,4,0,4);
	vector<TH1F*> h_nvertices = MakeNHistos("h_nvertices", 2,33,0,33);
	vector<TH1F*> h_dR_l1l2 = MakeNHistos("h_dR_l1l2", 2,20,0,5);
	
#ifdef DEBUG
	std::cout<<"made vectors of TH1F pointers to histos"<<std::endl;
#endif

	int nhistos = 10;	//max is 10
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
	//std::vector<TH1F*> histos1(nhistos); // TTbar
	//std::vector<TH1F*> histos3(nhistos); // WZ
	//std::vector<TH1F*> histos4(nhistos); // ZZ
	//std::vector<TH1F*> histos5(nhistos); // WJets
	std::vector<TH1F*> histos_data(nhistos); // data

	///load the vectors tied to all MC and real data samples into one vector<vector<TH1F*> > object
	vector<vector<TH1F*> > histos(2);
	histos[0] = histos0;//DY
	//histos[1] = histos1;//TTBar
	//histos[2] = histos3;//WZ
	//histos[3] = histos4;//ZZ
	//histos[4] = histos5;//WJets
	histos[1] = histos_data;

	///link the histos made by calls to MakeNHistos() above to specific elements of vector<vector<TH1F*> > container 
	for(int j=0; j<2; j++){
		histos[j][0] = h_Mll[j];
		histos[j][1] = h_l1pt[j];
		histos[j][2] = h_l2pt[j];
		histos[j][3] = h_l1eta[j];
		histos[j][4] = h_l2eta[j];
		histos[j][5] = h_l1phi[j];
		histos[j][6] = h_l2phi[j];
		histos[j][7] = h_nleptons[j];
		histos[j][8] = h_nvertices[j];
		histos[j][9] = h_dR_l1l2[j];
	}
	
#ifdef DEBUG
	std::cout<<"linked MC processes and real data to histos of specific quantities"<<std::endl;
#endif

	///calculate PU weights and store them in persistent vectors of floats
	std::vector<float> PUW0 = PileUpWeights(tree0,tree_data);
	//std::vector<float> PUW1 = PileUpWeights(tree1,tree_data);
	//std::vector<float> PUW3 = PileUpWeights(tree3,tree_data);
	//std::vector<float> PUW4 = PileUpWeights(tree4,tree_data);
	//std::vector<float> PUW5 = PileUpWeights(tree5,tree_data);
	std::vector<float> PUW_data = PileUpWeights(tree_data,tree_data);

	///fill the histos with content from TTrees
	Fill_Histo(histos[0],tree0,PUW0,false,false); // DY
	//Fill_Histo(histos[1],tree1,PUW1,false,false); // TTbar
	//Fill_Histo(histos[2],tree3,PUW3,false,false); // WZ
	//Fill_Histo(histos[3],tree4,PUW4,false,false); // ZZ
	//Fill_Histo(histos[4],tree5,PUW5,false,false); // WJets

	Fill_Histo(histos[1],tree_data,PUW_data,false,false);	///real data

	//Float_t intLumi = 444;
	// Scale = xsection*luminosity/events
	for(std::vector<TH1F*>::size_type i = 0; i != nhistos; i++){
#ifdef DEBUG
		std::cout<<"setting line and fill colors, and adding histos to THStack objects"<<std::endl;
#endif
				
		histos[0][i]->Scale(1/histos[0][i]->Integral());	///aMC@NLO DYJets 50ns
		histos[0][i]->SetFillColor(5);
		
		histos[1][i]->Scale(1/histos[1][i]->Integral());	///madgraph DYJets 50ns
		histos[1][i]->SetMarkerStyle(20);

		ths[i]->Add(histos[0][i]);
		
		///rescale the max y value to make room for the legend
		//Double_t oldMax = ( (ths[i]->GetMaximum() > histos[5][i]->GetMaximum() ) ? ths[i]->GetMaximum() : histos[5][i]->GetMaximum() );
		//ths[i]->SetMaximum(20*oldMax);
		//histos[5][i]->SetMaximum(20*oldMax);
		//ths[i]->SetMinimum(0.2);
		//histos[5][i]->SetMinimum(0.2);
	}


	///make a legend with appropriate labels for MC processes
	TLegend *leg = new TLegend( 0.62, 0.66, 0.89, 0.89 ) ;
	leg->SetNColumns(2);
	leg->AddEntry( histos[1][0], "DYJets Madgraph" ,"p") ;
	leg->AddEntry( histos[0][0], "DYJets AMCNLO" ) ; 
	//leg->AddEntry( histos[1][0], "ttbar" ) ;
	//leg->AddEntry( histos[4][0], "WJets" ) ;  
	//leg->AddEntry( histos[2][0], "WZ" ) ; 
	//leg->AddEntry( histos[3][0], "ZZ" ) ;
	leg->SetFillColor( kWhite ) ;

	TString xtitles[] = {"M_{EE} [GeV]","leading electron p_{T} [GeV]","subleading electron p_{T} [GeV]","leading electron #eta","subleading electron #eta","leading electron #phi","subleading electron #phi","number of electrons","number of vertices","#DeltaR lead ele sublead ele"};
	
	TString titles[] = {"CMS Preliminary GEN DiElectron Mass  #surds = 13 TeV 50ns","CMS Preliminary GEN Lead Electron p_{T}  #surds = 13 TeV 50ns","CMS Preliminary GEN Sublead Electron p_{T}  #surds = 13 TeV 50ns","CMS Preliminary GEN Lead Electron #eta  #surds = 13 TeV 50ns","CMS Preliminary GEN Sublead Electron #eta  #surds = 13 TeV 50ns","CMS Preliminary GEN leading electron #phi  #surds = 13 TeV 50ns","CMS Preliminary GEN Subleading electron #phi  #surds = 13 TeV 50ns","CMS Preliminary GEN number of electrons  #surds = 13 TeV 50ns","CMS Preliminary GEN number of vertices  #surds = 13 TeV 50ns","CMS Preliminary GEN #DeltaR lead ele Sublead ele  #surds = 13 TeV 50ns"};

	TString fnames[] = {"GEN_MEE","GEN_l1_pt","GEN_l2_pt","GEN_l1_eta","GEN_l2_eta","GEN_l1_phi","GEN_l2_phi","GEN_nleptons","GEN_nvertices","GEN_dR_l1l2"};

#ifdef DEBUG
	std::cout<<"declared titles and file names for histos"<<std::endl;
#endif

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
		Double_t oldMax = ( (ths[icanvas]->GetMaximum() > histos[1][icanvas]->GetMaximum() ) ? ths[icanvas]->GetMaximum() : histos[1][icanvas]->GetMaximum() );
		ths[icanvas]->SetMaximum(1.2*oldMax);
		histos[1][icanvas]->SetMaximum(1.2*oldMax);
	
		ths[icanvas]->Draw("hist");
		histos[1][icanvas]->Draw("psame");
		mycanvas->Update();
	
		///set the x axis title
		ths[icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
		histos[1][icanvas]->GetXaxis()->SetTitle(xtitles[icanvas].Data());
	
		///set the y axis title using the bin size from the data histo, and the x axis title
		TString yLabel = DetermineYaxisName(histos[1][icanvas], xtitles[icanvas]);
		ths[icanvas]->GetYaxis()->SetTitle(yLabel );
		histos[1][icanvas]->GetYaxis()->SetTitle(yLabel );

		///set the histogram title to show a name
		ths[icanvas]->SetTitle(titles[icanvas]);
		histos[1][icanvas]->SetTitle(titles[icanvas]);

		///draw the legend, and a text box above the legend with the COM energy and integrated lumi
		leg->Draw();
		
		mycanvas->Update();
		TString tag = "_50ns";
		
		TString fn = "tempPlots/compareDYJets50ns/";
		TString fn_pdf = fn + fnames[icanvas].Data() + tag + "_checkZee.pdf";
		TString fn_png = fn + fnames[icanvas].Data() + tag + "_checkZee.png";
		mycanvas->Print(fn_pdf.Data());
		mycanvas->Print(fn_png.Data());
		mycanvas->SetLogy();

		///rescale vertical axis for log scale plot
		ths[icanvas]->SetMaximum(10*oldMax);
		histos[1][icanvas]->SetMaximum(10*oldMax);
		ths[icanvas]->SetMinimum(0.1);
		histos[1][icanvas]->SetMinimum(0.1);

		leg->Draw();
		mycanvas->Update();
		TString fn_log_pdf = fn + fnames[icanvas].Data() + tag + "_checkZee_log.pdf";
		TString fn_log_png = fn + fnames[icanvas].Data() + tag + "_checkZee_log.png";
		mycanvas->Print(fn_log_pdf.Data());
		mycanvas->Print(fn_log_png.Data());
		mycanvas->Close();

	}
}///end dyJetsMC_checkZee_compare()


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
	Float_t ptEle[2];
	Float_t etaEle[2];
	Float_t phiEle[2];
	Float_t dileptonMass;
	
	Float_t dR_l1l2;
	Int_t nLeptons;
	Int_t nVertices;
	Float_t evWeightSign;

	tree->SetBranchAddress("ptEle",ptEle);
	tree->SetBranchAddress("etaEle",etaEle);
	tree->SetBranchAddress("phiEle",phiEle);
	tree->SetBranchAddress("dileptonMass",&dileptonMass);
	
	tree->SetBranchAddress("dR_leadingLeptonSubleadingLepton",&dR_l1l2);
	tree->SetBranchAddress("nLeptons",&nLeptons);
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

				h1[0]->Fill(dileptonMass,reweight);
				h1[1]->Fill(ptEle[0],reweight);
				h1[2]->Fill(ptEle[1],reweight);
				h1[3]->Fill(etaEle[0],reweight);
				h1[4]->Fill(etaEle[1],reweight);
				h1[5]->Fill(phiEle[0],reweight);
				h1[6]->Fill(phiEle[1],reweight);
				h1[7]->Fill(nLeptons,reweight);	
				h1[8]->Fill(nVertices,reweight);	
				h1[9]->Fill(dR_l1l2,reweight);
			}///end if(!is_data)
			else {
				h1[0]->Fill(dileptonMass);
				h1[1]->Fill(ptEle[0]);
				h1[2]->Fill(ptEle[1]);
				h1[3]->Fill(etaEle[0]);
				h1[4]->Fill(etaEle[1]);
				h1[5]->Fill(phiEle[0]);
				h1[6]->Fill(phiEle[1]);
				h1[7]->Fill(nLeptons);	
				h1[8]->Fill(nVertices);	
				h1[9]->Fill(dR_l1l2);
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
