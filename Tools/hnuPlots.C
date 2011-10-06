#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"

#include "tdrstyle.C"

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>

const int NCOLORS = 5;
const int colors[] = {
	kCyan+2,
    kYellow,
    kGreen-2,
    kOrange+1,
    kMagenta-1
};

bool compHistInt(std::pair<string, TH1*> h1, std::pair<string, TH1*> h2)
{
   return h1.second->Integral() < h2.second->Integral();
}

void plot(std::vector<std::pair<std::string, TFile*> >& vpbg, std::pair<std::string, TFile*>& pdata, char * plot = "hNu/cut6_Mu1HighPt/mMuMu", char * xaxislabel = "M_{#mu#mu} (GeV)", char * yaxislabel = "counts", int rebin = -1)
{
    //gROOT->SetStyle("Plain");
    setTDRStyle();

	using namespace std;
	
	vector<pair<string, TH1*> > bghists;
	int iColor = 0;
	for(vector<pair<string, TFile*> >::const_iterator ibgf = vpbg.begin(); ibgf != vpbg.end(); ibgf++)
	{
		bghists.push_back(pair<string, TH1*>(ibgf->first, (TH1*)ibgf->second->Get(plot)->Clone()));
		bghists.back().second->SetFillColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetLineColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetMarkerColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetLineWidth(0);
		if(rebin > 1) bghists.back().second->Rebin(rebin);
		iColor++;
	}    
	sort(bghists.begin(), bghists.end(), compHistInt);
    
    pair<string, TH1*> hd(pdata.first, (TH1*)pdata.second->Get(plot)->Clone());
    if(rebin > 1) hd.second->Rebin(rebin);
   
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->cd();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);
    
    TLegend *leg = new TLegend(0.60, 0.70, 0.94, 0.94);
    
    THStack *hbg = new THStack("Background", "background");
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++) 
    {
    	hbg->Add(ihbg->second);
    }
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.end()-1; ihbg != bghists.begin()-1; ihbg--) 
    {
    	leg->AddEntry(ihbg->second, ihbg->first.c_str());
    }
    leg->AddEntry(hd.second, hd.first.c_str());
    
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, hd.second->GetBinLowEdge(1), hd.second->GetBinLowEdge(hd.second->GetNbinsX())+hd.second->GetBinWidth(hd.second->GetNbinsX()));
    dummy->GetXaxis()->SetTitle(xaxislabel);
    dummy->GetYaxis()->SetRangeUser(0.001, std::max(hbg->GetMaximum(),hd.second->GetMaximum())*1.2);
    dummy->GetYaxis()->SetTitle(yaxislabel);
    dummy->GetYaxis()->SetTitleOffset(1.0);
    dummy->SetStats(0);
    dummy->SetTitle(0);
    
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    
    //TLatex* mark2 = new TLatex(0.9, 0.5, "240 pb^{-1} at 7 TeV") ; 
    //mark2->SetTextSize(0.04);
    //mark2->SetTextFont(42);
    
    hd.second->SetMarkerColor(kBlack);
    hd.second->SetMarkerStyle(20);
    hd.second->SetLineWidth(2.0);
    
    fixOverlay();
    dummy->Draw();
    fixOverlay();
    hbg->Draw("hist same");
    fixOverlay();
    hd.second->Draw("same");
    fixOverlay();
    leg->Draw("same");
    mark.DrawLatex(0.65, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.62, 0.66, "1804 pb^{-1} at 7 TeV");
    //mark2->Draw("same");
}

void plot2011Wr()
{
	using namespace std;
	
	//background legend label, TFile
	vector<pair<string, TFile*> > bg;
	bg.push_back(pair<string, TFile*>("t#bar{t}",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_1804ipbKFACTOR.root")));
    bg.push_back(pair<string, TFile*>("Z+Jets", new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_1804ipbKFACTOR.root")));
    bg.push_back(pair<string, TFile*>("W+Jets", new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_1804ipb.root")));
    bg.push_back(pair<string, TFile*>("VV",     new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/vv.root")));
    bg.push_back(pair<string, TFile*>("tW",     new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/tW.root")));

    //data
    pair<string, TFile*> data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/1804ipb/run2011_mu24mu40_1804ipb.root"));
    
    plot(bg, data, "hNu/cut5_diLmass/mWR", "M_{W_{R}} (GeV)", "Events/40 GeV", -1);
}

void plot2011Arb(char * plotname, char * ylabel)
{
	using namespace std;
	
	//background (legend label, TFile)
	vector<pair<string, TFile*> > bg;
	bg.push_back(pair<string, TFile*>("#lower[-10.0]{t#bar{t}} background",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_1804ipb.root")));
    bg.push_back(pair<string, TFile*>("Z+Jets background", new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_1804ipb.root")));
    bg.push_back(pair<string, TFile*>("W+Jets background", new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_1804ipb.root")));
    bg.push_back(pair<string, TFile*>("VV background",     new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/vv.root")));
    bg.push_back(pair<string, TFile*>("tW background",     new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/test/rootfiles/tW.root")));

    //data
    pair<string, TFile*> data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/1804ipb/run2011_mu24mu40_1804ipb.root"));
    
    plot(bg, data, plotname, ylabel, "Events/40 GeV", -1);
}














