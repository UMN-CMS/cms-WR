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

class HnuPlots
{
public:
    struct FileStruct
    {
        std::string label;
        TFile* file;
        std::string histpath;
        
        FileStruct(std::string s, TFile* f, std::string h);
    };

    //HnuPlots();
    HnuPlots(FileStruct& fdata, std::vector<FileStruct>& vfbg, std::vector<FileStruct>& vfsig);
    void plot();
    void plotMCComp();
    void setRebin(int rbval);
    void setXAxisTitle(std::string label);
    void setYAxisTitle(std::string label);
    void setAutoSort(bool as);
    void setLog(bool log);
    
private:
    //bool compHistInt(std::pair<std::string, TH1*> h1, std::pair<std::string, TH1*> h2);
    std::vector<std::pair<std::string, TH1*> > bghists;
    std::vector<std::pair<std::string, TH1*> > sighists;
    std::pair<std::string, TH1*> datahist;

    int rebin;
    std::string xaxislabel;
    std::string yaxislabel;
    bool autosort, islog;
};

bool compHistInt(std::pair<std::string, TH1*> h1, std::pair<std::string, TH1*> h2)
{
   return h1.second->Integral() < h2.second->Integral();
}

HnuPlots::HnuPlots(HnuPlots::FileStruct& fdata, std::vector<HnuPlots::FileStruct>& vfbg, std::vector<HnuPlots::FileStruct>& vfsig)
{
    using namespace std;
    
    int iColor = 0;
	for(vector<HnuPlots::FileStruct>::const_iterator ibgf = vfbg.begin(); ibgf != vfbg.end(); ibgf++)
	{
		bghists.push_back(pair<string, TH1*>(ibgf->label, (TH1*)ibgf->file->Get(ibgf->histpath.c_str())->Clone()));
		bghists.back().second->SetFillColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetLineColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetMarkerColor(colors[iColor%NCOLORS]);
		bghists.back().second->SetLineWidth(0);
		iColor++;
	}

	iColor = 2;
	
   	for(vector<HnuPlots::FileStruct>::const_iterator isigf = vfsig.begin(); isigf != vfsig.end(); isigf++)
	{
		sighists.push_back(pair<string, TH1*>(isigf->label, (TH1*)isigf->file->Get(isigf->histpath.c_str())->Clone()));
		sighists.back().second->SetFillColor(kWhite+iColor);
		sighists.back().second->SetLineColor(kWhite+iColor);
		sighists.back().second->SetMarkerColor(kWhite+iColor);
		sighists.back().second->SetLineWidth(0);
		iColor++;
	}

    datahist = pair<string, TH1*>(fdata.label, (TH1*)fdata.file->Get(fdata.histpath.c_str())->Clone());
    
    rebin = -1;
    xaxislabel = "";
    yaxislabel = "";
    autosort = false;
    islog = false;
}

HnuPlots::FileStruct::FileStruct(std::string l, TFile* f, std::string h)
{
    label = l;
    file = f;
    histpath = h;
}

void HnuPlots::plot()
{
    using namespace std;
	
    //gROOT->SetStyle("Plain");
    setTDRStyle();

	if(rebin > 1) 
    {std::cout << "HERE" <<std::endl;
        datahist.second->Rebin(rebin);
        for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
        	ihbg->second->Rebin(rebin);
        }  
    }
    
    if(autosort) sort(bghists.begin(), bghists.end(), compHistInt);
   
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->cd();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);
    
    TLegend *leg = new TLegend(0.60, 0.70, 0.94, 0.94);
    
    THStack *hbg = new THStack("Background", "background");
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.end()-1; ihbg != bghists.begin()-1; ihbg--) 
    {
    	hbg->Add(ihbg->second);
    }
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++) 
    {
        float integral = ihbg->second->Integral(0,ihbg->second->GetNbinsX()+1);
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)",ihbg->first.c_str(), integral);
    	leg->AddEntry(ihbg->second, hllabel);
    }
    float integral = datahist.second->Integral(0,datahist.second->GetNbinsX()+1);
    char hllabel[128];
    sprintf(hllabel, "%s (%.0f)",datahist.first.c_str(), integral);
    leg->AddEntry(datahist.second, hllabel);
    
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.second->GetBinLowEdge(1), datahist.second->GetBinLowEdge(datahist.second->GetNbinsX())+datahist.second->GetBinWidth(datahist.second->GetNbinsX()));
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    if(islog)
    {
        dummy->GetYaxis()->SetRangeUser(0.03, std::max(hbg->GetMaximum(),datahist.second->GetMaximum())*8);
        c1->SetLogy(1);
    }
    else
    {
        dummy->GetYaxis()->SetRangeUser(0.001, std::max(hbg->GetMaximum(),datahist.second->GetMaximum())*1.2);
    }
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
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
    
    datahist.second->SetMarkerColor(kBlack);
    datahist.second->SetMarkerStyle(20);
    datahist.second->SetLineWidth(2.0);
    
    fixOverlay();
    dummy->Draw();
    fixOverlay();
    hbg->Draw("hist same");
    fixOverlay();
    datahist.second->Draw("same");
    fixOverlay();
    leg->Draw("same");
    mark.DrawLatex(0.65, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.62, 0.66, "2140 pb^{-1} at 7 TeV");
    //mark2->Draw("same");
}

void HnuPlots::plotMCComp()
{
    using namespace std;
	
    //gROOT->SetStyle("Plain");
    setTDRStyle();

	if(rebin > 1) 
    {
        for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
        	ihbg->second->Rebin(rebin);
        }  
    }
    
    //sort(bghists.begin(), bghists.end(), compHistInt);
   
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->cd();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);
    
    TLegend *leg = new TLegend(0.60, 0.70, 0.94, 0.94);
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.end()-1; ihbg != bghists.begin()-1; ihbg--) 
    {
    	leg->AddEntry(ihbg->second, ihbg->first.c_str());
    }
    //leg->AddEntry(datahist.second, datahist.first.c_str());
    
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.second->GetBinLowEdge(1), datahist.second->GetBinLowEdge(datahist.second->GetNbinsX())+datahist.second->GetBinWidth(datahist.second->GetNbinsX()));
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetYaxis()->SetRangeUser(0.001, std::max(bghists[0].second->GetMaximum(),bghists[1].second->GetMaximum())*1.2);
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.0);
    dummy->SetStats(0);
    dummy->SetTitle(0);
    
    fixOverlay();
    dummy->Draw();
    const int compcolors[] = {kRed, kBlue, kBlack, kGreen+2};
    int i = 0;
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++) 
    {
        ihbg->second->SetLineColor(compcolors[i%3]);
        ihbg->second->SetFillColor(0);
        ihbg->second->SetMarkerStyle(21 + i);
        ihbg->second->SetMarkerColor(compcolors[i%3]);
    	ihbg->second->Draw("same");
    	fixOverlay();
    	i++;
    }
    leg->Draw("same");
}

void HnuPlots::setRebin(int rbval)
{
    rebin = rbval;
}

void HnuPlots::setXAxisTitle(std::string label)
{
    xaxislabel = label;
}

void HnuPlots::setYAxisTitle(std::string label)
{
    yaxislabel = label;
}

void HnuPlots::setAutoSort(bool as)
{
    autosort = as;
}

void HnuPlots::setLog(bool log)
{
    islog = log;
}
	
const std::string cutlevels[] =
{
    "cut0_none",
    "cut1_LLJJpt",
    "cut2_TrigMatches",
    "cut3_Vertex",
    "cut4_Mu1HighPt",
    "cut5_diLmass",
    "cut6_mWRmass"
};

void plot2011Wr(int cutlevel = 5)
{
	using namespace std;

	//background legend label, TFile
	vector<HnuPlots::FileStruct> bg, sig;
	//HnuPlots::FileStruct
	//bg.push_back(HnuPlots::FileStruct("t#bar{t}",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2140ipbKFACTOR.root"),"hNu/" + cutlevels[cutlevel] + "/mWR"));
    //bg.push_back(HnuPlots::FileStruct("Z+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2140ipbKFACTOR.root"),        "hNu/" + cutlevels[cutlevel] + "/mWR"));
    bg.push_back(HnuPlots::FileStruct("W+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2140ipb.root"),                    "hNu/" + cutlevels[cutlevel] + "/mWR"));
    bg.push_back(HnuPlots::FileStruct("tW",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/tW_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/mWR"));
    bg.push_back(HnuPlots::FileStruct("VV",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/vv_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/mWR"));
    bg.push_back(HnuPlots::FileStruct("Z+Jets Madgraph",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuAnalysis_DYJetsToLL_madgraph_2140ipb.root"),                  "hNuMu40/" + cutlevels[cutlevel] + "/mWR"));
    bg.push_back(HnuPlots::FileStruct("t#bar{t} Madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/ttbar-PFJets_2140ipb.root"),                                         "hNu/" + cutlevels[cutlevel] + "/mWR"));
    
    //data
    HnuPlots::FileStruct data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/2140ipb/RecoEff/oct18/run2011A-recoEff-2140ipb.root"), "hNu/" + cutlevels[cutlevel] + "/mWR");
    
    HnuPlots hps(data, bg, sig);
    hps.setXAxisTitle("M_{W_{R}} [GeV]");
    hps.setYAxisTitle("Events");
    hps.setRebin(2);
    hps.plot();
}

void plotMCVar(int cutlevel, std::string name, std::string xaxis = "M_{W_{R}} [GeV]")
{
    using namespace std;
	
	//background legend label, TFile
	vector<HnuPlots::FileStruct> bg, sig;
	//HnuPlots::FileStruct
	//bg.push_back(HnuPlots::FileStruct("DY madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuAnalysis_DYJetsToLL_madgraph_2140ipb.root"),                    "hNuMu40/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("DY sherpa",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2140ipbKFACTOR.root"), "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("t#bar{t} madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuTopAnalysis_TTJets_TuneZ2_7TeV-madgraph-tauola_oct21_2140ipb.root"),                    "hNuTopHiLumi/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("t#bar{t} powheg",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_Top_2140ipb.root"),                    "hNuTopHiLumi/" + cutlevels[cutlevel] + "/" + name));
    bg.push_back(HnuPlots::FileStruct("t#bar{t} madgraph",  new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6.root"),                    "hNu/" + cutlevels[cutlevel] + "/" + name));
    bg.push_back(HnuPlots::FileStruct("t#bar{t} powheg",    new TFile("/local/cms/user/dahmes/wr2011/bgMC/Summer11/aug30/ttbar-PFJets.root"),                    "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("e#mu Data",    new TFile("/local/cms/user/dahmes/wr2011/data_2140ipb/run2011A-top.root"),                    "hNuTop/" + cutlevels[cutlevel] + "/" + name));        
    //bg.push_back(HnuPlots::FileStruct("W+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2140ipb.root"),                    "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("VV",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/vv_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("tW",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/tW_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/" + name));

    //data
    HnuPlots::FileStruct data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/2140ipb/RecoEff/oct18/run2011A-recoEff-2140ipb.root"), "hNu/" + cutlevels[cutlevel] + "/" + name);
    
    HnuPlots hps(data, bg, sig);
    hps.setXAxisTitle(xaxis.c_str());
    hps.setYAxisTitle("Events");
    hps.setRebin(2);
    hps.plotMCComp();
}

void plot2011(int cutlevel = 5, std::string plot = "mWR", bool log = true)
{
	using namespace std;

	//background legend label, TFile
	vector<HnuPlots::FileStruct> bg, sig;
	//bg.push_back(HnuPlots::FileStruct("t#bar{t}",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2140ipb.root"),"hNu/" + cutlevels[cutlevel] + "/" + plot));
    bg.push_back(HnuPlots::FileStruct("Z+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2140ipbKFACTOR.root"),        "hNu/" + cutlevels[cutlevel] + "/" + plot));
    //bg.push_back(HnuPlots::FileStruct("t#bar{t} Madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/ttbar-PFJets_2140ipb.root"),                                         "hNu/" + cutlevels[cutlevel] + "/" + plot));
    //bg.push_back(HnuPlots::FileStruct("Z+Jets Madgraph",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuAnalysis_DYJetsToLL_madgraph_2140ipb.root"),                  "hNuMu40/" + cutlevels[cutlevel] + "/" + plot));
    //bg.push_back(HnuPlots::FileStruct("tW",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/tW_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/" + plot));
    //bg.push_back(HnuPlots::FileStruct("VV",        new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/vv_2140ipb.root"),                                                            "hNu/" + cutlevels[cutlevel] + "/" + plot));
    //bg.push_back(HnuPlots::FileStruct("W+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2140ipb.root"),                    "hNu/" + cutlevels[cutlevel] + "/" + plot));
    bg.push_back(HnuPlots::FileStruct("Background",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/Zotherbg_2140ipb.root"),"hNu/" + cutlevels[cutlevel] + "/" + plot));
	
    //data
    HnuPlots::FileStruct data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/2140ipb/RecoEff/oct20/run2011A-recoEff-2140ipb.root"), "hNu/" + cutlevels[cutlevel] + "/" + plot);
    
    HnuPlots hps(data, bg, sig);
    hps.setXAxisTitle("M_{#mu#mu} [GeV]");
    hps.setYAxisTitle("Events");
    hps.setRebin(3);
    hps.setLog(log);
    hps.plot();
}














