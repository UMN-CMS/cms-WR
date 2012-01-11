#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"

#include "tdrstyle.C"

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cmath>


const int colors[] = {
    kCyan + 2,
    kYellow,
    kGreen - 2,
    kOrange + 1,
    kMagenta - 1,
    kRed,
    kBlue,
    kGreen
};
const int NCOLORS = sizeof(colors) / sizeof(int);

int terribleroothackcounter = 0;

class HnuPlots
{
public:

    struct FileStruct
    {
        std::string label;
        TFile* file;
        std::string histpath;
        double intLumi, cs, kfactor;
        std::string normhistpath;
        double clow, chigh;
        bool projX;

        FileStruct(){ }
        FileStruct(std::string l, TFile* f, std::string h);
        FileStruct(std::string l, TFile* f, std::string h, double iL, double c, double kf, std::string nh, double cl = 0.0, double ch = 0.0, bool px = true);
    } ;

    //HnuPlots();
    HnuPlots(FileStruct& fdata, std::vector<FileStruct>& vfbg, std::vector<FileStruct>& vfsig);
    HnuPlots(FileStruct& fdata, std::vector<std::vector<FileStruct> >& vfbg, std::vector<std::vector<FileStruct> >& vfsig, double iL);
    void plot();
    void plotMCComp();
    void plotMCShape(std::string bgfilename);
    void plotQCDFit();
    void plotNorm(double lower, double upper);
    void cutFlow();
    void integrals(double min, double max);
    void setRebin(int rbval);
    void setXAxisTitle(std::string label);
    void setYAxisTitle(std::string label);
    void setAutoSort(bool as);
    void setLog(bool log);
    void setFormLabel(std::string);
    void setXRange(double min, double max);
    void setSavePlots(bool sp);

private:
    std::vector<std::pair<std::string, TH1*> > bghists;
    std::vector<std::pair<std::string, TH1*> > sighists;
    std::pair<std::string, TH1*> datahist;

    int rebin;
    double iLumi, xmin, xmax;
    std::string xaxislabel;
    std::string yaxislabel;
    std::string formlabel;
    bool autosort, islog, saveplots;

    std::string getHistogramXAxisTitle(std::string name);

    TH1* project(TH2* h2d, double cl, double ch, bool porjx = true);
    int projcount;
} ;

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
        TH1 * h = (TH1*)ibgf->file->Get(ibgf->histpath.c_str());
        if(h)
        {
            bghists.push_back(pair<string, TH1*>(ibgf->label, (TH1*)h->Clone()));
            bghists.back().second->SetFillColor(colors[iColor % NCOLORS]);
            bghists.back().second->SetLineColor(colors[iColor % NCOLORS]);
            bghists.back().second->SetMarkerColor(colors[iColor % NCOLORS]);
            bghists.back().second->SetLineWidth(0);
            iColor++;
        }
        else std::cout << "failed to get File:hist - " << ibgf->file->GetName() << " : " << ibgf->histpath << std::endl;
    }

    iColor = 2;

    for(vector<HnuPlots::FileStruct>::const_iterator isigf = vfsig.begin(); isigf != vfsig.end(); isigf++)
    {
        TH1 * h = (TH1*)isigf->file->Get(isigf->histpath.c_str());
        if(h)
        {
            sighists.push_back(pair<string, TH1*>(isigf->label, (TH1*)h->Clone()));
            sighists.back().second->SetFillColor(kWhite + iColor);
            sighists.back().second->SetLineColor(kWhite + iColor);
            sighists.back().second->SetMarkerColor(kWhite + iColor);
            sighists.back().second->SetLineWidth(0);
            iColor++;
        }
        else std::cout << "failed to get File:hist - " << isigf->file->GetName() << " : " << isigf->histpath << std::endl;
    }

    TH1 *h = (TH1*)fdata.file->Get(fdata.histpath.c_str());
    if(h) datahist = pair<string, TH1*>(fdata.label, (TH1*)h->Clone());
    else std::cout << "failed to get File:hist - " << fdata.file->GetName() << " : " << fdata.histpath << std::endl;

    rebin = -1;
    xaxislabel = "";
    yaxislabel = "";
    formlabel = "hNu";
    autosort = false;
    islog = false;
    iLumi = 1.0;
    xmin = xmax = 0.0;
    saveplots = true;
    projcount = 0;
}

HnuPlots::HnuPlots(FileStruct& fdata, std::vector<std::vector<HnuPlots::FileStruct> >& vfbg, std::vector<std::vector<HnuPlots::FileStruct> >& vfsig, double iL)
{
    using namespace std;

    bool first;

    int iColor = 0;
    for(vector<vector<HnuPlots::FileStruct> >::const_iterator ivbg = vfbg.begin(); ivbg != vfbg.end(); ivbg++)
    {
        first = true;
        for(vector<HnuPlots::FileStruct>::const_iterator ibgf = ivbg->begin(); ibgf != ivbg->end(); ibgf++)
        {
            TH1 * h = NULL;
            if(fabs(ibgf->clow) < 1e-300 && fabs(ibgf->chigh) < 1e-300) h = (TH1*)ibgf->file->Get(ibgf->histpath.c_str());
            else
            {
                TH2 *h2 = (TH2*)ibgf->file->Get(ibgf->histpath.c_str());
                if(!h2) std::cout << "failed to get File:hist - " << ibgf->file->GetName() << " : " << ibgf->histpath << std::endl;
                h = project(h2, ibgf->clow, ibgf->chigh);
            }
            TH1 * hn = (TH1*)ibgf->file->Get(ibgf->normhistpath.c_str());
            if(h && hn)
            {
                h->Scale(ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->GetBinContent(1));
                //cout << ibgf->label << " : " << ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->GetBinContent(1) << endl;
                if(first)
                {
                    bghists.push_back(pair<string, TH1*>(ibgf->label, (TH1*)h->Clone()));
                    first = false;
                }
                else
                {
                    bghists.back().second->Add(h);
                }
            }
            else
            {
                if(!h) std::cout << "failed to get File:hist - " << ibgf->file->GetName() << " : " << ibgf->histpath << std::endl;
                if(!hn) std::cout << "failed to get File:hist - " << ibgf->file->GetName() << " : " << ibgf->normhistpath << std::endl;
            }
        }
        bghists.back().second->SetFillColor(colors[iColor % NCOLORS]);
        bghists.back().second->SetLineColor(colors[iColor % NCOLORS]);
        bghists.back().second->SetMarkerColor(colors[iColor % NCOLORS]);
        bghists.back().second->SetLineWidth(0);
        iColor++;
    }

    iColor = 2;
    for(vector<vector<HnuPlots::FileStruct> >::const_iterator ivsig = vfsig.begin(); ivsig != vfsig.end(); ivsig++)
    {
        first = true;
        for(vector<HnuPlots::FileStruct>::const_iterator isigf = ivsig->begin(); isigf != ivsig->end(); isigf++)
        {
            TH1 * h = (TH1*)isigf->file->Get(isigf->histpath.c_str());
            TH1 * hn = (TH1*)isigf->file->Get(isigf->normhistpath.c_str());
            if(h && hn)
            {
                h->Scale(isigf->intLumi * isigf->cs * isigf->kfactor / hn->GetBinContent(1));
                if(first)
                {
                    sighists.push_back(pair<string, TH1*>(isigf->label, (TH1*)h->Clone()));
                    first = false;
                }
                else
                {
                    sighists.back().second->Add(h);
                }
            }
            else
            {
                if(!h) std::cout << "failed to get File:hist - " << isigf->file->GetName() << " : " << isigf->histpath << std::endl;
                if(!hn) std::cout << "failed to get File:hist - " << isigf->file->GetName() << " : " << isigf->normhistpath << std::endl;
            }
        }
        //sighists.back().second->SetFillColor(kWhite + iColor);
        //sighists.back().second->SetFillStyle(0);
        sighists.back().second->SetLineColor(kWhite + iColor);
        sighists.back().second->SetMarkerColor(kWhite + iColor);
        sighists.back().second->SetLineWidth(0);
        iColor += 2;
    }

    TH1 *h;
    if(fabs(fdata.clow) < 1e-300 && fabs(fdata.chigh) < 1e-300) h = (TH1*)fdata.file->Get(fdata.histpath.c_str());
    else
    {
        TH2 *h2 = (TH2*)fdata.file->Get(fdata.histpath.c_str());
        if(!h2) std::cout << "failed to get File:hist - " << fdata.file->GetName() << " : " << fdata.histpath << std::endl;
        h = project(h2, fdata.clow, fdata.chigh);
    }
    if(h) datahist = pair<string, TH1*>(fdata.label, (TH1*)h->Clone());
    else std::cout << "failed to get File:hist - " << fdata.file->GetName() << " : " << fdata.histpath << std::endl;

    rebin = -1;
    xaxislabel = "";
    yaxislabel = "";
    formlabel = "hNu";
    autosort = false;
    islog = false;
    iLumi = iL;
    xmin = xmax = 0.0;
    saveplots = true;
    projcount = 0;
}

TH1* HnuPlots::project(TH2* h2d, double cl, double ch, bool porjx)
{
    char tmp[128];
    TH1* h = NULL;

    if(porjx)
    {
        sprintf(tmp, "%s_projx_%d", h2d->GetName(), projcount++);
        h = new TH1D(tmp, tmp, h2d->GetXaxis()->GetNbins() + 1, h2d->GetXaxis()->GetBinLowEdge(1), h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetNbins()));
        for(int i = h2d->GetYaxis()->FindBin(cl); i <= h2d->GetYaxis()->FindBin(ch); i++)
        {
            for(int j = 0; j < h2d->GetXaxis()->GetNbins() + 1; j++)
            {
                h->SetBinContent(j, h->GetBinContent(j) + h2d->GetBinContent(j, i));
                double be1 = h->GetBinError(j), be2 = h2d->GetBinError(j, i);
                h->SetBinError(j, sqrt(be1*be1+be2*be2));
            }
        }
    }
    else
    {
        sprintf(tmp, "%s_projy_%d", h2d->GetName(), projcount++);
        h = new TH1D(tmp, tmp, h2d->GetYaxis()->GetNbins() + 1, h2d->GetYaxis()->GetBinLowEdge(1), h2d->GetYaxis()->GetBinUpEdge(h2d->GetYaxis()->GetNbins()));
        for(int i = h2d->GetXaxis()->FindBin(cl); i <= h2d->GetXaxis()->FindBin(ch); i++)
        {
            for(int j = 0; j < h2d->GetYaxis()->GetNbins() + 1; j++)
            {
                h->SetBinContent(j, h->GetBinContent(j) + h2d->GetBinContent(i, j));
                double be1 = h->GetBinError(j), be2 = h2d->GetBinError(i, j);
                h->SetBinError(j, sqrt(be1*be1+be2*be2));
            }
        }
    }

    return h;
}

HnuPlots::FileStruct::FileStruct(std::string l, TFile* f, std::string h)
{
    label = l;
    file = f;
    histpath = h;
    intLumi = cs = kfactor = 1.0;
    clow = chigh = 0.0;
}

HnuPlots::FileStruct::FileStruct(std::string l, TFile* f, std::string h, double iL, double c, double kf, std::string nh, double cl, double ch, bool px)
{
    label = l;
    file = f;
    histpath = h;
    intLumi = iL;
    cs = c;
    kfactor = kf;
    normhistpath = nh;
    clow = cl;
    chigh = ch;
    projX = px;
}

void HnuPlots::plot()
{
    using namespace std;

    //gROOT->SetStyle("Plain");
    setTDRStyle();

    if(rebin > 1)
    {
        datahist.second->Rebin(rebin);
        for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->second->Rebin(rebin);
        }
    }

    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 7 TeV", iLumi / 1000);

    if(autosort) sort(bghists.begin(), bghists.end(), compHistInt);

    if(!xaxislabel.compare("please auto set the axis")) xaxislabel = getHistogramXAxisTitle(datahist.second->GetName());
    if(!yaxislabel.compare("please auto set the axis"))
    {
        char temp[128];
        sprintf(temp, "Events / %.0f GeV", datahist.second->GetBinWidth(1));
        yaxislabel = temp;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 900);
    c1->Divide(1, 2);
    c1->cd(1);
    gPad->SetPad("p1", "p1", 0, 2.0 / 9.0, 1, 1, kWhite, 0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.01);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);

    TLegend *leg = new TLegend(0.55, 0.70, 0.94, 0.94);

    float dataintegral = datahist.second->Integral(0, datahist.second->GetNbinsX() + 1);
    char datahllabel[128];
    sprintf(datahllabel, "%s (%.0f)", datahist.first.c_str(), dataintegral);
    leg->AddEntry(datahist.second, datahllabel);

    THStack *hbg = new THStack("Background", "background");
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.end() - 1; ihbg != bghists.begin() - 1; ihbg--)
    {
        hbg->Add(ihbg->second);
    }
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        float integral = ihbg->second->Integral(0, ihbg->second->GetNbinsX() + 1);
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)", ihbg->first.c_str(), integral);
        leg->AddEntry(ihbg->second, hllabel);
    }
    for(vector<pair<string, TH1*> >::const_iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
    {
        float integral = ihsig->second->Integral(0, ihsig->second->GetNbinsX() + 1);
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)", ihsig->first.c_str(), integral);
        leg->AddEntry(ihsig->second, hllabel);
    }

    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.second->GetBinLowEdge(1), datahist.second->GetBinLowEdge(datahist.second->GetNbinsX()) + datahist.second->GetBinWidth(datahist.second->GetNbinsX()));
    if(xmin != xmax) dummy->GetXaxis()->SetRangeUser(xmin, xmax);
    //dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    if(islog)
    {
        dummy->GetYaxis()->SetRangeUser(std::max(0.001, 0.2 * std::min(hbg->GetMaximum(), datahist.second->GetMinimum(0.0001))), std::max(hbg->GetMaximum(), datahist.second->GetMaximum())*8);
        gPad->SetLogy(1);
    }
    else
    {
        dummy->GetYaxis()->SetRangeUser(0.001, std::max(hbg->GetMaximum(), datahist.second->GetMaximum())*1.2);
    }
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.0);
    dummy->SetStats(0);
    dummy->SetTitle(0);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    datahist.second->SetMarkerColor(kBlack);
    datahist.second->SetMarkerStyle(20);
    datahist.second->SetLineWidth(2.0);

    fixOverlay();
    dummy->Draw();
    fixOverlay();
    hbg->Draw("hist same");
    fixOverlay();
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator isig = sighists.begin(); isig != sighists.end(); isig++)
    {
        isig->second->Draw("hist same");
    }
    datahist.second->Draw("same");
    fixOverlay();
    leg->Draw("same");
    mark.DrawLatex(0.65, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.65, 0.66, lumistamp);

    c1->cd(2);
    gPad->SetPad("p2", "p2", 0, 0, 1, 2.0 / 9.0, kWhite, 0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.50);
    TH1* chbg = 0;
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        if(ihbg == bghists.begin()) chbg = (TH1*)ihbg->second->Clone();
        else chbg->Add(ihbg->second);
    }
    TH1* chdata = (TH1*)datahist.second->Clone();
    double ratio = 0, error = 0;
    for(int i = 1; i <= chdata->GetNbinsX(); i++)
    {
        if(chdata->GetBinContent(i) != 0 && chbg->GetBinContent(i) != 0)
        {
            ratio = chdata->GetBinContent(i) / chbg->GetBinContent(i);
            error = ratio * sqrt(pow(chdata->GetBinError(i) / chdata->GetBinContent(i), 2) + pow(chbg->GetBinError(i) / chbg->GetBinContent(i), 2));
            chdata->SetBinContent(i, ratio);
            chdata->SetBinError(i, error);
        }
        else
        {
            chdata->SetBinContent(i, -1);
            chdata->SetBinError(i, 0);
        }
    }

    TH1 *dummy2 = new TH1F("dummy2", "dummy2", 1000, datahist.second->GetBinLowEdge(1), datahist.second->GetBinLowEdge(datahist.second->GetNbinsX()) + datahist.second->GetBinWidth(datahist.second->GetNbinsX()));
    dummy2->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy2->GetXaxis()->SetTitleOffset(1.05);
    dummy2->GetXaxis()->SetTitleSize(0.22);
    dummy2->GetXaxis()->SetLabelSize(0.22);
    dummy2->GetYaxis()->SetRangeUser(0, max(2.0, min(3.5, chdata->GetMaximum()*1.2)));
    dummy2->GetYaxis()->SetNdivisions(3, 5, 0);
    dummy2->GetYaxis()->SetTitle("Data/MC");
    dummy2->GetYaxis()->SetTitleOffset(0.25);
    dummy2->GetYaxis()->SetTitleSize(0.16);
    dummy2->GetYaxis()->SetLabelSize(0.22);
    if(xmin != xmax) dummy2->GetXaxis()->SetRangeUser(xmin, xmax);
    dummy2->SetStats(0);
    dummy2->SetTitle(0);
    if(dummy2->GetNdivisions() % 100 > 5) dummy2->GetXaxis()->SetNdivisions(6, 5, 0);

    TF1 *fline = new TF1("line", "pol0", datahist.second->GetBinLowEdge(1), datahist.second->GetBinLowEdge(datahist.second->GetNbinsX()) + datahist.second->GetBinWidth(datahist.second->GetNbinsX()));
    fline->SetParameter(0, 1);
    fline->SetLineColor(kRed);

    dummy2->Draw();
    fline->Draw("same");
    chdata->Draw("same");

    if(saveplots)
    {
        char ofn[128], tmp[32];
        int cutlevel = 11111;
        sscanf(strstr(datahist.second->GetTitle(), "cut:"), "cut:%d", &cutlevel);
        sprintf(tmp, "cut:%da", cutlevel);
        if(strstr(datahist.second->GetTitle(), tmp) != NULL) sprintf(ofn, "%s_cut%da_%s.pdf", formlabel.c_str(), cutlevel, datahist.second->GetName());
        else sprintf(ofn, "%s_cut%d_%s.pdf", formlabel.c_str(), cutlevel, datahist.second->GetName());
        c1->Print(ofn);
    }
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

    if(true)
    {
        for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin() + 1; ihbg != bghists.end(); ihbg++)
        {
            ihbg->second->Scale(bghists.begin()->second->Integral() / ihbg->second->Integral());
        }
    }

    //sort(bghists.begin(), bghists.end(), compHistInt);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);

    TLegend *leg = new TLegend(0.60, 0.70, 0.94, 0.94);
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.end() - 1; ihbg != bghists.begin() - 1; ihbg--)
    {
        leg->AddEntry(ihbg->second, ihbg->first.c_str());
    }
    //leg->AddEntry(datahist.second, datahist.first.c_str());

    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, bghists[0].second->GetBinLowEdge(1), bghists[0].second->GetBinLowEdge(bghists[0].second->GetNbinsX()) + bghists[0].second->GetBinWidth(bghists[0].second->GetNbinsX()));
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetYaxis()->SetRangeUser(0.001, std::max(bghists[0].second->GetMaximum(), bghists[1].second->GetMaximum())*1.2);
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.0);
    dummy->SetStats(0);
    dummy->SetTitle(0);

    fixOverlay();
    dummy->Draw();
    const int compcolors[] = {kRed, kBlue, kBlack, kGreen + 2};
    int i = 0;
    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        ihbg->second->SetLineColor(compcolors[i % 3]);
        ihbg->second->SetFillColor(0);
        ihbg->second->SetMarkerStyle(21 + i);
        ihbg->second->SetMarkerColor(compcolors[i % 3]);
        ihbg->second->Draw("same");
        fixOverlay();
        i++;
    }
    leg->Draw("same");
}

void HnuPlots::plotMCShape(std::string bgfilename)
{
    using namespace std;

    setTDRStyle();

    bool drawtrialfuncs = false;

    FILE *fbackfits = fopen(bgfilename.c_str(), "w");

    vector<pair<string, string> > fffs;
    fffs.push_back(make_pair("#it{e^{a+bM}}", "exp([0]+[1]*x)"));
    fffs.push_back(make_pair("#it{e^{a+bMlog(M)}}", "exp([0]+[1]*x*log(x))"));
    fffs.push_back(make_pair("#it{e^{a+bM+cM^{2}}}", "exp([0]+[1]*x+[2]*x*x)"));
    fffs.push_back(make_pair("#it{e^{a+bM+cM^{3}}}", "exp([0]+[1]*x+[2]*x*x*x)"));
    //fffs.push_back(make_pair("#it{e^{a+bM}/M^{c}}", "exp([0]+[1]*x)/pow(x,[2])"));
    fffs.push_back(make_pair("#it{e^{a+bM}+c}", "exp([0]+[1]*x)+[2]"));
    //fffs.push_back(make_pair("#it{e^{a+bM}+cM^{2}}", "exp([0]+[1]*x)+[2]*x*x"));
    //fffs.push_back(make_pair("#it{e^{a+bM^{C}}}", "exp([0]+[1]*pow(x,[2]))"));

    TCanvas * cans[bghists.size()];
    int fcount = 0;

    for(vector<pair<string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        cout << "Model: " << ihbg->first << endl;
        TF1 **f = new TF1*[fffs.size()], *f1, *f2, *f3;
        TLegend *leg = new TLegend(0.62, 0.64, 0.94, 0.94);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);

        char tmp[32];
        sprintf(tmp, "c%d", fcount);
        cans[fcount] = new TCanvas(tmp, tmp, 800, 800);
        cans[fcount]->cd();
        cans[fcount]->SetLeftMargin(0.15);
        cans[fcount]->SetRightMargin(0.06);
        cans[fcount]->SetTopMargin(0.06);
        cans[fcount++]->SetLogy(islog);

        TH1* h = (TH1*)ihbg->second->Clone();
        h->SetLineColor(kBlack);
        h->SetMarkerStyle(21);
        h->SetMarkerColor(kBlack);
        h->SetFillColor(kWhite);
        h->Rebin(rebin);
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitle();
        h->GetYaxis()->SetTitle();

        leg->AddEntry(h, ihbg->first.c_str());

        for(unsigned int ifunc = 0; ifunc < fffs.size(); ifunc++)
        {
            char tmp2[32];
            sprintf(tmp2, "f%d", ifunc);
            f[ifunc] = new TF1(tmp, fffs[ifunc].second.c_str(), 600, 2800);
            if(f[ifunc]->GetNpar() > 2 && strstr(fffs[ifunc].first.c_str(), "#it{e^{a+bM}+c}") != NULL)
            {
                f[ifunc]->SetParLimits(2, 0.0, 1000);
                f[ifunc]->SetParameter(2, 2);
            }
            h->Fit(f[ifunc], "LMNQR");
            cout << fffs[ifunc].second << " Chi^2/NDOF:" << f[ifunc]->GetChisquare() << " / " << f[ifunc]->GetNDF() << endl;
            cout << "par0: " << f[ifunc]->GetParameter(0) << "\t par1: " << f[ifunc]->GetParameter(1) << "\tpar1error: " << f[ifunc]->GetParError(1) << "\t par2: " << f[ifunc]->GetParameter(2) << endl;
            f[ifunc]->SetLineColor(colors[ifunc % NCOLORS]);
            f[ifunc]->SetLineWidth(2);
            f[ifunc]->SetMarkerStyle(0);
        }

        f1 = new TF1("f1", "exp([0]+[1]*x)", 600, 2800);
        f2 = new TF1("f2", "exp([0]+[1]*x)", 600, 2800);
        f3 = new TF1("f3", "exp([0]+[1]*x)", 600, 2800);
        h->Fit(f1, "LMNQ", "", 600, 2800);
        h->Fit(f2, "LMNQ", "", 600, 1400);
        h->Fit(f3, "LMNQ", "", 1400, 2800);
        double s1 = f1->GetParameter(1), s2 = f2->GetParameter(1), s3 = f3->GetParameter(1);
        //double e1 = f1->GetParError(1), e2 = f2->GetParError(1), e3 = f3->GetParError(1);
        //double diff1 = fabs(s1 - s2), diff2 = fabs(s2 - s3), diff1e = e1 * e1 + e2*e2, diff2e = e3 * e3 + e2*e2;
        //double slopeError = fabs((diff1 / diff1e + diff2 / diff2e) / (1 / diff1e + 1 / diff2e));

        TF1 *fu = new TF1("fu", "[0]*exp([1]*(x-[2]))", 600, 2800);
        fu->SetParameters(f1->Eval(600), max(s1, max(s2, s3)), 600);
        TF1 *fd = new TF1("fd", "[0]*exp([1]*(x-[2]))", 600, 2800);
        fd->SetParameters(f1->Eval(600), min(s1, min(s2, s3)), 600);

        TGraphAsymmErrors *tgshape = new TGraphAsymmErrors(), *tgslope = new TGraphAsymmErrors(), *tgtotal = new TGraphAsymmErrors();
        tgshape->SetFillColor(kYellow);
        tgslope->SetFillColor(kBlue + 1);
        tgtotal->SetFillColor(kGreen - 1);
        tgshape->SetLineColor(kYellow);
        tgslope->SetLineColor(kBlue + 1);
        tgtotal->SetLineColor(kGreen - 1);
        tgshape->SetMarkerStyle(0);
        tgslope->SetMarkerStyle(0);
        tgtotal->SetMarkerStyle(0);
        leg->AddEntry(tgshape, "shape systematic");
        leg->AddEntry(tgslope, "slope systematic");
        leg->AddEntry(tgtotal, "total systematic");

        int istart = h->FindBin(600);
        double integral = 0;
        if(!ihbg->first.compare("Other")) fprintf(fbackfits, "%s\t%s\t", "Other", "bgest");
        else if(!ihbg->first.compare("t#bar{t}")) fprintf(fbackfits, "%s\t%s\t", "TT", "bgest");
        else if(!ihbg->first.compare("Z+Jets")) fprintf(fbackfits, "%s\t%s\t", "ZJ", "bgest");
        else fprintf(fbackfits, "%s\t%s\t", ihbg->first.c_str(), "bgest");
        for(int i = istart; i <= h->GetNbinsX(); i++)
        {
            double max = 0.0, min = 1e300, x = h->GetBinCenter(i), y, ynom = f[0]->Eval(x);
            tgshape->SetPoint(i - istart, x, ynom);
            tgslope->SetPoint(i - istart, x, ynom);
            tgtotal->SetPoint(i - istart, x, ynom);
            for(unsigned int ifunc = 0; ifunc < fffs.size(); ifunc++)
            {
                y = f[ifunc]->Eval(x);
                max = std::max(max, y);
                min = std::min(min, y);
            }
            double edo = 0, eup = 0;
            if(!ihbg->first.compare("Other"))
            {
                tgshape->SetPointError(i - istart, 0.0, 0.0, 0.0, 0.0);
                tgslope->SetPointError(i - istart, 0.0, 0.0, ynom - fd->Eval(x), fu->Eval(x) - ynom);
                edo = ynom - fd->Eval(x);
                eup = fu->Eval(x) - ynom;
                tgtotal->SetPointError(i - istart, 0.0, 0.0, edo, eup);
            }
            else
            {
                tgshape->SetPointError(i - istart, 0.0, 0.0, ynom - min, max - ynom);
                tgslope->SetPointError(i - istart, 0.0, 0.0, ynom - fd->Eval(x), fu->Eval(x) - ynom);
                edo = sqrt(pow(ynom - fd->Eval(x), 2) + pow(ynom - min, 2));
                eup = sqrt(pow(fu->Eval(x) - ynom, 2) + pow(max - ynom, 2));
                tgtotal->SetPointError(i - istart, 0.0, 0.0, edo, eup);
            }
            if(fmod(h->GetBinLowEdge(i), 200) == 0.0)
            {
                if(h->GetBinLowEdge(i) > 600) fprintf(fbackfits, "%f\t", integral);
                integral = 0.0;
            }
            integral += (ynom + eup + ynom - edo) / 2;
        }
        fprintf(fbackfits, "%f\n", integral);

        double i1 = 0.0, i2 = 0.0, x, y;
        if(!ihbg->first.compare("Other")) fprintf(fbackfits, "%s\t%s\t", "Other", "shape");
        else if(!ihbg->first.compare("t#bar{t}")) fprintf(fbackfits, "%s\t%s\t", "TT", "shape");
        else if(!ihbg->first.compare("Z+Jets")) fprintf(fbackfits, "%s\t%s\t", "ZJ", "shape");
        else fprintf(fbackfits, "%s\t%s\t", ihbg->first.c_str(), "shape");
        for(int i = istart; i <= h->GetNbinsX(); i++)
        {
            if(fmod(h->GetBinLowEdge(i), 200) == 0.0)
            {
                if(h->GetBinLowEdge(i) > 600) fprintf(fbackfits, "%f\t", 2 * i1 / (i1 + i2));
                i1 = 0.0;
                i2 = 0.0;
            }
            tgtotal->GetPoint(i - istart, x, y);
            i1 += y + tgtotal->GetErrorYhigh(i - istart);
            i2 += y - tgtotal->GetErrorYlow(i - istart);
        }
        fprintf(fbackfits, "%f\n", 2 * i1 / (i1 + i2));

        TH1 *dummy = new TH1F("dummy", "dummy", 1000, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()));
        dummy->GetXaxis()->SetTitle("M_{#mu#mujj} [GeV]");
        dummy->GetYaxis()->SetRangeUser(std::max(0.001, 0.5 * h->GetMinimum()), h->GetMaximum()*3);
        dummy->GetYaxis()->SetTitle("Events");
        dummy->GetYaxis()->SetTitleOffset(1.0);
        dummy->SetStats(0);
        dummy->SetTitle(0);

        TLatex mark;
        mark.SetTextSize(0.04);
        mark.SetTextFont(42);
        mark.SetNDC(true);

        char tmp4[128];
        sprintf(tmp4, "Integral (> 600 GeV): %0.3f", h->Integral(h->FindBin(600), h->GetNbinsX() + 1));

        dummy->Draw();
        tgtotal->Draw("E3 same");
        tgshape->Draw("E3 same");
        tgslope->Draw("E3Top same");
        h->Draw("same");

        TF1 *fu2 = new TF1("fu", "[0]*exp([1]*(x-[2]))", 600, 2800);
        fu2->SetParameters(f[0]->Eval(600), f[0]->GetParameter(1) + f[0]->GetParError(1), 600);
        TF1 *fd2 = new TF1("fd", "[0]*exp([1]*(x-[2]))", 600, 2800);
        fd2->SetParameters(f[0]->Eval(600), f[0]->GetParameter(1) - f[0]->GetParError(1), 600);
        fu2->SetLineWidth(2);
        fu2->SetLineStyle(2);
        fu2->SetMarkerStyle(0);
        fd2->SetLineWidth(2);
        fd2->SetLineStyle(2);
        leg->AddEntry(fu2, "slope stat err");
        fu2->Draw("same");
        fd2->Draw("same");

        if(drawtrialfuncs)
            for(unsigned int ifunc = 0; ifunc < fffs.size(); ifunc++)
            {
                f[ifunc]->Draw("same");
                leg->AddEntry(f[ifunc], fffs[ifunc].first.c_str());
            }
        else
        {
            f[0]->Draw("same");
            f[0]->SetFillColor(kWhite);
            //char tmp3[128];
            //sprintf(tmp3, "Fit (slope: %0.2e)", f[0]->GetParameter(1));
            leg->AddEntry(f[0], "Exp Fit");
        }
        leg->Draw("same");
        mark.DrawLatex(0.65, 0.95, "CMS Preliminary");
        mark.SetTextSize(0.025);
        mark.DrawLatex(0.645, 0.615, tmp4);

        delete [] f;
    }
    fclose(fbackfits);
}

void HnuPlots::plotQCDFit()
{
    setTDRStyle();

    TH1* hd = (TH1*)datahist.second->Clone();
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {
        hd->Add(ihbg->second, -1);
    }

    printf("Integral              : %f\n", hd->Integral(0, hd->GetNbinsX() + 1));
    printf("Integral above 600 GeV: %f\n", hd->Integral(hd->FindBin(600), hd->GetNbinsX() + 1));

    TF1* ff = new TF1("ff", "expo", 600, 2000);
    ff->FixParameter(1, -0.004937);
    hd->Fit(ff, "LLRQNMB");

    printf("Slope: %f +/- %f\n", ff->GetParameter(1), ff->GetParError(1));

    TText* CMSprelim = new TText(0.16, 0.96, "CMS Preliminary") ;
    CMSprelim->SetNDC();
    CMSprelim->SetTextSize(0.03);

    char tmp_iLumi[128];
    sprintf(tmp_iLumi, "%0.1f fb^{-1} at #sqrt{s} = 7 TeV", iLumi / 1000);
    TLatex* intLumi = new TLatex(0.72, 0.96, tmp_iLumi);
    intLumi->SetNDC();
    intLumi->SetTextSize(0.03);
    intLumi->SetLineWidth(2);

    TBox* box = new TBox(0, 0.001, 600, 2.);
    box->SetFillColor(kBlue);
    box->SetFillStyle(3003);

    TCanvas* ca = new TCanvas("ca", "ca");
    ca->Clear();
    ca->cd(0)->SetLogy();
    hd->SetLineColor(kBlack);
    hd->GetXaxis()->SetTitle("M(#mu#mujj) [GeV]");
    hd->GetYaxis()->SetTitle("Events / (200 GeV)");
    hd->GetYaxis()->SetTitleOffset(1.5);
    hd->SetMaximum(2.0);
    hd->SetMinimum(0.001);
    hd->SetStats(true);
    hd->Draw();
    ff->Draw("same");
    fixOverlay();
    //fitbox->Draw("same");
    box->Draw("same");
    CMSprelim->Draw("same");
    intLumi->Draw("same");

    //    TCanvas *c1 = new TCanvas("ca", "ca", 600, 600);
    //    c1->SetLogy(true);
    //    hd->Draw();
    //    ff->Draw("same");
}

void HnuPlots::plotNorm(double lower, double upper)
{
    if(bghists.size() < 2)
    {
        std::cout << "TOO FEW BG HISTOGRAMS!" << std::endl;
        return;
    }
    //TH1* hdata = (TH1*)datahist.second->Clone();
    TH1* hrenorm = (TH1*)bghists[0].second->Clone();
    TH1* hbg = (TH1*)bghists[1].second->Clone();
    TH1* hdatamod = (TH1*)datahist.second->Clone();
    hdatamod->Add(hbg, -1);

    double dataerror, mcerror;
    double dataint = hdatamod->IntegralAndError(hdatamod->FindBin(lower), hdatamod->FindBin(upper), dataerror);
    double mcint = hrenorm->IntegralAndError(hrenorm->FindBin(lower), hrenorm->FindBin(upper), mcerror);

    std::cout << "Normalization factor: " << dataint << "/" << mcint << " = " << dataint / mcint << " +/- "
            << dataint / mcint * sqrt(dataerror * dataerror / (dataint * dataint) + mcerror * mcerror / (mcint * mcint)) << std::endl;

    bghists[0].second->Scale(dataint / mcint);
    plot();
}

void HnuPlots::cutFlow()
{
    std::vector<std::pair<std::string, TH1*> > hists;
    hists.push_back(datahist);
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {
        hists.push_back(*ihbg);
    }
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        hists.push_back(*ihbg);
    }

    printf("%18s & ", "cut level");
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator i = hists.begin(); i != hists.end(); i++)
    {
        if(i != hists.end() - 1) printf("%10s & ", i->first.c_str());
        else printf("%10s \\\\ \\hline\n", i->first.c_str());
    }
    for(int icl = 1; icl < 9; icl++)
    {
        for(std::vector<std::pair<std::string, TH1*> >::const_iterator i = hists.begin(); i != hists.end(); i++)
        {
            if(i == hists.begin()) printf("%18s & ", i->second->GetXaxis()->GetBinLabel(icl));
            if(i->second->GetBinContent(icl) > 3 || i == hists.begin())
            {
                if(i != hists.end() - 1) printf("%10.0f & ", i->second->GetBinContent(icl)*4680.0 / 2510);
                else printf("%10.0f \\\\ \\hline\n", i->second->GetBinContent(icl));
            }
            else
            {
                if(i != hists.end() - 1) printf("%10.3f & ", i->second->GetBinContent(icl)*4680.0 / 2510);
                else printf("%10.3f \\\\ \\hline\n", i->second->GetBinContent(icl)*4680.0 / 2510);
            }
        }
    }
}

void HnuPlots::integrals(double min, double max)
{
    printf("Integrals\n%s: %f\n", datahist.first.c_str(), datahist.second->Integral(datahist.second->FindBin(min), datahist.second->FindBin(max)));
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {
        printf("%s: %f\n", ihbg->first.c_str(), ihbg->second->Integral(ihbg->second->FindBin(min), ihbg->second->FindBin(max)));
    }
    for(std::vector<std::pair<std::string, TH1*> >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        printf("%s: %f\n", ihbg->first.c_str(), ihbg->second->Integral(ihbg->second->FindBin(min), ihbg->second->FindBin(max)));
    }
}

std::string HnuPlots::getHistogramXAxisTitle(std::string name)
{
    if(!name.compare("mWR")) return "M_{#mu#mujj} [GeV]";
    if(!name.compare("mMuMu")) return "M_{#mu#mu} [GeV]";
    if(!name.compare("mMuMuZoom")) return "M_{#mu#mu} [GeV]";
    if(!name.compare("mNuR1")) return "M_{N_{#mu1}} [GeV]";
    if(!name.compare("mNuR2")) return "M_{N_{#mu2}} [GeV]";
    if(!name.compare("mJJ")) return "M_{jj} [GeV]";
    return "Sorry, no approperiate label found";
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

void HnuPlots::setFormLabel(std::string fl)
{
    formlabel = fl;
}

void HnuPlots::setXRange(double min, double max)
{
    xmin = min;
    xmax = max;
}

void HnuPlots::setSavePlots(bool sp)
{
    saveplots = sp;
}

const std::string cutlevels[] = {
    "cut0_none",
    "cut1_LLJJpt",
    "cut2_TrigMatches",
    "cut3_Vertex",
    "cut4_Mu1HighPt",
    "cut5_diLmass",
    "cut6_mWRmass",
    "cutX_LLpt",
    "cut5a_loDiLmass",
    ""
};

void plotMCVar(int cutlevel, std::string name, std::string xaxis = "M_{W_{R}} [GeV]")
{
    using namespace std;

    //background legend label, TFile
    vector<HnuPlots::FileStruct> bg, sig;
    //HnuPlots::FileStruct
    //bg.push_back(HnuPlots::FileStruct("DY madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuAnalysis_DYJetsToLL_madgraph_2140ipb.root"),                    "hNuMu40/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("DY sherpa",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles//local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2140ipbKFACTOR.root"), "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("t#bar{t} madgraph",  new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles/heavyNuTopAnalysis_TTJets_TuneZ2_7TeV-madgraph-tauola_oct21_2140ipb.root"),                    "hNuTopHiLumi/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("t#bar{t} powheg",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles//local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_Top_2140ipb.root"),                    "hNuTopHiLumi/" + cutlevels[cutlevel] + "/" + name));
    bg.push_back(HnuPlots::FileStruct("t#bar{t} madgraph", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6.root"), "hNu/" + cutlevels[cutlevel] + "/" + name));
    bg.push_back(HnuPlots::FileStruct("t#bar{t} powheg", new TFile("/local/cms/user/dahmes/wr2011/bgMC/Summer11/aug30/ttbar-PFJets.root"), "hNu/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("e#mu Data",    new TFile("/local/cms/user/dahmes/wr2011/data_2140ipb/run2011A-top.root"),                    "hNuTop/" + cutlevels[cutlevel] + "/" + name));
    //bg.push_back(HnuPlots::FileStruct("W+Jets",    new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/Tools/rootfiles//local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2140ipb.root"),                    "hNu/" + cutlevels[cutlevel] + "/" + name));
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

void plotFNALVar(int cutlevel, std::string name, std::string xaxis = "M_{W_{R}} [GeV]")
{
    using namespace std;

    //background legend label, TFile
    vector<HnuPlots::FileStruct> bg, sig;
    //HnuPlots::FileStruct
    bg.push_back(HnuPlots::FileStruct("PF2PAT + Refit",  new TFile("/local/cms/user/pastika/heavyNuAnalysis_2011/Summer11/heavyNuBasePFAnalysis_noAug5.root"), "hNuFNALwPF2PAT/diLcut_no70to100/" + name));
    bg.push_back(HnuPlots::FileStruct("PFJet + Tracker",  new TFile("/local/cms/user/pastika/heavyNuAnalysis_2011/Summer11/heavyNuBasePFAnalysis_noAug5.root"), "hNuFNAL/diLcut_no70to100/" + name));
    //data
    HnuPlots::FileStruct data("Data", new TFile("/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/42X/2140ipb/RecoEff/oct18/run2011A-recoEff-2140ipb.root"), "hNu/" + cutlevels[cutlevel] + "/" + name);

    HnuPlots hps(data, bg, sig);
    hps.setXAxisTitle(xaxis.c_str());
    hps.setYAxisTitle("Events");
    hps.setRebin(2);
    hps.plotMCComp();
}

void setBgandData(bool is2011A, bool is2011B, HnuPlots::FileStruct& data, std::vector<std::vector<HnuPlots::FileStruct> >& bg, double& lumi, int cutlevel = 5, std::string plot = "mWR")
{
    char fdata[128];

    //background
    std::vector<HnuPlots::FileStruct> bgTT, bgZJ, bgOther;
    if(is2011A)
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_dec14.root"),                   "hNuMu24/" + cutlevels[cutlevel] + "/" + plot,  216.2,  3048, 1.545,  "hNuMu24/cutlevel"));
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_dec14.root"),                   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  3048, 1.545,  "hNuMu40/cutlevel"));
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),           "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,   17.4, 1.05,  "hNuMu24/cutlevel"));
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),           "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  17.4, 1.05,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.3, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,     43, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,    43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,   18.2, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.9, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    if(is2011B)
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_dec14.root"),                   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 3048,  1.534,  "hNuMu40/cutlevel"));
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_dec14.root"),           "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,  17.4, 0.94,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"),    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,    43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,  18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011b-dec23.root");
        lumi += 2510.5;
    }
    if(is2011A && is2011B) sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
    bg.push_back(bgTT);
    bg.push_back(bgZJ);
    bg.push_back(bgOther);

    //data
    data.label = "Data";
    data.file = new TFile(fdata);
    data.histpath = "hNu/" + cutlevels[cutlevel] + "/" + plot;
}

void plot2011(bool is2011A = true, bool is2011B = true, int cutlevel = 5, std::string plot = "mWR", bool log = true, double xmin = 0.0, double xmax = 0.0)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig, vsig2;
    setBgandData(is2011A, is2011B, data, bg, lumi, cutlevel, plot);

    //signal
    vsig.push_back(HnuPlots::FileStruct("M_{W_{R}} = 1.8 TeV", new TFile("/local/cms/user/jmmans/heavyNu_signal/signal_1800_1000.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.01252, 1.27,  "hNuMu40/cutlevel"));
    vsig2.push_back(HnuPlots::FileStruct("M_{W_{R}} = 1.0 TeV", new TFile("/local/cms/user/jmmans/heavyNu_signal/signal_1000_600.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.3573, 1.31,  "hNuMu40/cutlevel"));
    //sig.push_back(vsig2);
    sig.push_back(vsig);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("please auto set the axis");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(2);
    hps.setLog(log);
    hps.setXRange(xmin, xmax);
    if(is2011A) hps.setFormLabel("hNu_2011A");
    if(is2011B) hps.setFormLabel("hNu_2011B");
    if(is2011A && is2011B) hps.setFormLabel("hNu_2011Combined");
    hps.plot();
}

void plotMCFits(bool is2011A = true, int cutlevel = 5, bool log = true)
{
    using namespace std;

    char plot[] = "mWR";
    double lumi;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    setBgandData(is2011A, !is2011A, data, bg, lumi, cutlevel, plot);

    std::string bgestfname = "bgest";
    if(is2011A) bgestfname += "2011A.txt";
    else  bgestfname += "2011B.txt";

    HnuPlots hps(data, bg, sig, lumi);
    hps.setRebin(1);
    hps.setLog(log);
    if(is2011A) hps.setFormLabel("bgfits_2011A");
    else hps.setFormLabel("bgfits_2011B");
    hps.plotMCShape(bgestfname);
}

void plotTTBarNorm(bool is2011A = true, bool is2011B = true, int cutlevel = 5, bool log = true)
{
    using namespace std;

    char plot[] = "mWR", fdata[256];
    double lumi = 0.0;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgTT, bgOther;
    if(is2011A)
    {
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),              "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  17.4, 1.0,  "hNuMu24/cutlevel"));
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),              "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 17.4, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                      "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  3048,  1.545, "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                      "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 3048,  1.545, "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  5.3, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  43, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  18.2, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,  5.9, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    if(is2011B)
    {
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_dec14.root"),              "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 17.4, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_QCD_jan6.root"),                      "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 3048,  1.534, "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"),    "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"), "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other",    new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011b-dec23.root");
        lumi += 2510.5;
    }
    if(is2011A && is2011B) sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
    bg.push_back(bgTT);
    bg.push_back(bgOther);

    //data
    HnuPlots::FileStruct data("Data", new TFile(fdata), "hNuTop/" + cutlevels[cutlevel] + "/" + plot);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("M_{e#mujj} [GeV]");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(2);
    hps.setLog(log);
    if(is2011A && !is2011B) hps.setFormLabel("ttnorm_2011A");
    else if(!is2011A && is2011B) hps.setFormLabel("ttnorm_2011B");
    else hps.setFormLabel("ttnorm_2011Combined");
    hps.plotNorm(20.0, 3000.0);
}

void plotZJNorm(bool is2011A = true, bool is2011B = true, int cutlevel = 4, bool log = true)
{
    using namespace std;

    char plot[] = "mMuMuZoom", fdata[256];
    double lumi = 0.0;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther;
    if(is2011A)
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_dec14.root"),                   "hNuMu24/" + cutlevels[cutlevel] + "/" + plot,  216.2,  3048, 1.0, "hNuMu24/cutlevel"));
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_dec14.root"),                   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  3048, 1.0, "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),         "hNuMu24/" + cutlevels[cutlevel] + "/" + plot,  216.2,  17.4, 1.051,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),         "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  17.4, 1.051,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.3, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.3, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,     43, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,    43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,   18.2, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,  18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu24/" + cutlevels[cutlevel] + "/" + plot, 216.2,    5.9, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1956.7,   5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    if(is2011B)
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_dec14.root"),                   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 3048,  1.0, "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_dec14.root"),         "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5, 17.4, 0.940,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"),    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"), "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,    43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,  18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 2510.5,   5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011b-dec23.root");
        lumi += 2510.5;
    }
    if(is2011A && is2011B) sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
    bg.push_back(bgZJ);
    bg.push_back(bgOther);

    //data
    HnuPlots::FileStruct data("Data", new TFile(fdata), "hNu/" + cutlevels[cutlevel] + "/" + plot);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("please auto set the axis");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(20);
    hps.setLog(log);
    if(is2011A && !is2011B) hps.setFormLabel("zjnorm_2011A");
    else if(!is2011A && is2011B) hps.setFormLabel("zjnorm_2011B");
    else hps.setFormLabel("zjnorm_2011Combined");
    hps.setXRange(60.0, 500.0);
    hps.plotNorm(60.0, 120.0);
}

void plotCutFlow(bool is2011A = true, bool is2011B = true)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig;
    setBgandData(is2011A, is2011B, data, bg, lumi, 9, "cutlevel");

    //signal
    vsig.push_back(HnuPlots::FileStruct("M_{W_{R}} = 1.8 TeV", new TFile("/local/cms/user/jmmans/heavyNu_signal/signal_1800_1000.root"), "hNuMu40/cutlevel", lumi, 0.01252, 1.27,  "hNuMu40/cutlevel"));
    sig.push_back(vsig);

    HnuPlots hps(data, bg, sig, lumi);
    hps.cutFlow();
}

void plotQCD(bool is2011A = true, std::string cutlevel = "diLmassCuts", std::string plot = "mWR", bool log = true, double xmin = 0.0, double xmax = 0.0)
{
    using namespace std;

    char fdata[128];
    double lumi = 0.0;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgTT, bgZJ, bgOther, vsig;
    if(is2011A)
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                   "hNuQCDMu24/" + cutlevel + "/" + plot,  216.2, 3160,  1.49 * 0.02269,  "hNuMu24/cutlevel"));
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                   "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 3160,  1.49 * 0.02269,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_QCD_dec27.root"),           "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  17.4, 1.06,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("t#bar{t}", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_QCD_dec27.root"),           "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 17.4, 1.06,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root"),    "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root"),    "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root"), "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root"), "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  43, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  18.2, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.9, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    else
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_QCD_jan6.root"),                    "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 3160, 1.48 * 0.02269 * 1.16,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_QCD_dec27.root"),              "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 17.4, 0.95,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2011B_QCD_jan9.root"),                            "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 31314, 1.48 * 0.02269 * 1.16,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_QCD_dec27.root"),    "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_QCD_dec27.root"), "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root"),                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
        lumi += 2510.5 + 216.2 + 1956.7;
    }
    //bg.push_back(bgTT);
    bg.push_back(bgZJ);
    bg.push_back(bgOther);

    //data
    //HnuPlots::FileStruct data("Data", new TFile("data-run2011a-run2011b-nov11.root"), "hNu/" + cutlevels[cutlevel] + "/" + plot);
    HnuPlots::FileStruct data("Data", new TFile(fdata), "hNuQCD/" + cutlevel + "/" + plot);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("please auto set the axis");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(5);
    hps.setLog(log);
    hps.setXRange(xmin, xmax);
    if(is2011A) hps.setFormLabel("QCD_2011A");
    else hps.setFormLabel("QCD_2011B");
    hps.setSavePlots(false);
    hps.plot(); //Norm(85, 95);
    hps.integrals(200, 2000);
    hps.plotQCDFit();
}

void plotRussianTTBarNorm(bool is2011A = true, bool is2011B = true, int cutlevel = 5, double cl = 0.0, double ch = 2801.0, std::string plot = "mWRvsminLPt", bool log = true)
{
    using namespace std;

    char fdata[256];
    double lumi = 0.0;

    const string rcuts[] = {
        "Mu1Pt30GeV",
        "Mu1Pt40GeV",
        "Mu1Pt50GeV",
        "Mu1Pt60GeV",
        "Mu1Pt80GeV",
        "Mu1Pt100GeV"
    };

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgTT, bgOther;
    if(is2011A)
    {
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),              "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  17.4, 1.0,  "hNuMu24/cutlevel", cl, ch));
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_dec14.root"),              "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 17.4, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                   "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  3048,  1.545, "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root"),                   "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 3048,  1.545, "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  5.3, 1.0,  "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"),    "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 5.3, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_dec14.root"), "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  43, 1.0,  "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 43, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  18.2, 1.0,  "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu24/" + rcuts[cutlevel] + "/" + plot, 216.2,  5.9, 1.0,  "hNuMu24/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel", cl, ch));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    if(is2011B)
    {
        bgTT.push_back(   HnuPlots::FileStruct("t#bar{t}",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_dec14.root"),              "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 17.4, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Background", new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_QCD_jan6.root"),                   "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 3048,  1.534, "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"),    "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 5.3, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_dec14.root"), "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 5.3,  1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 43, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 18.2, 1.0,  "hNuMu40/cutlevel", cl, ch));
        bgOther.push_back(HnuPlots::FileStruct("Other",      new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_dec14.root"),                "hNuTopMu40/" + rcuts[cutlevel] + "/" + plot, 2510.5, 5.9, 1.0,  "hNuMu40/cutlevel", cl, ch));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011b-dec23.root");
        lumi += 2510.5;
    }
    if(is2011A && is2011B) sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
    bg.push_back(bgTT);
    bg.push_back(bgOther);

    //data
    HnuPlots::FileStruct data("Data", new TFile(fdata), "hNuTop/" + rcuts[cutlevel] + "/" + plot);
    data.clow = cl;
    data.chigh = ch;

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("M_{e#mujj} [GeV]");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(2);
    hps.setLog(log);
    hps.setSavePlots(false);
    if(is2011A && !is2011B) hps.setFormLabel("ttnorm_2011A");
    else if(!is2011A && is2011B) hps.setFormLabel("ttnorm_2011B");
    else hps.setFormLabel("ttnorm_2011Combined");
    hps.plotNorm(20.0, 3000.0);
}








