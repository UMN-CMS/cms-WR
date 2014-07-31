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
#include "TMath.h"
#include "TExec.h"
#include "TSpline.h"
#include "TGraphSmooth.h"
#include "Math/VectorUtil.h"

#include "tdrstyle.C"
//#include "fitHNBackground.cc"
#include "../AnalysisModules/src/HeavyNuTree.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include <TLorentzVector.h>

const int colors[] = {
    //kCyan + 2,
    //kOrange + 1,
    //kGreen - 3,
    kAzure,
    kOrange + 7,
    kGreen + 2,
    //kBlue,
    //kRed,
    //kGreen,
    //kBlue + 2,
    //kRed,
    //kGreen - 2,
    //kBlue + 2,
    //kRed + 2,
    //kGreen + 2,
    kYellow + 4,
    kMagenta - 1,
    kRed,
    kBlue,
    kGreen
};
const int NCOLORS = sizeof(colors) / sizeof(int);

const int shcolors[] = {
    kRed,
    //kOrange+7,
    kBlue,
    kGreen,
    kYellow,
    kMagenta,
    kOrange
};
const int NSHCOLORS = sizeof(shcolors) / sizeof(int);

const int hatchs[] = {
    3454,
    3445,
    3002,
    3013,
    3006,
    3007,
    3017,
    3018
};
const int NHATCHS = sizeof(hatchs) / sizeof(int);

double bins[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 4.0};
double systBins[] =            {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 4.0};
//double systBins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 4000.0};
//double bins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 4000.0};
//double bins[] = {600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 2200.0, 4000.0};
//double bins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 2000.0, 4000.0};

class HnuPlots
{
public:

    struct FileStruct
    {
        std::string label;
        std::string file;
        std::string histpath;
        double intLumi, cs, kfactor;
        std::string normhistpath;
        double clow, chigh;
        bool projX;
        int normbin;
        bool unh;
        double normll, normul;
        bool loadtuple, histFromTuple;
        int cutlevel;
        double thll, thul, thb;
        HeavyNuTree::HNuSlopeFitInfo *tpll, *tpul;
        bool smooth_hist;

        FileStruct();
        FileStruct(std::string l, std::string f, std::string h);                            //                   "",             0.0,             0.0,           true,          1,           true,              0.0,              0.0,                  lt,                hft,               0,               0,                -1,                                   &ll,                                    &ul
        FileStruct(std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh = "", double cl = 0.0, double ch = 0.0, bool px = true, int nb = 1, bool un = true, double nll = 0.0, double nul = 0.0, bool loadtuple = false, bool lhft = false, double ll = 0.0, double ul = 1.0, double bn = 1.0, HeavyNuTree::HNuSlopeFitInfo* tll = 0, HeavyNuTree::HNuSlopeFitInfo* tul = 0, bool smooth = false);
        //FileStruct(std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh, double nll = 0.0, double nul = 0.0, bool loadtuple = false, bool lhft = false, double ll = 0.0, double ul = 1.0, double bn = 1.0);
        //, double cl = 0.0, double ch = 0.0, bool px = true, int nb = 1, bool un = true,
        FileStruct(bool loadtuple, std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh, int cutlevel, int nb = 1);
        ~FileStruct();
    } ;

    struct HistStruct
    {
        std::string label;
        TH1 *hist, *normhist;
        double normll, normul;
        std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> > fittree;
        bool smooth_hist;

        HistStruct();
        HistStruct(std::string l, TH1* h, TH1* nh = NULL, double nll = 0.0, double nul = 0.0, bool smooth = false);
    } ;

    struct Limits
    {
        double thll, thul, nb;

        Limits() {}
        Limits(double ll, double ul, double b) : thll(ll), thul(ul), nb(b){ }
    } ;

    HnuPlots(){ }
    HnuPlots(FileStruct& fdata, std::vector<std::vector<FileStruct> >& vfbg, std::vector<std::vector<FileStruct> >& vfsig, double iL);
    void plot();
    void plot1D();
    void plot2D();
    void plotMCComp(bool rescale = false);
    void plotMCShape(std::string bgfilename);
    //void plotMCShapeUnbined(vector<FitHNBackground::outputData*> opdtmp);
    void plotQCDFit();
    void plotNorm(double lower, double upper, bool flip = false);
    void plotRatios();
    void scaleByShape(double llow, double lhigh, int npar = 2);
    void cutFlow();
    void sigEff();
    void integrals(double min, double max, double* passnum = NULL, double* err = NULL);
    void sigRMS();
    void sigStatErr();
    void sigMatch();
    void mcBgShape(int cutlevel = 5, std::string sample = "");
    void mcSystCalc(int forceBins = false, std::map<std::string, std::vector<double> >* systMap = NULL);
    void loadSystFile(std::string systfile, std::string ratefile, bool includeBG = false);
    void setRebin(int rbval);
    void setXAxisTitle(std::string label);
    void setYAxisTitle(std::string label);
    void setAutoSort(bool as);
    void setLog(bool log);
    void setCompPlot(bool cp);
    void setFormLabel(std::string);
    void setXRange(double min, double max);
    void setYRange(double min, double max);
    void setSavePlots(bool sp);
    void autoSetHistogramAxisTitle(int mode = 0);

private:
    std::vector<HistStruct> bghists;
    std::vector<HistStruct> sighists;
    HistStruct datahist;

    std::vector<std::vector<float> > systematics;
    std::vector<float> shapeerr;

    int rebin, nhist;
    double iLumi, xmin, xmax, ymin, ymax;
    std::string xaxislabel;
    std::string yaxislabel;
    std::string formlabel;
    bool autosort, islog, saveplots, plotSMoData;
    double sigscale;

    TH1* project(TH2* h2d, double cl, double ch, bool porjx = true);
    int projcount;

    TH1* histFromTuple(std::string histpath, double nb, std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >& bgtvec, HeavyNuTree::HNuSlopeFitInfo *ll = NULL, HeavyNuTree::HNuSlopeFitInfo *ul = NULL, bool smooth = false);
    bool runFilter(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iE);
    void histFromDataCard(std::map<std::pair<std::string, std::string>, std::vector<float> >& uncerts);
    double getTupleVar(std::string var, const HeavyNuTree::HNuSlopeFitInfo& ll);
    bool dynamicalCut(double var, double cut, char cutType);

    class fitfunction
    {
    public:
        std::vector<HistStruct> bghs;
        int npar;
        double operator()(double * x, double * par);
    } ;
} ;

bool compHistInt(HnuPlots::HistStruct h1, HnuPlots::HistStruct h2)
{
    return h1.hist->Integral() < h2.hist->Integral();
}

double HnuPlots::fitfunction::operator()(double * x, double * par)
{
    double retval = 0.0;
    int n = 0;
    for(std::vector<HnuPlots::HistStruct>::const_iterator ihbg = bghs.begin(); ihbg != bghs.end(); ++ihbg)
    {
        retval += ((n < npar)?par[n]:1.0) * ihbg->hist->GetBinContent(ihbg->hist->FindBin(x[0]));
        n++;
    }
    return retval;
}

HnuPlots::HnuPlots(FileStruct& fdata, std::vector<std::vector<HnuPlots::FileStruct> >& vfbg, std::vector<std::vector<HnuPlots::FileStruct> >& vfsig, double iL)
{
    using namespace std;

    bool first;

    TH1::AddDirectory(kFALSE);

    int iColor = 0;
    for(vector<vector<HnuPlots::FileStruct> >::const_iterator ivbg = vfbg.begin(); ivbg != vfbg.end(); ++ivbg)
    {
        first = true;
        std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> > bgtvec;
        for(vector<HnuPlots::FileStruct>::const_iterator ibgf = ivbg->begin(); ibgf != ivbg->end(); ++ibgf)
        {
            TFile * file = new TFile(ibgf->file.c_str());
            //ROOT-FORTRAN magic so that my histograms are not associated with their TFiles
            //gDirectory->cd("");

            std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> > ibgtvec;

            TH1 * h = NULL;
            TH1 * hn = (TH1*)file->Get(ibgf->normhistpath.c_str());
            //if(hn) hn = (TH1*)hn->Clone();
            // Load Tupple if needed
            if(ibgf->loadtuple)
            {
                HeavyNuTree* hnt;
                TDirectory* tdir = (TDirectory*)file->Get((ibgf->histpath.substr(0, ibgf->histpath.find("/"))).c_str());
                hnt = new HeavyNuTree(*tdir, false);

                double scale = 0.0;
                if(ibgf->unh && hn && ibgf->normbin >= 0) scale = ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->GetBinContent(ibgf->normbin);
                else if(ibgf->unh && hn) scale = ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->Integral(0, hn->GetNbinsX() + 1);
                else scale = ibgf->intLumi * ibgf->cs * ibgf->kfactor;

                do
                {
                    if((ibgf->file.compare("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_2/Mu-Run2012ABCD-22Jan2013-v1.root") == 0) && (hnt->event_.mlljj > 1000 && hnt->event_.mlljj < 1200)) hnt->event_.weight *= 1.2294;
                    //if((ibgf->label.compare("t#bar{t}") == 0) && (hnt->event_.mlljj > 1000 && hnt->event_.mlljj < 1200)) hnt->event_.weight *= 1.2294;

                    bgtvec.push_back(std::make_pair(hnt->event_, scale));
                    ibgtvec.push_back(std::make_pair(hnt->event_, scale));
                }
                while(hnt->GetNextEvent());
                delete [] hnt;
            }

            // gethistogram
            if(ibgf->loadtuple && ibgf->histFromTuple)
            {
                h = histFromTuple(ibgf->histpath.substr(ibgf->histpath.rfind("/") + 1, ibgf->histpath.size()), ibgf->thb, ibgtvec, ibgf->tpll, ibgf->tpul);
            }
            else if(fabs(ibgf->clow) < 1e-300 && fabs(ibgf->chigh) < 1e-300)
            {
                h = (TH1*)file->Get(ibgf->histpath.c_str());
                //if(h) h = (TH1*)h->Clone();
            }
            else
            {
                TH2 *h2 = (TH2*)file->Get(ibgf->histpath.c_str());
                if(!h2) std::cout << "failed to get File:hist - " << file->GetName() << " : " << ibgf->histpath << std::endl;
                //else h2 = (TH2*)h2->Clone();
                h = project(h2, ibgf->clow, ibgf->chigh);
            }

            // Scale histogram to proper lumonisity
            if(h)
            {
                if(ibgf->unh && hn && ibgf->normbin >= 0) h->Scale(ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->GetBinContent(ibgf->normbin));
                else if(ibgf->unh && hn) h->Scale(ibgf->intLumi * ibgf->cs * ibgf->kfactor / hn->Integral(0, hn->GetNbinsX() + 1));
                else h->Scale(ibgf->intLumi * ibgf->cs * ibgf->kfactor);
                if(first)
                {
                    bghists.push_back(HnuPlots::HistStruct(ibgf->label, (TH1*)h->Clone(), hn?((TH1*)hn->Clone()):NULL, ibgf->normll, ibgf->normul));
                    first = false;
                }
                else
                {
                    bghists.back().hist->Add(h);
                    if(bghists.back().normhist && hn) bghists.back().normhist->Add(hn);
                }
            }
            else
            {
                if(!h) std::cout << "failed to get File:hist - " << file->GetName() << " : " << ibgf->histpath << std::endl;
                if(!hn) std::cout << "failed to get File:hist - " << file->GetName() << " : " << ibgf->normhistpath << std::endl;
            }

            file->Close();
        }
        if(bgtvec.size() > 0) bghists.back().fittree = bgtvec;
        bghists.back().hist->SetFillColor(colors[iColor % NCOLORS]);
        bghists.back().hist->SetFillStyle(hatchs[iColor % NHATCHS]);
        bghists.back().hist->SetLineColor(colors[iColor % NCOLORS]);
        bghists.back().hist->SetMarkerColor(colors[iColor % NCOLORS]);
        bghists.back().hist->SetMarkerStyle(0);
        bghists.back().hist->SetLineWidth(1);
        iColor++;

        //TFile * fdfd = new TFile("fbg.root", "RECREATE");
        //TH1 * hcopy = (TH1*)bghists.back().hist->Clone("mWR");
        //TDirectory * td1 = fdfd->mkdir("hNuMu");
        //TDirectory * td2 = td1->mkdir("cut5_diLmass");
        //td2->cd();
        //hcopy->Write();
        //fdfd->Close();
    }

    iColor = 0;
    for(vector<vector<HnuPlots::FileStruct> >::const_iterator ivsig = vfsig.begin(); ivsig != vfsig.end(); ++ivsig)
    {
        first = true;
        for(vector<HnuPlots::FileStruct>::const_iterator isigf = ivsig->begin(); isigf != ivsig->end(); ++isigf)
        {
            TFile * file = new TFile(isigf->file.c_str());
            //ROOT-FORTRAN magic so that my histograms are not associated with their TFiles
            //gDirectory->cd("Rint:/");

            std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> > sigtvec;

            TH1 * hn = (TH1*)file->Get(isigf->normhistpath.c_str());
            //if(hn) hn = (TH1*)hn->Clone();
            // Load Tupple if needed
            if(isigf->loadtuple)
            {
                HeavyNuTree* hnt;
                TDirectory* tdir = (TDirectory*)file->Get((isigf->histpath.substr(0, isigf->histpath.find("/"))).c_str());
                std::cout << (isigf->histpath.substr(0, isigf->histpath.find("/"))) << std::endl;
                hnt = new HeavyNuTree(*tdir, false);

                double scale = 0.0;
                if(isigf->unh && hn && isigf->normbin >= 0) scale = isigf->intLumi * isigf->cs * isigf->kfactor / hn->GetBinContent(isigf->normbin);
                else if(isigf->unh && hn) scale = isigf->intLumi * isigf->cs * isigf->kfactor / hn->Integral(0, hn->GetNbinsX() + 1);
                else scale = isigf->intLumi * isigf->cs * isigf->kfactor;

                do
                {
                    sigtvec.push_back(std::make_pair(hnt->event_, scale));
                }
                while(hnt->GetNextEvent());
                delete [] hnt;
            }

            TH1 * h = 0;
            if(isigf->loadtuple && isigf->histFromTuple)
            {
                std::cout << isigf->histpath.substr(isigf->histpath.rfind("/") + 1, isigf->histpath.size()) << std::endl;
                h = histFromTuple(isigf->histpath.substr(isigf->histpath.rfind("/") + 1, isigf->histpath.size()), isigf->thb, sigtvec, isigf->tpll, isigf->tpul, isigf->smooth_hist);
            }
            else
            {
                h = (TH1*)file->Get(isigf->histpath.c_str());
            }
            //if(h) h = (TH1*)h->Clone();
            if(h)
            {
                if(hn) h->Scale(isigf->intLumi * isigf->cs * isigf->kfactor / hn->GetBinContent(isigf->normbin));
                else   h->Scale(isigf->intLumi * isigf->cs * isigf->kfactor);
                if(first)
                {
                    sighists.push_back(HnuPlots::HistStruct(isigf->label, (TH1*)h->Clone(), NULL,  0.0,  0.0, isigf->smooth_hist));
                    first = false;
                }
                else
                {
                    sighists.back().hist->Add(h);
                }
            }
            else
            {
                if(!h) std::cout << "failed to get File:hist - " << file->GetName() << " : " << isigf->histpath << std::endl;
                if(!hn) std::cout << "failed to get File:hist - " << file->GetName() << " : " << isigf->normhistpath << std::endl;
            }

            file->Close();
        }
        //sighists.back().second->SetFillColor(kWhite + iColor);
        //sighists.back().second->SetFillStyle(0);
        sighists.back().hist->SetLineColor(shcolors[iColor % NSHCOLORS]);
        sighists.back().hist->SetMarkerColor(shcolors[iColor % NSHCOLORS]);
        sighists.back().hist->SetLineWidth(2.5);
        iColor++;
    }

    if(fdata.file.size())
    {
        TFile * file = new TFile(fdata.file.c_str());
        //ROOT-FORTRAN magic so that my histograms are not associated with their TFiles
        //gDirectory->cd("Rint:/");

        std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> > dtvec;

        if(fdata.loadtuple)
        {
            HeavyNuTree* hnt;
            TDirectory* tdir = (TDirectory*)file->Get((fdata.histpath.substr(0, fdata.histpath.find("/"))).c_str());
            hnt = new HeavyNuTree(*tdir, false);

            do
            {
                dtvec.push_back(std::make_pair(hnt->event_, 1.0));
            }
            while(hnt->GetNextEvent());
            delete [] hnt;
        }

        TH1 *h;
        if(fdata.loadtuple && fdata.histFromTuple)
        {
            h = histFromTuple(fdata.histpath.substr(fdata.histpath.rfind("/") + 1, fdata.histpath.size()), fdata.thb, dtvec, fdata.tpll, fdata.tpul);
            //TFile * fdfd = new TFile("fdata.root", "RECREATE");
            //TH1 * hcopy = (TH1*)h->Clone("mWR");
            //TDirectory * td1 = fdfd->mkdir("hNuE");
            //TDirectory * td2 = td1->mkdir("cut6_mWRmass");
            //td2->cd();
            //hcopy->Write();
            //fdfd->Close();
        }
        else if(fabs(fdata.clow) < 1e-300 && fabs(fdata.chigh) < 1e-300)
        {
            h = (TH1*)file->Get(fdata.histpath.c_str());
            //h = (TH1*)h->Clone();
        }
        else
        {
            //cout << "I AM HERE" << endl;
            TH2 *h2 = (TH2*)file->Get(fdata.histpath.c_str());
            if(!h2) std::cout << "failed to get File:hist - " << file->GetName() << " : " << fdata.histpath << std::endl;
            //else h2 = (TH2*)h2->Clone();
            h = project(h2, fdata.clow, fdata.chigh);
        }
        if(h) datahist = HnuPlots::HistStruct(fdata.label, (TH1*)h->Clone());

        else std::cout << "failed to get File:hist - " << file->GetName() << " : " << fdata.histpath << std::endl;
        datahist.hist->SetLineColor(kBlack);
        datahist.hist->SetMarkerColor(kBlack);
        datahist.hist->SetMarkerStyle(23);
        datahist.hist->SetFillStyle(0);
        datahist.hist->SetFillColor(kWhite);

        file->Close();
    }

    rebin = -1;
    xaxislabel = "";
    yaxislabel = "";
    formlabel = "hNu";
    autosort = false;
    islog = false;
    iLumi = iL;
    xmin = xmax = 0.0;
    ymin = ymax = 0.0;
    saveplots = true;
    plotSMoData = true;
    projcount = 0;
}

HnuPlots::HistStruct::HistStruct()
{

    label = "";
    hist = NULL;
    normhist = NULL;
}

HnuPlots::HistStruct::HistStruct(std::string l, TH1* h, TH1* nh, double nll, double nul, bool smooth)
{

    label = l;
    hist = h;
    normhist = nh;
    normll = nll;
    normul = nul;
    smooth_hist = smooth;
}

TH1* HnuPlots::project(TH2* h2d, double cl, double ch, bool projx)
{
    char tmp[128];
    TH1* h = NULL;

    if(projx)
    {
        sprintf(tmp, "%s_projx_%d", h2d->GetName(), projcount++);
        h = new TH1D(tmp, tmp, h2d->GetXaxis()->GetNbins() + 1, h2d->GetXaxis()->GetBinLowEdge(1), h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetNbins()));
        for(int i = h2d->GetYaxis()->FindBin(cl); i <= h2d->GetYaxis()->FindBin(ch); i++)
        {
            for(int j = 0; j < h2d->GetXaxis()->GetNbins() + 1; j++)
            {
                h->SetBinContent(j, h->GetBinContent(j) + h2d->GetBinContent(j, i));
                double be1 = h->GetBinError(j), be2 = h2d->GetBinError(j, i);
                h->SetBinError(j, sqrt(be1 * be1 + be2 * be2));
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
                h->SetBinError(j, sqrt(be1 * be1 + be2 * be2));
            }
        }
    }

    return h;
}

TH1* HnuPlots::histFromTuple(std::string histValues, double nb, std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >& bgtvec, HeavyNuTree::HNuSlopeFitInfo *ll, HeavyNuTree::HNuSlopeFitInfo *ul, bool smooth)
{
    std::vector<std::string> histQs;
    bool invertCuts = false;
    std::map<std::pair<std::string, char>, double> cuts;
    size_t cutStart = histValues.find(";");
    std::string histValName = histValues.substr(0, cutStart);
    // read variable names to plot
    for(size_t pos = 0, npos = 0; npos != size_t(-1);pos = npos + 1)
    {
        npos = histValName.find(':', pos + 1);
        histQs.push_back(histValName.substr(pos, npos - pos));
    }
    // read cut inversion
    if((invertCuts = ((histValues.find('!') != size_t(-1))))) histValues.erase(histValues.find('!'), 1);
    // read cuts to implament
    for(size_t pos = cutStart + 1, npos = 0; npos != size_t(-1);pos = npos + 1)
    {
        npos = histValues.find(';', pos + 1);
        std::string tmp = histValues.substr(pos, npos - pos);
        size_t sepPos = 0;
        char cutType = ' ';
        if((sepPos = tmp.find('>')) != size_t(-1)) cutType = '>';
        else if((sepPos = tmp.find('<')) != size_t(-1)) cutType = '<';
        else continue;
        std::string t1 = tmp.substr(0, sepPos), t2 = tmp.substr(sepPos + 1, (size_t(-1)));
        t1.erase(remove(t1.begin(),t1.end(),' '),t1.end());
        t2.erase(remove(t2.begin(),t2.end(),' '),t2.end());
        char vname[32];
        sscanf(t1.c_str(), "%s", vname);
        tmp = vname;
        double cutvalue;
        sscanf(t2.c_str(), "%lf", &cutvalue);
        cuts[std::make_pair(tmp, cutType)] = cutvalue;
    }

    std::vector<Limits> vlim;
    char hname[128], htitle[128];
    //sprintf(hname, "%s_%s_%d", histValues.c_str(), label.c_str(), nhist++);
    sprintf(hname, "%s", histValues.c_str());
    if(histQs.size() == 2)sprintf(htitle, "%s;%s;%s", histQs[0].c_str(), histQs[0].c_str(), histQs[1].c_str());
    for(std::vector<std::string>::const_iterator ihlabel = histQs.begin(); ihlabel != histQs.end(); ++ihlabel)
    {
        std::string histV = *ihlabel;
        if(nb < 0)
        {
            if(histV.compare("mWR") == 0) vlim.push_back(Limits(0.0, 4.0, 100));
            else if(histV.compare("st") == 0) vlim.push_back(Limits(0.0, 4000.0, 100));
            else if(histV.compare("mLL") == 0 || histV.compare("mJJ") == 0) vlim.push_back(Limits(0.0, 2000.0, 100));
            else if(histV.compare("mLLZoom") == 0) vlim.push_back(Limits(0.0, 2000.0, 1000));
            else if(histV.compare("mLLNorm") == 0) vlim.push_back(Limits(0.06, 0.5, 220));
            else if(histV.compare("mNuR1") == 0 || histV.compare("mNuR2") == 0) vlim.push_back(Limits(0.0, 3000.0, 150));
            else if(histV.compare("mOuR1") == 0 || histV.compare("mOuR2") == 0) vlim.push_back(Limits(0.0, 2500.0, 125));
            else if(histV.compare("ptL1") == 0 || histV.compare("ptL2") == 0) vlim.push_back(Limits(0.0, 1000.0, 100));
            else if(histV.compare("etaL1") == 0 || histV.compare("etaL2") == 0) vlim.push_back(Limits(-2.5, 2.5, 50));
            else if(histV.compare("phiL1") == 0 || histV.compare("phiL2") == 0) vlim.push_back(Limits(-3.1415926535, 3.1415926535, 30));
            else if(histV.compare("ptJ1") == 0 || histV.compare("ptJ2") == 0) vlim.push_back(Limits(0.0, 1000.0, 100));
            else if(histV.compare("etaJ1") == 0 || histV.compare("etaJ2") == 0) vlim.push_back(Limits(-5.0, 5.0, 100));
            else if(histV.compare("phiJ1") == 0 || histV.compare("phiJ2") == 0) vlim.push_back(Limits(-3.1415926535, 3.1415926535, 30));
            else if(histV.compare("pL1") == 0 || histV.compare("pL2") == 0) vlim.push_back(Limits(0.0, 2000.0, 100));
            else if(histV.compare("n_vertex") == 0) vlim.push_back(Limits(0.0, 50.0, 100));
            else if(histV.compare("seL1") == 0 || histV.compare("seL2") == 0 || histV.compare("rhL1") == 0 || histV.compare("rhL2") == 0) vlim.push_back(Limits(0.0, 2000.0, 200));
            else if(histV.compare("gR9L1") == 0 || histV.compare("gR9L2") == 0) vlim.push_back(Limits(0.0, 2000.0, 200));
            else if(histV.compare("rhoScL1") == 0 || histV.compare("rhoScL2") == 0) vlim.push_back(Limits(0.0, 2.0, 100));
            else if(histV.compare("met") == 0) vlim.push_back(Limits(0.0, 1000.0, 100));
            else if(histV.compare("dEtaL") == 0 || histV.compare("dEtaJ") == 0) vlim.push_back(Limits(-5.0, 5.0, 100));
            else if(histV.compare("dPhiL") == 0 || histV.compare("dPhiJ") == 0) vlim.push_back(Limits(-3.1415926535, 3.1415926535, 100));
            else if(histV.compare("run") == 0) vlim.push_back(Limits(190000, 210000, 1000));
            else if(histV.compare("cutlevel") == 0) vlim.push_back(Limits(-1.5, 9.5, 11));
            else if(histV.compare("idL1") == 0 || histV.compare("idL2") == 0 || histV.compare("bJ1") == 0 || histV.compare("bJ2") == 0) vlim.push_back(Limits(-0.5, 1.5, 2));
            else if(histV.compare("jmult") == 0) vlim.push_back(Limits(-0.5, 9.5, 10));
            else if(histV.compare("bmult") == 0) vlim.push_back(Limits(-0.5, 9.5, 10));
            else if(histV.compare("mLQmin") == 0 || histV.compare("mLQavg") == 0 || histV.compare("mLQmax") == 0) vlim.push_back(Limits(0, 1500, 150));
            else if(histV.compare("mLQdiff") == 0) vlim.push_back(Limits(0, 1000, 200));
            //else if(histV.compare("mLQavg") == 0) vlim.push_back(Limits(0, 1500, 150));
            else printf("No limits set for variable: %s", histV.c_str());
        }
    }

    TH1 *hist = 0;

    if(vlim.size() == 1)       
    {
        if(!smooth) hist = new TH1D(hname, hname, vlim[0].nb, vlim[0].thll, vlim[0].thul);
        else        hist = new TH1D(hname, hname, 100, vlim[0].thll, vlim[0].thul);
        
    }
    else if (vlim.size() == 2) 
    {
        hist = new TH2D(hname, htitle, vlim[0].nb, vlim[0].thll, vlim[0].thul, vlim[1].nb, vlim[1].thll, vlim[1].thul);
        smooth = false;
    }
    else
    {
        printf("!!!Too many histogram dimmensions!!!\n");
        return 0;
    }

    for(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iT = bgtvec.begin(); iT != bgtvec.end(); ++iT)
    {
        std::vector<double> values;

        // prepair tuple with basic cuts
        if(ll)
        {
            if(iT->first.cutlevel != 17 && (iT->first.cutlevel < ll->cutlevel)) continue;
            if((iT->first.mlljj    < ll->mlljj   )) continue;
            if((iT->first.mll      < ll->mll     )) continue;
            if((iT->first.l1pt     < ll->l1pt    )) continue;
            if((iT->first.l1eta    < ll->l1eta   )) continue;
            if((iT->first.l1phi    < ll->l1phi   )) continue;
            if((iT->first.l2pt     < ll->l2pt    )) continue;
            if((iT->first.l2eta    < ll->l2eta   )) continue;
            if((iT->first.l2phi    < ll->l2phi   )) continue;
            if((iT->first.j1pt     < ll->j1pt    )) continue;
            if((iT->first.j1eta    < ll->j1eta   )) continue;
            if((iT->first.j1phi    < ll->j1phi   )) continue;
            if((iT->first.j2pt     < ll->j2pt    )) continue;
            if((iT->first.j2eta    < ll->j2eta   )) continue;
            if((iT->first.j2phi    < ll->j2phi   )) continue;
        }
        if(ul)
        {
            if(iT->first.cutlevel != 17 && (iT->first.cutlevel > ul->cutlevel)) continue;
            if((iT->first.mlljj    > ul->mlljj   )) continue;
            if((iT->first.mll      > ul->mll     )) continue;
            if((iT->first.l1pt     > ul->l1pt    )) continue;
            if((iT->first.l1eta    > ul->l1eta   )) continue;
            if((iT->first.l1phi    > ul->l1phi   )) continue;
            if((iT->first.l2pt     > ul->l2pt    )) continue;
            if((iT->first.l2eta    > ul->l2eta   )) continue;
            if((iT->first.l2phi    > ul->l2phi   )) continue;
            if((iT->first.j1pt     > ul->j1pt    )) continue;
            if((iT->first.j1eta    > ul->j1eta   )) continue;
            if((iT->first.j1phi    > ul->j1phi   )) continue;
            if((iT->first.j2pt     > ul->j2pt    )) continue;
            if((iT->first.j2eta    > ul->j2eta   )) continue;
            if((iT->first.j2phi    > ul->j2phi   )) continue;
        }
        if(iT->first.weight > 1000 || iT->first.weight < 0.001) continue;  //sanity check on weights

        //if(!runFilter(iT)) continue;

        // dynamical cuts are applied here
        bool passCut = true;
        for(std::map<std::pair<std::string, char>, double>::const_iterator iC = cuts.begin(); iC != cuts.end(); ++iC)
        {
            passCut = passCut && dynamicalCut(getTupleVar(iC->first.first, iT->first), iC->second, iC->first.second);
        }
        if((!passCut && !invertCuts) || (passCut && invertCuts)) continue;

        //int nbjet = 0;
        //if(iT->first.j1B + iT->first.j2B < nbjet) continue;
        //if(iT->first.j1B + iT->first.j2B >= 1) printf("b run: %d\n", iT->first.run);

        // prepair appropriate variables for fill
        for(std::vector<std::string>::const_iterator ihlabel = histQs.begin(); ihlabel != histQs.end(); ++ihlabel)
        {
            values.push_back(getTupleVar(*ihlabel, iT->first));
        }

        // fill histograms with prepaired values and proper weights
        if(values.size() == 1)
        {
            if(histQs.begin()->compare("cutlevel") == 0)
            {
                for(int i = -1; i <= values[0]; i++)
                {
                    hist->Fill(i, iT->first.weight);
                }
            }
            else hist->Fill(values[0], iT->first.weight);
        }
        else if(values.size() == 2) ((TH2*)hist)->Fill(values[0], values[1], iT->first.weight);

        //if(iT->first.cutlevel >= 5 && iT->first.mlljj > 700) printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d:%d:%d\n", iT->first.mlljj, iT->first.mll, getTupleVar("mJJ", iT->first), iT->first.l1pt, iT->first.l2pt, iT->first.j1pt, iT->first.j2pt, iT->first.run, iT->first.ls, iT->first.event);

    }
    
    return hist;
}

double HnuPlots::getTupleVar(std::string var, const HeavyNuTree::HNuSlopeFitInfo& tpls)
{
    double pL1 = tpls.l1pt * cosh(tpls.l1eta);
    double pL2 = tpls.l2pt * cosh(tpls.l2eta);
    double pLmin = std::min(pL1, pL2);
    double pLmax = std::max(pL1, pL2);

    double gR9L1 = tpls.sE1 - tpls.rhE1;
    double gR9L2 = tpls.sE2 - tpls.rhE2;
    double rhoScL1 = tpls.rhE1 / tpls.sE1;
    double rhoScL2 = tpls.rhE2 / tpls.sE2;

    if(var.compare("mWR") == 0)        return (tpls.mlljj/1000.0);
    else if(var.compare("mLL") == 0 || var.compare("mLLZoom") == 0)   return (tpls.mll);
    else if(var.compare("mLLNorm") == 0)  return (tpls.mll/1000.0);
    else if(var.compare("mNuR1") == 0)
    {
        TLorentzVector J1, J2, L;
        J1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        J2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        L.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        return ((J1 + J2 + L).M());
    }
    else if(var.compare("mNuR2") == 0)
    {
        TLorentzVector J1, J2, L;
        J1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        J2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        L.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        return ((J1 + J2 + L).M());
    }
    else if(var.compare("mOuR1") == 0)
    {
        TLorentzVector L1, L2, J;
        L1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        L2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        J.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        return ((L1 + L2 + J).M());
    }
    else if(var.compare("mOuR2") == 0)
    {
        TLorentzVector L1, L2, J;
        L1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        L2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        J.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        return ((L1 + L2 + J).M());
    }
    else if(var.compare("ptL1") == 0)  return (tpls.l1pt );
    else if(var.compare("etaL1") == 0) return (tpls.l1eta);
    else if(var.compare("phiL1") == 0) return (tpls.l1phi);
    else if(var.compare("ptL2") == 0)  return (tpls.l2pt );
    else if(var.compare("etaL2") == 0) return (tpls.l2eta);
    else if(var.compare("phiL2") == 0) return (tpls.l2phi);
    else if(var.compare("ptJ1") == 0)  return (tpls.j1pt );
    else if(var.compare("etaJ1") == 0) return (tpls.j1eta);
    else if(var.compare("phiJ1") == 0) return (tpls.j1phi);
    else if(var.compare("ptJ2") == 0)  return (tpls.j2pt );
    else if(var.compare("etaJ2") == 0) return (tpls.j2eta);
    else if(var.compare("phiJ2") == 0) return (tpls.j2phi);
    else if(var.compare("bJ1") == 0) return (tpls.j1B);
    else if(var.compare("bJ2") == 0) return (tpls.j2B);
    else if(var.compare("pL1") == 0)   return (pLmax);
    else if(var.compare("pL2") == 0)   return (pLmin);
    else if(var.compare("SS") == 0)    return (tpls.cL1 == tpls.cL2);
    else if(var.compare("idL1") == 0)  return (tpls.l1id);
    else if(var.compare("idL2") == 0)  return (tpls.l2id);
    else if(var.compare("mJJ") == 0)
    {
        TLorentzVector J1, J2;
        J1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        J2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        return ((J1 + J2).M() );
    }
    else if(var.compare("n_vertex") == 0) return (tpls.n_primaryVertex);
    else if(var.compare("seL1") == 0)     return (std::max(tpls.sE1, tpls.sE2));
    else if(var.compare("seL2") == 0)     return (std::min(tpls.sE1, tpls.sE2));
    else if(var.compare("rhL1") == 0)     return (std::max(tpls.rhE1, tpls.rhE2));
    else if(var.compare("rhL2") == 0)     return (std::min(tpls.rhE1, tpls.rhE2));
    else if(var.compare("gR9L1") == 0)    return ((pL1 > pL2)?gR9L1:gR9L2);
    else if(var.compare("gR9L2") == 0)    return ((pL1 < pL2)?gR9L1:gR9L2);
    else if(var.compare("rhoScL1") == 0)  return ((pL1 > pL2)?rhoScL1:rhoScL2);
    else if(var.compare("rhoScL2") == 0)  return ((pL1 < pL2)?rhoScL1:rhoScL2);
    else if(var.compare("met") == 0)      return (tpls.met);
    else if(var.compare("dEtaL") == 0)    return (tpls.l1eta - tpls.l2eta);
    else if(var.compare("dEtaJ") == 0)    return (tpls.j1eta - tpls.j2eta);
    else if(var.compare("dPhiL") == 0)
    {
        TLorentzVector L1, L2;
        L1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        L2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        return (ROOT::Math::VectorUtil::DeltaPhi(L1, L2));
    }
    else if(var.compare("dPhiJ") == 0)
    {
        TLorentzVector J1, J2;
        J1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        J2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        return (ROOT::Math::VectorUtil::DeltaPhi(J1, J2));
    }
    else if(var.compare("run") == 0)      return ((double)tpls.run);
    else if(var.compare("cutlevel") == 0) return ((double)tpls.cutlevel);
    else if(var.compare("jmult") == 0) return ((double)tpls.jmult);
    else if(var.compare("bmult") == 0) return ((double)tpls.bmult);
    else if(var.compare("st") == 0) return (tpls.j1pt + tpls.j2pt + tpls.l1pt + tpls.l2pt);
    else if(var.compare("mLQmin") == 0)
    {
        TLorentzVector j1, j2, l1, l2;
        j1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        j2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        l1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        l2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        double m11 = (j1 + l1).M();
        double m21 = (j2 + l1).M();
        double m12 = (j1 + l2).M();
        double m22 = (j2 + l2).M();
        
        if(fabs(m11 - m22) < fabs(m12 - m21))
        {
            return min(m11, m22);
        }
        else return min(m12, m21);
    }
    else if(var.compare("mLQmax") == 0)
    {
        TLorentzVector j1, j2, l1, l2;
        j1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        j2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        l1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        l2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        double m11 = (j1 + l1).M();
        double m21 = (j2 + l1).M();
        double m12 = (j1 + l2).M();
        double m22 = (j2 + l2).M();
        
        if(fabs(m11 - m22) < fabs(m12 - m21))
        {
            return max(m11, m22);
        }
        else return max(m12, m21);
    }
    else if(var.compare("mLQavg") == 0)
    {
        TLorentzVector j1, j2, l1, l2;
        j1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        j2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        l1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        l2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        double m11 = (j1 + l1).M();
        double m21 = (j2 + l1).M();
        double m12 = (j1 + l2).M();
        double m22 = (j2 + l2).M();
        
        if(fabs(m11 - m22) < fabs(m12 - m21))
        {
            return (m11 + m22) / 2;
        }
        else return (m12 + m21) / 2;
    }
    else if(var.compare("mLQdiff") == 0)
    {
        TLorentzVector j1, j2, l1, l2;
        j1.SetPtEtaPhiM(tpls.j1pt, tpls.j1eta, tpls.j1phi, 0);
        j2.SetPtEtaPhiM(tpls.j2pt, tpls.j2eta, tpls.j2phi, 0);
        l1.SetPtEtaPhiM(tpls.l1pt, tpls.l1eta, tpls.l1phi, 0);
        l2.SetPtEtaPhiM(tpls.l2pt, tpls.l2eta, tpls.l2phi, 0);
        double m11 = (j1 + l1).M();
        double m21 = (j2 + l1).M();
        double m12 = (j1 + l2).M();
        double m22 = (j2 + l2).M();
        
        return min(fabs(m11 - m22), fabs(m12 - m21));
    }
    else printf("Variable not found: %s, returning -999.0\n", var.c_str());

    return -999.0;
}

bool HnuPlots::dynamicalCut(double var, double cut, char cutType)
{
    if     (cutType == '<') return var < cut;
    else if(cutType == '>') return var > cut;
    else printf("Unrecognized cut type, %c\n", cutType);

    return false;
}

bool HnuPlots::runFilter(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iE)
{
    char pkey[128];
    sprintf(pkey,"%d:%d:%d", iE->first.run, iE->first.ls, iE->first.event);
    std::string key(pkey);

    std::set<std::string> events;

    events.insert("199812:652:755626716");
    events.insert("198210:40:27320573");
    events.insert("207905:1341:1662736109");
    events.insert("202060:238:287005652");
    events.insert("205666:384:601051081");
    events.insert("207372:52:85843102");
    events.insert("194912:1518:1945772859");
    events.insert("199409:315:416621165");
    events.insert("208427:552:824605654");
    events.insert("195915:556:836688041");
    events.insert("206210:286:264896604");
    events.insert("207273:106:124694968");

    return (events.find(key) != events.end());
}

HnuPlots::FileStruct::FileStruct()
{
    intLumi = cs = kfactor = 1.0;
    clow = chigh = 0.0;
    normbin = 1;
    file = "";
    loadtuple = false;
    histFromTuple = false;
}

HnuPlots::FileStruct::FileStruct(std::string l, std::string f, std::string h)
{
    intLumi = cs = kfactor = 1.0;
    clow = chigh = 0.0;
    normbin = 1;
    file = "";
    loadtuple = false;
    histFromTuple = false;
    label = l;
    file = f;
    histpath = h;
}

HnuPlots::FileStruct::FileStruct(std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh, double cl, double ch, bool px, int nb, bool un, double nll, double nul, bool lt, bool lhft, double ll, double ul, double bn, HeavyNuTree::HNuSlopeFitInfo* tll, HeavyNuTree::HNuSlopeFitInfo * tul, bool smooth)
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
    normbin = nb;
    unh = un;
    normll = nll;
    normul = nul;
    loadtuple = lt;
    histFromTuple = lhft;
    thll = ll;
    thul = ul;
    thb = bn;
    tpll = tll;
    tpul = tul;
    smooth_hist = smooth;
}

HnuPlots::FileStruct::FileStruct(bool lt, std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh, int cl, int nb)
{

    loadtuple = lt;
    label = l;
    file = f;
    histpath = h;
    intLumi = iL;
    cs = c;
    kfactor = kf;
    normhistpath = nh;
    cutlevel = cl;
    normbin = nb;
}

HnuPlots::FileStruct::~FileStruct(){ }

void HnuPlots::plot()
{
    if(datahist.hist->IsA()->InheritsFrom(TH2::Class())) plot2D();
    else                                                 plot1D();
}

void HnuPlots::plot1D()
{
    using namespace std;

    //gROOT->SetStyle("Plain");
    setTDRStyle();

    if(rebin > 1)
    {
        datahist.hist->Rebin(rebin);
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist->Rebin(rebin);
        }
        for(vector<HnuPlots::HistStruct>::const_iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
        {
            if(!ihsig->smooth_hist) ihsig->hist->Rebin(rebin);
        }
    }
    else if(rebin < 0)
    {
        datahist.hist = datahist.hist->Rebin(sizeof(bins) / sizeof(double) - 1, "mWR_limitbins" , bins);
        for(int i = 1; i <= datahist.hist->GetNbinsX(); i++)
        {
            datahist.hist->SetBinContent(i, 0.2 * datahist.hist->GetBinContent(i) / datahist.hist->GetBinWidth(i));
            datahist.hist->SetBinError(i, 0.2 * datahist.hist->GetBinError(i) / datahist.hist->GetBinWidth(i));
        }
        for(vector<HnuPlots::HistStruct>::iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist = ihbg->hist->Rebin(sizeof(bins) / sizeof(double) - 1, "blarg" , bins);
            for(int i = 1; i <= ihbg->hist->GetNbinsX(); i++)
            {
                ihbg->hist->SetBinContent(i, 0.2 * ihbg->hist->GetBinContent(i) / ihbg->hist->GetBinWidth(i));
            }
        }
        for(vector<HnuPlots::HistStruct>::iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
        {
            if(!ihsig->smooth_hist)
            {
                ihsig->hist = ihsig->hist->Rebin(sizeof(bins) / sizeof(double) - 1, "splat", bins);
                for(int i = 1; i <= ihsig->hist->GetNbinsX(); i++)
                {
                    ihsig->hist->SetBinContent(i, 0.2 * ihsig->hist->GetBinContent(i) / ihsig->hist->GetBinWidth(i));
                }
            }
        }
    }

    //BLAHBLAHBLAH
    //TF1 *tf = new TF1("tf","expo", 800, 10000);
    //tf->SetLineColor(kBlack);
    //tf->SetLineWidth(2);
    //tf->SetLineStyle(2);
    //datahist.hist->Fit(tf, "LQN", "", 800, 10000);

    char lumistamp[128];
    //sprintf(lumistamp, "%.1f fb^{-1} at 8 TeV", iLumi / 1000);
    //sprintf(lumistamp, "CMS    #sqrt{s} = 8 TeV    L = %0.1f fb^{-1}", iLumi / 1000);
    sprintf(lumistamp, "%0.1f fb^{-1} (8 TeV)", iLumi / 1000);

    if(autosort) sort(bghists.begin(), bghists.end(), compHistInt);

    bool isGeV = true;
    if(!yaxislabel.compare("please auto set the axis"))
    {
        char temp[128];
        if(xaxislabel.find("GeV") < xaxislabel.size())
        {
            if(rebin >= 0) sprintf(temp, "Events / %.0f GeV", datahist.hist->GetBinWidth(1));
            else 
            {
                sprintf(temp, "Events / 200 GeV");
            }
            yaxislabel = temp;
        }
        else if(xaxislabel.find("TeV") < xaxislabel.size())
        {
            if(rebin >= 0) 
            {
                if(datahist.hist->GetBinWidth(1) >= 0.1)  sprintf(temp, "Events / %.1f TeV", datahist.hist->GetBinWidth(1));
                else                                      sprintf(temp, "Events / %.2f TeV", datahist.hist->GetBinWidth(1));
            }
            else 
            {
                sprintf(temp, "Events / 0.2 TeV");
            }
            yaxislabel = temp;
            isGeV = false;
        }
        else
        {
            sprintf(temp, "Events / %.2f", datahist.hist->GetBinWidth(1));
            yaxislabel = temp;
        }
    }

    TCanvas *c1;
    double fontScale = 1.0;
    if(plotSMoData)
    {
        c1 = new TCanvas("c1", "c1", 800, 900);
        c1->Divide(1, 2);
        c1->cd(1);
        gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
        gPad->SetBottomMargin(0.01);
        fontScale = 1.0;
    }
    else
    {
        c1 = new TCanvas("c1", "c1", 800, 800);
        c1->Divide(1, 1);
        c1->cd(1);
        gPad->SetPad("p1", "p1", 0, 0, 1, 1, kWhite, 0, 0);
        gPad->SetBottomMargin(0.15);
        //fontScale = 8.0 / 9;
        fontScale = 6.5 / 8;
    }
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.06 * (8.0 / 6.5) * fontScale);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);
    
    datahist.hist->SetMarkerColor(kBlack);
    datahist.hist->SetMarkerStyle(20);
    datahist.hist->SetLineWidth(2.0);

    //TLegend *leg = new TLegend(0.52, 0.67, 0.94, 0.91);
    TLegend *leg = new TLegend(0.45, 0.61, 0.89, 0.91);
    leg->SetFillStyle(0); //Color(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    float dataintegral = 0.0;
    if(rebin >= 0) dataintegral = datahist.hist->Integral(0, datahist.hist->GetNbinsX() + 1);
    else dataintegral = datahist.hist->Integral(1, datahist.hist->GetNbinsX(), "width")/(isGeV?200:0.2);
    char datahllabel[128];
    sprintf(datahllabel, "%s (%.0f)", datahist.label.c_str(), dataintegral);
    leg->AddEntry(datahist.hist, datahllabel, "ep");

    THStack *hbg = new THStack("Background", "background");
    //TH1* sig = 0;
    //if(sighists.size() > 0) sig = (TH1*)(sighists[0].hist->Clone());
    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.end() - 1; ihbg != bghists.begin() - 1; ihbg--)
    {
        hbg->Add(ihbg->hist);
        //    if(sig) sig->Add(ihbg->hist);
    }
    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        float integral = 0.0;
        if(rebin >= 0) integral = ihbg->hist->Integral(0, ihbg->hist->GetNbinsX() + 1);
        else integral = ihbg->hist->Integral(1, ihbg->hist->GetNbinsX(), "width")/(isGeV?200:0.2);
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)", ihbg->label.c_str(), floor(integral + 0.5));
        leg->AddEntry(ihbg->hist, hllabel, "F");
    }
    //double sigMaxMin = datahist.hist->GetMaximum();
    for(vector<HnuPlots::HistStruct>::const_iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
    {
        if(!ihsig->smooth_hist)
        {
            float integral = 0.0;
            if(rebin >= 0) integral = ihsig->hist->Integral(0, ihsig->hist->GetNbinsX() + 1);
            else integral = ihsig->hist->Integral(1, ihsig->hist->GetNbinsX(), "width") / (isGeV?200:0.2);
            char hllabel[128];
            sprintf(hllabel, "%s (%.0f)", ihsig->label.c_str(), floor(integral + 0.5));
            leg->AddEntry(ihsig->hist, hllabel, "L");
        }
        //sigMaxMin = std::min(ihsig->hist->GetMaximum(), sigMaxMin);
    }

    //BLAHBLAHBLAH
    //leg->AddEntry(tf, "Exponential Fit");

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
    if(xmin != xmax) dummy->GetXaxis()->SetRangeUser(xmin, xmax);
    //dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    if(islog)
    {
        if(ymin == ymax) dummy->GetYaxis()->SetRangeUser(std::max(0.0001, 0.2 * std::min(hbg->GetMaximum(), 0.3 * datahist.hist->GetMinimum(0.0001))), std::max(hbg->GetMaximum(), datahist.hist->GetMaximum())*4);
        else             dummy->GetYaxis()->SetRangeUser(ymin, ymax);
        gPad->SetLogy(1);
    }
    else
    {
        dummy->GetYaxis()->SetRangeUser(0.001, std::max(hbg->GetMaximum(), datahist.hist->GetMaximum())*1.2);
    }
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.1*1.05 / (fontScale));
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->SetStats(0);
    if(plotSMoData) dummy->GetXaxis()->SetTitle(0);
    else            dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    fixOverlay();
    dummy->Draw();
    fixOverlay();
    //if(sig) sig->Draw("hist same");
    //hbg->Draw("hist same");
    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        TH1 *hist = (TH1*)ihbg->hist->Clone("bob"), *hist2 = (TH1*)ihbg->hist->Clone("bob2");
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg2 = ihbg + 1; ihbg2 != bghists.end(); ihbg2++)
        {
            hist->Add(ihbg2->hist);
            hist2->Add(ihbg2->hist);
        }
        hist->SetFillStyle(1000);
        hist->SetFillColor(10); //this is non-transperent white, apparently kWhite is transparent!!!!!!!!
        hist->Draw("same hist");
        hist2->Draw("same hist");
    }
    fixOverlay();
    for(std::vector<HnuPlots::HistStruct >::const_iterator isig = sighists.begin(); isig != sighists.end(); isig++)
    {
        if(!isig->smooth_hist) isig->hist->Draw("hist same");
        else
        {
            TGraph *ttg = new TGraph();
            
            for(int i = 1; i <= isig->hist->GetNbinsX(); ++i)
            {
                ttg->SetPoint(i - 1, isig->hist->GetBinCenter(i), isig->hist->GetBinContent(i));
            }
            
            TGraph *ttg1 = (TGraph*)ttg->Clone();
            TGraph *ttg2 = (TGraph*)ttg->Clone();
            TGraph *ttg4 = (TGraph*)ttg->Clone();

            TGraphSmooth *gs = new TGraphSmooth("normal1");
            TGraph *tstg = gs->SmoothSuper(ttg1, "", 10, 0.001);
            TGraphSmooth *gs2 = new TGraphSmooth("normal2");
            TGraph *tstg2 = gs2->SmoothSuper(ttg2, "", 10, 0.12);
            TGraphSmooth *gs4 = new TGraphSmooth("normal3");
            TGraph *tstg3 = gs4->SmoothSuper(ttg4, "", 10, 0.10);
            
            //tstg->Draw("same L");
            //tstg2->Draw("same L");
            
            TF1 *tf1 = new TF1("tf1", "expo", 1.0, 2.0);
            TF1 *tf2 = new TF1("tf2", "expo", 2.8, 3.0);
            
            isig->hist->Fit(tf1, "RLNQ");
            isig->hist->Fit(tf2, "RLNQ");
            
            double scale = isig->hist->Integral(1, isig->hist->GetNbinsX(), "width") / (tstg->Integral()*isig->hist->GetBinWidth(1))*datahist.hist->GetBinWidth(1);
            for(int i = 0; i <= tstg->GetN(); ++i)
            {
                double x, y1, y2, y3;
                tstg->GetPoint(i, x, y1);
                tstg2->GetPoint(i, x, y2);
                tstg3->GetPoint(i, x, y3);
                
                std::cout << x << std::endl;
                
                //if(x < 2200)      ttg->SetPoint(i, x, y2 * scale * tstg->Eval(2190) / tstg2->Eval(2190));
                //else if(x > 2700) ttg->SetPoint(i, x, y3 * scale * tstg->Eval(2710) / tstg3->Eval(2710));
                if(x < 2.0)      ttg->SetPoint(i, x, tf1->Eval(x) * scale * tstg->Eval(1.990) / tf1->Eval(1.990));
                else if(x > 2.8) ttg->SetPoint(i, x, tf2->Eval(x) * scale * tstg->Eval(2.810) / tf2->Eval(2.810));
                else              ttg->SetPoint(i, x, y1 * scale);
            }
            
            ttg->SetLineWidth(2);
            ttg->SetLineStyle(2);
            ttg->SetLineColor(kRed + 2);
            
            leg->AddEntry(ttg, isig->label.c_str(), "L");
            
            ttg->Draw("same L");
        }
    }
    datahist.hist->Draw("same pe");
    fixOverlay();
    leg->Draw();
    
    TLatex mark;
    mark.SetTextSize(0.04 * 1.1 * 8 / 6.5 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    //mark.DrawLatex(0.17, 0.95, "CMS");
    //mark.DrawLatex(0.68, 0.95, lumistamp);
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.017, lumistamp);
    //mark.SetTextSize((0.04 * 7 / 6.5 * fontScale)*1.25);
    bool isCMS = false;
    if(isCMS)
    {
        mark.SetTextAlign(13);
        mark.SetTextSize(0.04 * 1.1 * 8 / 6.5 * 1.25 * fontScale);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.027, 1 - (gPad->GetTopMargin() + 0.027), "CMS"); // #scale[0.8]{#it{Preliminary}}");
    }
    else
    {
        mark.SetTextAlign(11);
        mark.SetTextSize(0.04 * 1.1 * 8 / 6.5 * 1.25 * fontScale);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.04 * 1.1 * 8 / 6.5 * fontScale);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.12, 1 - (gPad->GetTopMargin() - 0.017), "unpublished");
    }
    
    //mark.SetTextSize(0.04 * 7 / 6.5 * fontScale);
    //mark.DrawLatex(gPad->GetLeftMargin() + 0.025, 1 - (gPad->GetTopMargin() + 0.065), "#it{Preliminary}");
    fixOverlay();
    
    //BLAHBLAHBLAH
    //tf->Draw("same");

    if(plotSMoData)
    {
        bool isRatio = true;
        c1->cd(2);
        gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.37);
        TH1* chbg = 0;
        for(vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            if(ihbg == bghists.begin()) chbg = (TH1*)ihbg->hist->Clone();
            else chbg->Add(ihbg->hist);
        }
        TH1* chdata = (TH1*)datahist.hist->Clone();
        //TH1* chsig = (sig)?(TH1*)sig->Clone():0;
        //if(chsig) chsig->Divide(chbg);
        double ratio = 0, error = 0;
        for(int i = 1; i <= chdata->GetNbinsX(); i++)
        {
            if(isRatio)
            {
                if(chdata->GetBinContent(i) != 0 && chbg->GetBinContent(i) != 0)
                {
                    ratio = chdata->GetBinContent(i) / chbg->GetBinContent(i);
                    error = ratio * chdata->GetBinError(i) / chdata->GetBinContent(i);
                    chdata->SetBinContent(i, ratio);
                    chdata->SetBinError(i, error);
                }
                else
                {
                    chdata->SetBinContent(i, -1);
                    chdata->SetBinError(i, 0);
                }
            }
            else
            {
                if(chdata->GetBinError(i) > 1e-10) chdata->SetBinContent(i, (chdata->GetBinContent(i) - chbg->GetBinContent(i))/chdata->GetBinError(i));
                else                               chdata->SetBinContent(i, -100);
                chdata->SetBinError(i, 0);
            }
        }

        
        //double d2ymax = max(2.2, min(4.5, chdata->GetMaximum(25.0)*1.2));
        double d2ymin = max(0.4, chdata->GetMinimum(0) - 0.4);
        std::cout << "2dymin: " << d2ymin << std::endl;
        double d2ymax = min(4.5, max(1.7, chdata->GetMaximum(25.0)*1.3));
        TH1 *dummy2 = new TH1F("dummy2", "dummy2", 1000, datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
        dummy2->GetXaxis()->SetTitle(xaxislabel.c_str());
        dummy2->GetXaxis()->SetTitleOffset(1.05);
        if(isRatio)
        {
            dummy2->GetYaxis()->SetRangeUser(d2ymin, d2ymax);
            dummy2->GetYaxis()->SetTitle("Data/SM");
            dummy2->GetYaxis()->SetTitleOffset(0.42);
            dummy2->GetYaxis()->SetNdivisions(3, 5, 0, true);
            //dummy2->GetYaxis()->SetMoreLogLabels(true);
            //dummy2->GetYaxis()->SetNoExponent(true);
            //gPad->SetLogy(true);
        }
        else
        {
            dummy2->GetYaxis()->SetRangeUser(-3, 3);
            dummy2->GetYaxis()->SetTitle("#frac{Data - SM}{#sigma_{Data}}");
            dummy2->GetYaxis()->SetTitleOffset(0.42);
            dummy2->GetYaxis()->SetNdivisions(3, 5, 0);
        }
        dummy2->GetYaxis()->SetTitleSize(0.16 * 2 / 2.5);
        dummy2->GetYaxis()->SetLabelSize(0.20 * 2 / 2.5);
        dummy2->GetXaxis()->SetTitleSize(0.20 * 2 / 2.5);
        dummy2->GetXaxis()->SetLabelSize(0.20 * 2 / 2.5);
        if(xmin != xmax) dummy2->GetXaxis()->SetRangeUser(xmin, xmax);
        dummy2->SetStats(0);
        dummy2->SetTitle(0);
        if(dummy2->GetNdivisions() % 100 > 5) dummy2->GetXaxis()->SetNdivisions(6, 5, 0);

        TF1 * fline = new TF1("line", "pol0", datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
        fline->SetParameter(0, 1);
        if(isRatio) fline->SetLineColor(kBlack);
        else        fline->SetLineColor(kRed);

        dummy2->Draw();
        
        //Horrible manual axis replacement campaign begin

        //TF1 *f_h2_log10_x_axis = new TF1("f_h2_log10_y_axis", // name
        //                                 //"log10(x)", // formula
        //                                 "log10(x)", // formula
        //                                 d2ymin, // xmin
        //                                 d2ymax); // xmax
//
        //f_h2_log10_x_axis->SetLineColor(kBlue);
//
        //TGaxis *a = new TGaxis(dummy2->GetXaxis()->GetXmin(), // xmin
        //                       d2ymin, // ymin
        //                       dummy2->GetXaxis()->GetXmin(), // xmax
        //                       d2ymax, // ymax
        //                       "f_h2_log10_y_axis", // funcname
        //                       503, // ndiv (try 100006 or 506, don't try 1006)
        //                       "BS", // chopt (try "BS" or "UBS")
        //                       0.0); // gridlength
//
        //// a->SetTickSize(h2->GetTickLength("X")); // use "the same" size
        //a->SetTickSize(1.5 * dummy2->GetTickLength("Y")); // make it bigger
        //dummy2->SetTickLength(0.0, "Y"); // get rid of "original" ticks
//
        //if (!(TString(a->GetOption())).Contains("U"))
        //{
        //    a->SetLabelFont(dummy2->GetLabelFont("Y")); // use "the same" font
        //    a->SetLabelSize(dummy2->GetLabelSize("Y")); // use "the same" size
        //    dummy2->SetLabelSize(0.0, "Y"); // get rid of "original" labels
        //}
        //
        //TGaxis *a2 = (TGaxis*)a->Clone();
        //a2->SetX1(dummy2->GetXaxis()->GetXmax());
        //a2->SetX2(dummy2->GetXaxis()->GetXmax());
        //a2->SetOption("+ UBS");
        //
        ////Most horrible and terrible thing I have done to get a plot to look good
        //TText * axislabel = new TText();
        //axislabel->SetNDC(false);
        //axislabel->SetTextAlign(32);
        //axislabel->SetTextFont(a->GetLabelFont());
        //axislabel->SetTextSize(a->GetLabelSize());
        
        //Horrible manual axis replacement campaign end

        //if(chsig) chsig->Draw("hist same");

        //BLAHBLAHBLAH
        //TH1 *tfrh = (TH1*)datahist.hist->Clone("tfrh");
        //tfrh->SetLineStyle(2);
        //for(int ibin = 1; ibin < tfrh->GetNbinsX() && ibin < chbg->GetNbinsX(); ibin++)
        //{
        //    if(ibin < 4) tfrh->SetBinContent(ibin, 1);
        //    else if(chbg->GetBinContent(ibin) > 1e-15) tfrh->SetBinContent(ibin, tf->Eval(tfrh->GetBinCenter(ibin))/chbg->GetBinContent(ibin));
        //    else tfrh->SetBinContent(ibin, 1);
        //}

        //tfrh->Draw("same hist");

        if(true)//isRatio)
        {
            TExec * setex2 = new TExec("setex2", "gStyle->SetErrorX(0.5)");
            setex2->Draw();
            
            TH1* sysdata = (TH1*)datahist.hist->Clone();
            sysdata = sysdata->Rebin(sizeof(systBins) / sizeof(double) - 1, "systdata", systBins);
            TH1* sysbg = (TH1*)chbg->Clone();
            sysbg = sysbg->Rebin(sizeof(systBins) / sizeof(double) - 1, "systbg", systBins);

            TH1 **tgs = new TH1*[systematics.size()];
            int itg = 0;
            const int ebColors[] = {kRed - 9, kBlue - 1}, NEBCOLORS = sizeof(ebColors) / sizeof(int);
            //const int ebStyles[] = {3001, 3001}, NEBSTYLES = sizeof(ebColors) / sizeof(int);
            double chi2 = 0.0;
            double chi2_2 = 0.0;
            for(std::vector< std::vector<float> >::const_iterator sit = systematics.begin(); sit != systematics.end(); ++sit)
            {
                char hname[128];
                sprintf(hname, "hsyst_%d", itg);
                tgs[itg] = new TH1F(hname, hname, sizeof(systBins) / sizeof(double) - 1, systBins);
                for(int i = 1; i <= tgs[itg]->GetNbinsX(); i++)
                {
                    if(isRatio)
                    {
                        tgs[itg]->SetBinContent(i, 1.0);
                        if(i - 1 < (int)sit->size()) tgs[itg]->SetBinError(i, fabs(sit->at(i - 1)));
                        else if(sit->size() > 0) tgs[itg]->SetBinError(i, fabs(sit->back()));
                    }
                    else
                    {
                        tgs[itg]->SetBinContent(i, 0.0);
                        int ibin = sysdata->FindBin(tgs[itg]->GetBinCenter(i));
                        if(i - 1 < (int)sit->size())    
                        {
                            tgs[itg]->SetBinError(i, sysbg->GetBinContent(ibin)*(fabs(sit->at(i - 1)))/sysdata->GetBinError(ibin));
                            std::cout << ibin << "\t" << tgs[itg]->GetBinCenter(i) << "\t" << sysbg->GetBinContent(ibin) << "\t" << (fabs(sit->at(i - 1))) << "\t" << sysdata->GetBinError(ibin) << "\t" << sysdata->GetBinContent(ibin) << std::endl;
                        }
                        else if(sit->size() > 0)     tgs[itg]->SetBinError(i, sysbg->GetBinContent(ibin)*(fabs(sit->back()))/sysdata->GetBinError(ibin));
                    }

                    if(sit == systematics.begin())
                    {
                        double nbg = 0.0;
                        for(std::vector<HistStruct>::const_iterator bghr = bghists.begin(); bghr != bghists.end(); ++bghr)
                        {
                            nbg += bghr->hist->GetBinContent(i);
                        }
                        double syst_error = nbg * fabs(sit->at(i - 1));
                        double stat_error = datahist.hist->GetBinError(i);
                        chi2 += pow(nbg - datahist.hist->GetBinContent(i), 2)/(syst_error*syst_error + stat_error*stat_error);
                        chi2_2 += pow(nbg - datahist.hist->GetBinContent(i), 2)/(stat_error*stat_error);
                        //std::cout << nbg << "\t" <<
                    }
                }
                tgs[itg]->SetFillColor(ebColors[itg % NEBCOLORS]);
                //tgs[itg]->SetFillStyle(ebStyles[itg % NEBSTYLES]);
                tgs[itg]->SetMarkerStyle(0);
                tgs[itg]->Draw("E2 L same");
                itg++;
            }
            printf("Chi^2 = %f\n", chi2);
            printf("Chi^2 (stat only) = %f\n", chi2_2);

            if(systematics.size())
            {
                //TLine *sysStartLine = new TLine(600.0, std::max(-0.9, std::min(0.0, 1.0 - 1.2 * tgs[0]->GetBinError(tgs[0]->FindBin(3500)))), 600.0, 600.0);
                TLine *sysStartLine = 0;
                if(isRatio) sysStartLine = new TLine(0.6, d2ymin, 0.6, d2ymax);
                else        sysStartLine = new TLine(0.6, -3, 0.6, 3);
                sysStartLine->SetLineColor(kBlack);
                sysStartLine->SetLineStyle(2);
                sysStartLine->Draw();
            }


            TExec *setex = new TExec("setex", "gStyle->SetErrorX(0.0)");
            setex->Draw();
        }

        fline->Draw("same");
        if(!isRatio)
        {
            TF1 * fline2 = new TF1("line", "pol0", datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
            fline2->SetParameter(0, -1);
            fline2->SetLineColor(kRed);
            fline2->Draw("same");
        }
        
        //Draw the terrible axis hack
        //a->Draw();
        //a2->Draw();
        //gPad->Modified();
        //gPad->Update(); // make sure it's redrawn
        //axislabel->DrawText(-a->GetLabelOffset()*(gPad->GetX2() - gPad->GetX1()), 0.5, "0.5");

        chdata->Draw("same P");

        TLine *tl = new TLine();
        tl->SetLineColor(kBlack);
        tl->SetLineWidth(2);
        //if(isRatio)
        //{
        //    for(int i = 1; i <= chdata->GetNbinsX(); i++)
        //    {
        //        if(chdata->GetBinCenter(i) < xmin || chdata->GetBinCenter(i) > xmax) continue;
        //        if(chbg->GetBinContent(i) > 0.0001 && (d2ymax > chdata->GetBinContent(i) + chdata->GetBinError(i)))   tl->DrawLine(chdata->GetBinCenter(i), std::max(0.0, std::min(d2ymax, chdata->GetBinContent(i) + chdata->GetBinError(i))), chdata->GetBinCenter(i), std::max(0.0, chdata->GetBinContent(i) - chdata->GetBinError(i)));
        //        else if((chbg->GetBinContent(i) < 0.0001) && (datahist.hist->GetBinContent(i) > 0)) tl->DrawLine(chdata->GetBinCenter(i), 0.0, chdata->GetBinCenter(i), d2ymax);
        //    }
        //}
        fixOverlay();
    }

    if(saveplots)
    {
        char ofn[128], ofn2[128];
        //int cutlevel = 11111;
        //char * pos = strstr(datahist.hist->GetTitle(), "cut:");
        //int scans = pos?sscanf(pos, "cut:%d", &cutlevel):0;
        //if(scans == 0) cutlevel = 111111;
        //sprintf(tmp, "cut:%da", cutlevel);
        string tmp2 = datahist.hist->GetName();
        if(tmp2.find(';') < tmp2.size()) tmp2.erase(tmp2.find(';'));
        if(tmp2.find("_") < tmp2.size()) tmp2.erase(tmp2.find("_"));
        sprintf(ofn, "%s_%s%s.pdf", formlabel.c_str(), tmp2.c_str(), islog?"":"_linear");
        c1->Print(ofn);
        sprintf(ofn2, "%s_%s%s.png", formlabel.c_str(), tmp2.c_str(), islog?"":"_linear");
        c1->Print(ofn2);
    }
}

void HnuPlots::plot2D()
{
    setTDRStyle();

    TH2 *dummy = new TH2F("dummy", "dummy", 1000, datahist.hist->GetXaxis()->GetBinLowEdge(1), datahist.hist->GetXaxis()->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetXaxis()->GetBinWidth(datahist.hist->GetNbinsX()), 1000, datahist.hist->GetYaxis()->GetBinLowEdge(1), datahist.hist->GetYaxis()->GetBinLowEdge(datahist.hist->GetNbinsY()) + datahist.hist->GetYaxis()->GetBinWidth(datahist.hist->GetNbinsY()));
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.60);
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->SetStats(0);
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetXaxis()->SetTitleSize(0.042);
    dummy->GetXaxis()->SetLabelSize(0.042);
    dummy->GetYaxis()->SetTitleSize(0.042);
    dummy->GetYaxis()->SetLabelSize(0.042);
    if(dummy->GetXaxis()->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);
    if(dummy->GetYaxis()->GetNdivisions() % 100 > 5) dummy->GetYaxis()->SetNdivisions(6, 5, 0);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);

    dummy->Draw();

    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 8 TeV", iLumi / 1000);
    TLatex mark;
    mark.SetTextSize(0.042);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    mark.DrawLatex(0.15, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.68, 0.95, lumistamp);

    ((TH2*)datahist.hist)->Draw("same");
}

void HnuPlots::plotMCComp(bool rescale)
{
    using namespace std;

    //gROOT->SetStyle("Plain");
    setTDRStyle();

    if(rebin > 1)
    {
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist->Rebin(rebin);
        }
    }
    else if(rebin < 0)
    {
        for(vector<HnuPlots::HistStruct>::iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist = ihbg->hist->Rebin(sizeof(bins) / sizeof(double) - 1, "blarg" , bins);
            for(int i = 1; i <= ihbg->hist->GetNbinsX(); i++)
            {
                ihbg->hist->SetBinContent(i, ihbg->hist->GetBinContent(i) / ihbg->hist->GetBinWidth(i));
            }
        }
    }

    if(rescale)
    {
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin() + 1; ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist->Scale(bghists.begin()->hist->Integral() / ihbg->hist->Integral());
        }
    }

    //sort(bghists.begin(), bghists.end(), compHistInt);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    c1->SetLogy(islog);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);

    TLegend *leg = new TLegend(0.50, 0.75, 0.94, 0.94);
    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.end() - 1; ihbg != bghists.begin() - 1; ihbg--)
    {
        leg->AddEntry(ihbg->hist, ihbg->label.c_str());
    }
    //leg->AddEntry(datahist.second, datahist.first.c_str());

    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, bghists[0].hist->GetBinLowEdge(1), bghists[0].hist->GetBinLowEdge(bghists[0].hist->GetNbinsX()) + bghists[0].hist->GetBinWidth(bghists[0].hist->GetNbinsX()));
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetYaxis()->SetRangeUser(0.1, std::max(bghists[0].hist->GetMaximum(), bghists[1].hist->GetMaximum())*1.2);
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.0);
    dummy->SetStats(0);
    dummy->SetTitle(0);
    dummy->GetXaxis()->SetNdivisions(6, 5, 0);


    fixOverlay();
    dummy->Draw();
    const int compcolors[] = {kRed, kBlue, kBlack, kGreen + 2};
    int i = 0;
    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        ihbg->hist->SetLineColor(compcolors[i % 3]);
        ihbg->hist->SetFillColor(0);
        ihbg->hist->SetMarkerStyle(21 + i);
        ihbg->hist->SetMarkerColor(compcolors[i % 3]);
        ihbg->hist->Draw("hist same");
        fixOverlay();
        i++;
    }
    leg->Draw("same");

    if(bghists.size() >= 2)
    {
        for(int ibin = bghists[0].hist->FindBin(600); ibin <= bghists[0].hist->GetNbinsX(); ibin++)
        {
            printf("%f,", bghists[0].hist->GetBinContent(ibin) / bghists[1].hist->GetBinContent(ibin));
        }
        printf("\n");
        for(int ibin = bghists[0].hist->FindBin(600); ibin <= bghists[0].hist->GetNbinsX(); ibin++)
        {
            double b1 = bghists[0].hist->GetBinContent(ibin);
            double b2 = bghists[1].hist->GetBinContent(ibin);
            double e1 = bghists[0].hist->GetBinError(ibin);
            double e2 = bghists[1].hist->GetBinError(ibin);

            printf("%f,", (b1 / b2) * sqrt(e1*e1/(b1*b1) + e2*e2/(b2*b2)));
        }
        printf("\n");
    }
}

void HnuPlots::plotMCShape(std::string bgfilename)
{
    using namespace std;

    setTDRStyle();

    double fitmin = 800.0;

    bool drawtrialfuncs = true;

    FILE *fbackfits = fopen(bgfilename.c_str(), "w");

    vector<pair<string, string> > fffs;
    fffs.push_back(make_pair("#it{e^{a+bM}}", "exp([0]+[1]*x)"));
    //fffs.push_back(make_pair("#it{e^{a+bMlog(M)}}", "exp([0]+[1]*x*log(x))"));
    //fffs.push_back(make_pair("#it{e^{a+bM+cM^{2}}}", "exp([0]+[1]*x+[2]*x*x)"));
    //fffs.push_back(make_pair("#it{e^{a+bM+cM^{3}}}", "exp([0]+[1]*x+[2]*x*x*x)"));
    //fffs.push_back(make_pair("#it{e^{a+bM}/M^{c}}", "exp([0]+[1]*x)/pow(x,[2])"));
    //fffs.push_back(make_pair("#it{e^{a+bM}+c}", "exp([0]+[1]*x)+[2]"));
    //fffs.push_back(make_pair("#it{e^{a+bM}+cM^{2}}", "exp([0]+[1]*x)+[2]*x*x"));
    //fffs.push_back(make_pair("#it{e^{a+bM^{C}}}", "exp([0]+[1]*pow(x,[2]))"));

    TCanvas ** cans = new TCanvas*[bghists.size()];
    int fcount = 0;

    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        cout << "Model: " << ihbg->label << endl;
        TF1 **f = new TF1*[fffs.size()], *f1, *f2, *f3;
        TLegend *leg = new TLegend(0.55, 0.73, 0.94, 0.94);
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
        cans[fcount]->SetLogy(islog);

        TH1* h = (TH1*)ihbg->hist->Clone();
        h->SetLineColor(kBlack);
        h->SetMarkerStyle(21);
        h->SetMarkerColor(kBlack);
        h->SetFillColor(kWhite);
        h->Rebin(rebin);
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitle();
        h->GetYaxis()->SetTitle();

        leg->AddEntry(h, ihbg->label.c_str());

        for(unsigned int ifunc = 0; ifunc < fffs.size(); ifunc++)
        {
            char tmp2[32];
            sprintf(tmp2, "f%d", ifunc);
            f[ifunc] = new TF1(tmp, fffs[ifunc].second.c_str(), fitmin, 4000);
            if(f[ifunc]->GetNpar() > 2 && strstr(fffs[ifunc].first.c_str(), "#it{e^{a+bM}+c}") != NULL)
            {
                f[ifunc]->SetParLimits(2, 0.0, 1000);
                f[ifunc]->SetParameter(2, 2);
            }
            h->Fit(f[ifunc], "LMNQ", "", fitmin, 2800);
            cout << fffs[ifunc].second << " Chi^2/NDOF:" << f[ifunc]->GetChisquare() << " / " << f[ifunc]->GetNDF() << endl;
            cout << "par0: " << f[ifunc]->GetParameter(0) << "\t par1: " << f[ifunc]->GetParameter(1) << "\tpar1error: " << f[ifunc]->GetParError(1) << "\t par2: " << f[ifunc]->GetParameter(2) << endl;
            f[ifunc]->SetLineColor(colors[ifunc % NCOLORS]);
            f[ifunc]->SetLineWidth(2);
            f[ifunc]->SetMarkerStyle(0);
        }

        TF1 *fu2 = new TF1("fu2", "[0]*exp([1]*(x-[2]))", fitmin, 4000);
        fu2->SetParameters(f[0]->Eval(600)*(f[0]->GetParameter(1) + f[0]->GetParError(1)) / f[0]->GetParameter(1), f[0]->GetParameter(1) + f[0]->GetParError(1), 600);
        TF1 *fd2 = new TF1("fd2", "[0]*exp([1]*(x-[2]))", fitmin, 4000);
        fd2->SetParameters(f[0]->Eval(600)*(f[0]->GetParameter(1) - f[0]->GetParError(1)) / f[0]->GetParameter(1), f[0]->GetParameter(1) - f[0]->GetParError(1), 600);
        fu2->SetLineWidth(2);
        fu2->SetLineStyle(2);
        fu2->SetMarkerStyle(0);
        fd2->SetLineWidth(2);
        fd2->SetLineStyle(2);
        leg->AddEntry(fu2, "slope systematic");

        f1 = new TF1("f1", "exp([0]+[1]*x)", fitmin, 4000);
        f2 = new TF1("f2", "exp([0]+[1]*x)", fitmin, 4000);
        f3 = new TF1("f3", "exp([0]+[1]*x)", fitmin, 4000);
        h->Fit(f1, "LMNQ", "", fitmin, 2800);
        h->Fit(f2, "LMNQB", "", fitmin, 1400);
        h->Fit(f3, "LMNQB", "", 1400, 2800);
        double s1 = f1->GetParameter(1), s2 = f2->GetParameter(1), s3 = f3->GetParameter(1);
        //double e1 = f1->GetParError(1), e2 = f2->GetParError(1), e3 = f3->GetParError(1);
        //double diff1 = fabs(s1 - s2), diff2 = fabs(s2 - s3), diff1e = e1 * e1 + e2*e2, diff2e = e3 * e3 + e2*e2;
        //double slopeError = fabs((diff1 / diff1e + diff2 / diff2e) / (1 / diff1e + 1 / diff2e));

        TF1 *fu = new TF1("fu", "[0]*exp([1]*(x-[2]))", fitmin, 4000);
        fu->SetParameters(f1->Eval(600), max(s1, max(s2, s3)), 600);
        TF1 *fd = new TF1("fd", "[0]*exp([1]*(x-[2]))", fitmin, 4000);
        fd->SetParameters(f1->Eval(600), min(s1, min(s2, s3)), 600);

        TGraphAsymmErrors *tgshape = new TGraphAsymmErrors(), *tgslope = new TGraphAsymmErrors(), *tgtotal = new TGraphAsymmErrors();
        tgshape->SetFillColor(kGreen - 1);
        tgslope->SetFillColor(kBlue + 1);
        tgtotal->SetFillColor(kYellow);
        tgshape->SetLineColor(kGreen - 1);
        tgslope->SetLineColor(kBlue + 1);
        tgtotal->SetLineColor(kYellow);
        tgshape->SetMarkerStyle(0);
        tgslope->SetMarkerStyle(0);
        tgtotal->SetMarkerStyle(0);
        //leg->AddEntry(tgshape, "shape systematic");
        //leg->AddEntry(tgslope, "slope systematic");
        leg->AddEntry(tgtotal, "shape systematic");

        int istart = h->FindBin(600);
        double integral = 0;
        if(!ihbg->label.compare("Other")) fprintf(fbackfits, "%s,%s", "Other", "bgest");
        else if(!ihbg->label.compare("t#bar{t}")) fprintf(fbackfits, "%s,%s", "TT", "bgest");
        else if(!ihbg->label.compare("Z+Jets")) fprintf(fbackfits, "%s,%s", "ZJ", "bgest");
        else fprintf(fbackfits, "%s,%s", ihbg->label.c_str(), "bgest");
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
            if(false)//!ihbg->label.compare("Other"))
            {
                tgshape->SetPointError(i - istart, 0.0, 0.0, 0.0, 0.0);
                tgslope->SetPointError(i - istart, 0.0, 0.0, ynom - fd->Eval(x), fu->Eval(x) - ynom);
                edo = std::max(ynom - fd->Eval(x), fabs(f[0]->Eval(x) - fd2->Eval(x)));
                eup = std::max(fu->Eval(x) - ynom, fabs(fu2->Eval(x) - f[0]->Eval(x)));
                tgtotal->SetPointError(i - istart, 0.0, edo, eup);
            }
            else
            {
                tgshape->SetPointError(i - istart, 0.0, 0.0, ynom - min, max - ynom);
                tgslope->SetPointError(i - istart, 0.0, 0.0, ynom - fd->Eval(x), fu->Eval(x) - ynom);
                edo = std::max(sqrt(pow(ynom - fd->Eval(x), 2) + pow(ynom - min, 2)), fabs(f[0]->Eval(x) - fd2->Eval(x)));
                eup = std::max(sqrt(pow(fu->Eval(x) - ynom, 2) + pow(max - ynom, 2)), fabs(fu2->Eval(x) - f[0]->Eval(x)));
                tgtotal->SetPointError(i - istart, 0.0, 0.0, edo, eup);
            }
            if(fmod(h->GetBinLowEdge(i), 200) == 0.0)
            {
                if(h->GetBinLowEdge(i) > 600) fprintf(fbackfits, ",%f", integral);
                integral = 0.0;
            }
            //integral += (ynom + eup + ynom - edo) / 2;
            integral += ynom;
        }
        fprintf(fbackfits, ",%f\n", integral);

        double i1 = 0.0, i2 = 0.0, x, y;
        char histlabelforfile[128];
        if(!ihbg->label.compare("Other"))
        {
            fprintf(fbackfits, "%s,%s", "Other", "shape");
            sprintf(histlabelforfile, "%s", "Other");
        }
        else if(!ihbg->label.compare("t#bar{t}"))
        {
            fprintf(fbackfits, "%s,%s", "TT", "shape");
            sprintf(histlabelforfile, "%s", "TT");
        }
        else if(!ihbg->label.compare("Z+Jets"))
        {
            fprintf(fbackfits, "%s,%s", "ZJ", "shape");
            sprintf(histlabelforfile, "%s", "ZJ");
        }
        else
        {
            fprintf(fbackfits, "%s,%s", ihbg->label.c_str(), "shape");
            sprintf(histlabelforfile, "%s", ihbg->label.c_str());
        }
        for(int i = istart; i <= h->GetNbinsX(); i++)
        {
            if(fmod(h->GetBinLowEdge(i), 200) == 0.0)
            {
                if(h->GetBinLowEdge(i) > 600) fprintf(fbackfits, ",%f", 2 * i1 / (i1 + i2));
                i1 = 0.0;
                i2 = 0.0;
            }
            tgtotal->GetPoint(i - istart, x, y);
            i1 += y + tgtotal->GetErrorYhigh(i - istart);
            i2 += y - tgtotal->GetErrorYlow(i - istart);
        }
        fprintf(fbackfits, ",%f\n", 2 * i1 / (i1 + i2));

        TH1 *dummy = new TH1F("dummy", "dummy", 1000, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()));
        dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
        if(xmin != xmax) dummy->GetXaxis()->SetRangeUser(xmin, xmax);
        if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);
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
        //tgshape->Draw("E3 same");
        //tgslope->Draw("E3Top same");
        fu2->Draw("same");
        fd2->Draw("same");
        h->Draw("same");

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
            leg->AddEntry(f[0], "exp fit");
        }
        fixOverlay();
        leg->Draw("same");
        mark.DrawLatex(0.15, 0.95, "CMS Preliminary");
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} at 8 TeV", iLumi / 1000);
        mark.DrawLatex(0.65, 0.95, lumistamp);
        //mark.DrawLatex(0.645, 0.615, tmp4);

        if(saveplots)
        {

            char ofn[128], ofn2[128];
            sprintf(ofn, "%s_%s_%s%s.pdf", formlabel.c_str(), histlabelforfile, datahist.hist->GetName(), islog?"":"_linear");
            cans[fcount]->Print(ofn);
            sprintf(ofn2, "%s_%s_%s%s.png", formlabel.c_str(), histlabelforfile, datahist.hist->GetName(), islog?"":"_linear");
            cans[fcount]->Print(ofn2);
        }

        fcount++;
        delete [] f;
    }
    fclose(fbackfits);
}

void HnuPlots::plotQCDFit()
{
    setTDRStyle();

    TH1* hd = (TH1*)datahist.hist->Clone();
    for(std::vector<HnuPlots::HistStruct>::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {

        hd->Add(ihbg->hist, -1);
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

void HnuPlots::plotNorm(double lower, double upper, bool flip)
{
    //TH1* hdata = (TH1*)datahist.second->Clone();
    TH1* hrenorm = (TH1*)bghists[0].hist->Clone();
    TH1* hdatamod = (TH1*)datahist.hist->Clone();
    for(std::vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); )
    {
        ihbg++;
        if(ihbg != bghists.end())
        {
            TH1* hbg = (TH1*)ihbg->hist->Clone();
            hdatamod->Add(hbg, -1);
        }
    }

    double dataerror, mcerror;
    double dataint = hdatamod->IntegralAndError(hdatamod->FindBin(lower), hdatamod->FindBin(upper), dataerror);
    double mcint = hrenorm->IntegralAndError(hrenorm->FindBin(lower), hrenorm->FindBin(upper), mcerror);

    if(flip)
    {
        std::cout << "Normalization factor: " << mcint << "/" << dataint << " = " << mcint / dataint << " +/- "
                << mcint / dataint * sqrt(dataerror * dataerror / (dataint * dataint) + mcerror * mcerror / (mcint * mcint)) << std::endl;

        datahist.hist->Scale(mcint / dataint);
    }
    else
    {

        std::cout << "Normalization factor: " << dataint << "/" << mcint << " = " << dataint / mcint << " +/- "
                << dataint / mcint * sqrt(dataerror * dataerror / (dataint * dataint) + mcerror * mcerror / (mcint * mcint)) << std::endl;

        bghists[0].hist->Scale(dataint / mcint);
    }
    plot();
}

void HnuPlots::plotRatios()
{
    using namespace std;

    //gROOT->SetStyle("Plain");
    setTDRStyle();
    gStyle->SetHatchesSpacing(0.7);

    if(rebin > 1)
    {
        datahist.hist->Rebin(rebin);
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist->Rebin(rebin);
        }
    }
    else if(rebin < 0)
    {
        datahist.hist = datahist.hist->Rebin(sizeof(bins) / sizeof(double) - 1, "mWR_limitbins" , bins);
        for(int i = 1; i <= datahist.hist->GetNbinsX(); i++)
        {
            datahist.hist->SetBinContent(i, datahist.hist->GetBinContent(i) / datahist.hist->GetBinWidth(i));
            datahist.hist->SetBinError(i, datahist.hist->GetBinError(i) / datahist.hist->GetBinWidth(i));
        }
        for(vector<HnuPlots::HistStruct>::iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist = ihbg->hist->Rebin(sizeof(bins) / sizeof(double) - 1, "blarg" , bins);
            for(int i = 1; i <= ihbg->hist->GetNbinsX(); i++)
            {
                ihbg->hist->SetBinContent(i, ihbg->hist->GetBinContent(i) / ihbg->hist->GetBinWidth(i));
            }
        }
        for(vector<HnuPlots::HistStruct>::iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
        {
            ihsig->hist = ihsig->hist->Rebin(sizeof(bins) / sizeof(double) - 1, "splat", bins);
            for(int i = 1; i <= ihsig->hist->GetNbinsX(); i++)
            {
                ihsig->hist->SetBinContent(i, ihsig->hist->GetBinContent(i) / ihsig->hist->GetBinWidth(i));
            }
        }
    }

    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 8 TeV", iLumi / 1000);

    if(autosort) sort(bghists.begin(), bghists.end(), compHistInt);

    TCanvas *c1;
    double fontScale = 1.0;
    c1 = new TCanvas("c1", "c1", 800, 800);
    //c1->Divide(1, 1);
    //c1->cd(1);
    //gPad->SetPad("p1", "p1", 0, 0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.15);
    fontScale = 8.0 / 9;
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);

    //TLegend *leg = new TLegend(0.52, 0.67, 0.94, 0.91);
    TLegend *leg = new TLegend(0.18, 0.67, 0.54, 0.91);
    leg->SetFillStyle(0); //Color(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
    if(xmin != xmax) dummy->GetXaxis()->SetRangeUser(xmin, xmax);
    //dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetYaxis()->SetRangeUser(0,3);
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.05 / fontScale);
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->SetStats(0);
    dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    TLatex mark;
    mark.SetTextSize(0.04 * 7 / 6.5 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    fixOverlay();
    dummy->Draw();
    fixOverlay();

    mark.DrawLatex(0.15, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.68, 0.95, lumistamp);
    fixOverlay();

    for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
    {
        double ratio = 0, error = 0;
        TH1 *chdata = (TH1*)datahist.hist->Clone();
        for(int i = 1; i <= chdata->GetNbinsX(); i++)
        {
            if(chdata->GetBinContent(i) != 0 && ihbg->hist->GetBinContent(i) != 0)
            {
                ratio = chdata->GetBinContent(i) / ihbg->hist->GetBinContent(i);
                error = ratio * chdata->GetBinError(i) / chdata->GetBinContent(i);
                chdata->SetBinContent(i, ratio);
                chdata->SetBinError(i, error);
            }
            else
            {
                chdata->SetBinContent(i, -1);
                chdata->SetBinError(i, 0);
            }
            chdata->SetLineColor(ihbg->hist->GetLineColor());
        }
        for(int i = 1; i <= chdata->GetNbinsX(); i++) printf("%f\t", chdata->GetBinContent(i));
        printf("\n");
        chdata->Draw("same hist");
        chdata->SetMarkerStyle(0);
        leg->AddEntry(chdata, ihbg->label.c_str());
    }
    fixOverlay();
    leg->Draw("same");
    fixOverlay();

    if(false)//saveplots)
    {
        char ofn[128], ofn2[128], tmp[32];
        int cutlevel = 11111;
        char * pos = strstr(datahist.hist->GetTitle(), "cut:");
        int scans = pos?sscanf(pos, "cut:%d", &cutlevel):0;
        if(scans == 0) cutlevel = 111111;
        sprintf(tmp, "cut:%da", cutlevel);
        string tmp2 = datahist.hist->GetName();
        if(tmp2.find("_") < tmp2.size()) tmp2.erase(tmp2.find("_"), tmp2.size());
        if(strstr(datahist.hist->GetTitle(), tmp) != NULL) sprintf(ofn, "%s_cut%da_%s%s.pdf", formlabel.c_str(), cutlevel, tmp2.c_str(), islog?"":"_linear");
        else sprintf(ofn, "%s_cut%d_%s%s.pdf", formlabel.c_str(), cutlevel, tmp2.c_str(), islog?"":"_linear");
        c1->Print(ofn);
        if(strstr(datahist.hist->GetTitle(), tmp) != NULL) sprintf(ofn2, "%s_cut%da_%s%s.png", formlabel.c_str(), cutlevel, tmp2.c_str(), islog?"":"_linear");
        else sprintf(ofn2, "%s_cut%d_%s%s.png", formlabel.c_str(), cutlevel, tmp2.c_str(), islog?"":"_linear");
        c1->Print(ofn2);
    }
}

void HnuPlots::scaleByShape(double llow, double lhigh, int npar)
{
    HnuPlots::fitfunction fit;
    fit.npar = npar;
    fit.bghs = bghists;

    TF1 *ff = new TF1("fit", fit , 0, 4000, npar);

    TH1 * hd = (TH1*)datahist.hist->Clone("histtofit");
    hd->Fit(ff, "LNQM", "", llow, lhigh);

    for(int i = 0; i < fit.npar; i++)
    {
        std::cout << "Normalization factor - " << bghists[i].label << " : " << ff->GetParameter(i) << " +/-" << ff->GetParError(i) << std::endl;
        bghists[i].hist->Scale(ff->GetParameter(i));
    }

    plot();
}

void HnuPlots::cutFlow()
{
    std::vector<HnuPlots::HistStruct > hists;
    hists.push_back(datahist);
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ++ihbg)
    {
        hists.push_back(*ihbg);
    }
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ++ihbg)
    {
        hists.push_back(*ihbg);
    }

    for(std::vector<HnuPlots::HistStruct >::const_iterator i = hists.begin(); i != hists.end(); i++)
    {
        if(!i->label.compare("Data")) printf("%6s,cutlevel", "data");
        else if(!i->label.compare("t#bar{t}")) printf("%6s,cutlevel", "ttbar");
        else if(!i->label.compare("Z+Jets")) printf("%6s,cutlevel", "zjets");
        else if(!i->label.compare("Other")) printf("%6s,cutlevel", "other");
        else if(i->label.compare("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} =") > 0)
        {
            float mass;
            sscanf(i->label.c_str(), "M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = %f TeV", &mass);
            printf("signal_%i_%i,cutlevel", (int)(1000 * mass), (int)(1000 * mass) / 2);
        }
        else printf("%10s,cutlevel", i->label.c_str());
        for(int icl = 1; icl < 9; icl++)
        {
            //if(i == hists.begin()) printf("%18s,", i->hist->GetXaxis()->GetBinLabel(icl));

            double binval = i->hist->GetBinContent(icl);
            printf(",%12.3f", binval);
        }
        printf("\n");
    }
}

void HnuPlots::sigEff()
{
    //printf("%18s & ", "Mass Point");
    /*for(std::vector<HnuPlots::HistStruct >::const_iterator i = sighists.begin(); i != sighists.end(); i++)
    {
    if(i != sighists.end() - 1) printf("%10s & ", i->first.c_str());
    else printf("%10s \\\\ \\hline\n", i->first.c_str());
    }
    for(std::vector<HnuPlots::HistStruct >::const_iterator i = sighists.begin(); i != sighists.end(); i++)
    {
    if(i == sighists.begin()) printf("%18s & ", i->first.c_str());
    if(i != sighists.end() - 1) printf("%10.3f & ", (i->second->GetBinContent(6) / i->second->GetBinContent(1))/sigscale);
    else printf("%10.3f \\\\ \\hline\n", (i->second->GetBinContent(6) / i->second->GetBinContent(1))/sigscale);
    }*/
    /*for(std::vector<HnuPlots::HistStruct >::const_iterator i = sighists.begin(); i != sighists.end(); i++)
    {
    printf("%s,%s", i->first.c_str(), "sigeff");
    for(int j = 600; j < 4000; j+=200) printf(",%f", (i->second->GetBinContent(6) / i->second->GetBinContent(1))/sigscale);
    printf("\n");
    }*/
    for(std::vector<HnuPlots::HistStruct >::const_iterator i = sighists.begin(); i != sighists.end(); i++)
    {
        TH1 *h = (TH1*)i->hist->Clone();
        //double bins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 4000.0};
        TH1 *h2 = h->Rebin(sizeof(bins) / sizeof(double) - 1, "newhistforeffcalc", bins);
        // Add lost overflow into last bin
        h2->SetBinContent(h2->GetNbinsX(), h2->GetBinContent(h2->GetNbinsX()) + h2->GetBinContent(h2->GetNbinsX() + 1));
        printf("%s,%s", i->label.c_str(), "eff_2012");
        for(int j = h2->FindBin(bins[0] + 1); j <= h2->GetNbinsX(); j++)
        {
            //std::cout << "<" << h2->GetBinLowEdge(j) << ">";

            printf(",%f", h2->GetBinContent(j));
        }
        printf("\n");
    }
}

void HnuPlots::integrals(double min, double max, double* passnum, double* err)
{
    if(passnum)
    {
        TH1* hdatamod = (TH1*)datahist.hist->Clone();
        for(std::vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            TH1* hbg = (TH1*)ihbg->hist->Clone();
            hdatamod->Add(hbg, -1);
        }
        *passnum = hdatamod->IntegralAndError(hdatamod->FindBin(min), hdatamod->FindBin(max), *err);
    }
    else
    {
        printf("Integrals\n%s: %f\n", datahist.label.c_str(), datahist.hist->Integral(datahist.hist->FindBin(min), datahist.hist->FindBin(max)));
        for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
        {
            printf("%s: %f\n", ihbg->label.c_str(), ihbg->hist->Integral(ihbg->hist->FindBin(min), ihbg->hist->FindBin(max)));
        }
        for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {

            printf("%s: %f\n", ihbg->label.c_str(), ihbg->hist->Integral(ihbg->hist->FindBin(min), ihbg->hist->FindBin(max)));
        }
    }
}

void HnuPlots::sigRMS()
{
    printf("Signal RMS:\n");
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {

        printf("%s: %f\n", ihbg->label.c_str(), ihbg->hist->GetRMS());
    }
}

void HnuPlots::sigStatErr()
{
    printf("Signal Stat Err:\n");
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ihbg++)
    {

        double integral, err;
        integral = ihbg->hist->IntegralAndError(0, ihbg->hist->GetNbinsX() + 1, err);
        printf("%s: %f\n", ihbg->label.c_str(), 1 + err / integral);
    }
}

void HnuPlots::sigMatch()
{
    printf("Signal Matching ratios:\n");
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = sighists.begin(); ihbg != sighists.end(); ++ihbg)
    {

        printf("%s: %f\n", ihbg->label.c_str(), ihbg->hist->GetBinContent(3) / ihbg->hist->Integral(0, ihbg->hist->GetNbinsX() + 1));
    }
}

void HnuPlots::mcBgShape(int cutlevel, std::string sample)
{
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ++ihbg)
    {
        // bins contains the chosen binning for the limit setting code
        //double bins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 4000.0};
        TH1 *hshape = new TH1F("hshape", "hshape", sizeof(bins) / sizeof(double) - 1, bins);
        TH1 *hshapen = new TH1F("hshapen", "hshapen", sizeof(bins) / sizeof(double) - 1, bins);
        TH1 *hshaped = new TH1F("hshaped", "hshaped", sizeof(bins) / sizeof(double) - 1, bins);
        TH1 *herr = new TH1F("herr", "herr", sizeof(bins) / sizeof(double) - 1, bins);
        std::map<double, std::vector<HeavyNuTree::HNuSlopeFitInfo> > avgWgtMap;
        // Extract interesting data from tuple, second element of pair is total sub-sample k-factor
        for(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator it = ihbg->fittree.begin(); it != ihbg->fittree.end(); ++it)
        {
            if(it->first.cutlevel >= cutlevel && it->first.cutlevel < 10 && it->first.weight < 1e3 && it->first.weight > 1e-3)
            {
                hshape->Fill(it->first.mlljj, it->first.weight * it->second);
                herr->Fill(it->first.mlljj);
                avgWgtMap[it->second].push_back(it->first);
            }
        }
        // Use reorganized tuple info to caluclate average weight
        for(std::map<double, std::vector<HeavyNuTree::HNuSlopeFitInfo> >::const_iterator mapit = avgWgtMap.begin(); mapit != avgWgtMap.end(); ++mapit)
        {
            TH1 *hN = new TH1F("hN", "hN", sizeof(bins) / sizeof(double) - 1, bins);
            TH1 *hW = new TH1F("hW", "hW", sizeof(bins) / sizeof(double) - 1, bins);
            for(std::vector<HeavyNuTree::HNuSlopeFitInfo>::const_iterator sit = mapit->second.begin(); sit != mapit->second.end(); ++sit)
            {
                //hN->Fill(sit->mlljj);
                //hW->Fill(sit->mlljj, sit->weight);
                hW->Fill(sit->mlljj, sit->weight * mapit->first);
                hN->Fill(sit->mlljj, sit->weight * sit->weight * mapit->first * mapit->first);
            }
            // Here we calculate the poisson error and scale by the average weight to get the gamma error, this is then scaled by the final event expectation
            for(int i = 0; i <= hN->GetNbinsX(); i++)
            {
                if(hN->GetBinContent(i) > 0)// && mapit->first > 0)
                {
                    //hshapen->SetBinContent(i, hshapen->GetBinContent(i) + hN->GetBinContent(i) * fabs(mapit->first) * sqrt(hN->GetBinContent(i)) * fabs(mapit->first) * (hW->GetBinContent(i) / hN->GetBinContent(i)));
                    //hshaped->SetBinContent(i, hshaped->GetBinContent(i) + hN->GetBinContent(i) * fabs(mapit->first));
                    hshapen->SetBinContent(i, hshapen->GetBinContent(i) + hN->GetBinContent(i));
                    hshaped->SetBinContent(i, hshaped->GetBinContent(i) + hW->GetBinContent(i));
                }
            }
            hN->~TH1();
            hW->~TH1();
        }
        std::string label;
        if(!ihbg->label.compare("Other"))
        {
            label = "other" + sample;
        }
        else if(!ihbg->label.compare("t#bar{t}"))
        {
            label = "ttjets" + sample;
        }
        else if(!ihbg->label.compare("Z+Jets"))
        {
            label = "zjets" + sample;
        }
        else
        {
            label = ihbg->label + sample;
        }
        printf("%s,2012", label.c_str());
        for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", hshape->GetBinContent(i));
        printf("\n%s,GAMMASTATS", label.c_str());
        //for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", pow(hshapen->GetBinContent(i) / hshaped->GetBinContent(i), 2) / hshape->GetBinContent(i));
        for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", hshapen->GetBinContent(i) / hshaped->GetBinContent(i));
        printf("\n%s,n_eff", label.c_str());
        for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", hshape->GetBinContent(i) / (hshapen->GetBinContent(i) / hshaped->GetBinContent(i)));
        printf("\n%s,n_evt", label.c_str());
        for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", herr->GetBinContent(i));
        printf("\n%s,gamerr", label.c_str());
        for(int i = 1; i <= hshape->GetNbinsX(); i++) printf(",%f", 1 + sqrt(hshape->GetBinContent(i) / (hshapen->GetBinContent(i) / hshaped->GetBinContent(i))) * (hshapen->GetBinContent(i) / hshaped->GetBinContent(i)) / hshape->GetBinContent(i));
        printf("\n");

        hshape->~TH1();
        hshapen->~TH1();
        hshaped->~TH1();
        herr->~TH1();
    }

    //double bins2[] = {0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 4000.0};
    //TH1 *hshape2 = new TH1F("hshape", "hshape", sizeof(bins2) / sizeof(double) - 1, bins2);
    //TH1 *hshapen2 = new TH1F("hshapen", "hshapen", sizeof(bins2) / sizeof(double) - 1, bins2);
    //TH1 *hshaped2 = new TH1F("hshaped", "hshaped", sizeof(bins2) / sizeof(double) - 1, bins2);
    //TH1 *herr2 = new TH1F("herr", "herr", sizeof(bins2) / sizeof(double) - 1, bins2);
    //
    //for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ++ihbg)
    //{
    //    //std::cout << ihbg->label << "\t" << ihbg->fittree.size() << std::endl;
    //    for(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator it = ihbg->fittree.begin(); it != ihbg->fittree.end(); ++it)
    //    {
    //        if(it->first.cutlevel >= cutlevel && it->first.cutlevel < 10 && it->first.weight < 1000 && it->first.weight > 1e-10)
    //        {
    //            //if(it->first.mlljj < 200) std::cout << it->first.mlljj << "\t" << it->first.cutlevel << "\t" << it->first.weight << std::endl;
    //            hshape2->Fill(it->first.mlljj, it->first.weight * it->second);
    //            hshapen2->Fill(it->first.mlljj, it->first.weight * it->second * it->second);
    //            hshaped2->Fill(it->first.mlljj, it->second);
    //            herr2->Fill(it->first.mlljj, (it->first.weight * it->second));
    //        }
    //    }
    //}
    //
    ////printf("\n%s,realgamerr", "");
    ////for(int i = 1; i <= hshape2->GetNbinsX(); i++) printf(",%f", 1 + (sqrt(hshape2->GetBinContent(i) * hshapen2->GetBinContent(i) / hshaped2->GetBinContent(i)) / hshape2->GetBinContent(i)));
    ////printf("\n%s,gamerr", "");
    //if(shapeerr.size() > 0) shapeerr.erase(shapeerr.begin(), shapeerr.end());
    //for(int i = 1; i <= hshape2->GetNbinsX(); i++)
    //{
    //
    //    shapeerr.push_back(((hshape2->GetBinContent(i) > 1e-3)?(1 + (sqrt(herr2->GetBinContent(i))) / hshape2->GetBinContent(i)):1));
    //    std::cout << sqrt(herr2->GetBinContent(i)) << "\t" << hshape2->GetBinContent(i) << "\t" << ((hshape2->GetBinContent(i) > 1e-3 && hshape2->GetBinContent(i) < 1e10)?(1 + (sqrt(herr2->GetBinContent(i))) / hshape2->GetBinContent(i)):1) << std::endl;
    //}
    ////printf("\n");
    //
    //hshape2->~TH1();
    //hshapen2->~TH1();
    //hshaped2->~TH1();
    //herr2->~TH1();
}

void HnuPlots::loadSystFile(std::string systfile, std::string ratefile, bool includeBG)
{
    using namespace std;
    map<pair<string, string>, vector<float> > uncerts;

    FILE * inFile = fopen(systfile.c_str(), "r");
    char buff[4096];
    char *c;
    while(!feof(inFile) && (c = fgets(buff, 4096, inFile)) != NULL)
    {
        char sample[128], syst[128];
        float nums[20];
        int numRead;
        for(char* k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
        if(buff[0] != '#' && (numRead = sscanf(buff, "%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", sample, syst, nums, nums + 1, nums + 2, nums + 3, nums + 4, nums + 5, nums + 6, nums + 7, nums + 8, nums + 9, nums + 10, nums + 11, nums + 12, nums + 13, nums + 14, nums + 15, nums + 16, nums + 17, nums + 18, nums + 19)) >= 3)
        {
            vector<float> tmp;
            //cout << buff << sample << "\t" << syst << "\t";
            for(int i = 0; i < numRead - 2; i++)
            {
                tmp.push_back(nums[i] - 1);
                //std::cout << nums[i] << "\t";
            }
            uncerts[make_pair(string(sample), string(syst))] = tmp;
            //cout << std::endl;
        }
    }
    fclose(inFile);

    inFile = fopen(ratefile.c_str(), "r");
    while(!feof(inFile) && (c = fgets(buff, 4096, inFile)) != NULL)
    {
        char sample[128], syst[128];
        float nums[20];
        int numRead;
        for(char* k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
        if(buff[0] != '#' && (numRead = sscanf(buff, "%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", sample, syst, nums, nums + 1, nums + 2, nums + 3, nums + 4, nums + 5, nums + 6, nums + 7, nums + 8, nums + 9, nums + 10, nums + 11, nums + 12, nums + 13, nums + 14, nums + 15, nums + 16, nums + 17, nums + 18, nums + 19)) >= 3)
        {
            vector<float> tmp;
            //cout << buff << sample << "\t" << syst << "\t";
            for(int i = 0; i < numRead - 2; i++)
            {
                tmp.push_back(nums[i]);
                //std::cout << nums[i] << "\t";
            }
            uncerts[make_pair(string(sample), string(syst))] = tmp;
            //cout << std::endl;
        }
    }
    fclose(inFile);

    string samples[] = {"ttjets", "zjets", "other"};
    int nSamples = sizeof(samples) / sizeof(string);
    string systsd[] = {"gamerr", "shape", "norm", "pdf", "fact", "ren", "lumi", "jes", "jer", "mer", "muonid", "trig", "pu", "id", "escale"};
    int nSystsd = sizeof(systsd) / sizeof(string);
    string systss[] =           {"norm", "pdf", "fact", "ren", "lumi", "jes", "jer", "mer", "muonid", "trig", "pu", "id", "escale"};
    int nSystss = sizeof(systss) / sizeof(string);


    vector<float> domsysts, sdomsysts;
    vector<pair<double, double> > domSystTmp(uncerts[make_pair(samples[0], "2012")].size(), make_pair(0.0, 0.0));
    vector<pair<double, double> > sdomSystTmp(uncerts[make_pair(samples[0], "2012")].size(), make_pair(0.0, 0.0));
    for(int iSample = 0; iSample < nSamples; iSample++)
    {
        pair<string, string> wtag(samples[iSample], "2012");
        for(unsigned int ibin = 0; ibin < uncerts[wtag].size(); ibin++)
        {
            double weight = uncerts[wtag][ibin];
            domSystTmp[ibin].second += weight;
            for(int isystd = 0; isystd < nSystsd; isystd++)
            {
                pair<string, string> id(samples[iSample], systsd[isystd]);
                if(uncerts.find(id) == uncerts.end()) continue;
                double uncert = uncerts[id][(ibin < uncerts[id].size())?ibin:(uncerts[id].size() - 1)];
                domSystTmp[ibin].first += pow(weight * uncert, 2);
            }

            sdomSystTmp[ibin].second += weight;
            for(int isysts = 0; isysts < nSystss; isysts++)
            {
                pair<string, string> id(samples[iSample], systss[isysts]);
                if(uncerts.find(id) == uncerts.end()) continue;
                double uncert = uncerts[id][(ibin < uncerts[id].size())?ibin:(uncerts[id].size() - 1)];
                sdomSystTmp[ibin].first += pow(weight * uncert, 2);
            }
        }
    }
    for(vector<pair<double, double> >::const_iterator ibin = domSystTmp.begin(); ibin != domSystTmp.end(); ++ibin)
    {
        domsysts.push_back(sqrt(ibin->first) / ibin->second);
    }
    systematics.push_back(domsysts);

    for(vector<pair<double, double> >::const_iterator ibin = sdomSystTmp.begin(); ibin != sdomSystTmp.end(); ++ibin)
    {

        sdomsysts.push_back(sqrt(ibin->first) / ibin->second);
    }
    systematics.push_back(sdomsysts);

    if(includeBG) histFromDataCard(uncerts);
}

void HnuPlots::histFromDataCard(std::map<std::pair<std::string, std::string>, std::vector<float> >& uncerts)
{
    std::string samples[] = {"ttjets", "zjets", "other"};
    int nSamples = 3;

    for(int i = 0; i < nSamples; i++)
    {
        std::pair<std::string, std::string> wtag(samples[i], "2012");
        std::string samplename("samplename");
        TH1 *hist = new TH1D(samples[i].c_str(), "mWR", sizeof(bins)/sizeof(double) - 1, bins);
        for(int j = 1; j <= hist->GetNbinsX(); j++)
        {
            hist->SetBinContent(j, uncerts[wtag][j - 1]);
            std::cout << j << "\t" << uncerts[wtag][j - 1] << std::endl;
        }
        hist->SetFillColor(colors[i % NCOLORS]);
        hist->SetFillStyle(hatchs[i % NHATCHS]);
        hist->SetLineColor(colors[i % NCOLORS]);
        hist->SetMarkerColor(colors[i % NCOLORS]);
        hist->SetMarkerStyle(0);
        hist->SetLineWidth(1);

        if    (samples[i].find("ttjets") < samples[i].size()) samplename = "t#bar{t}";
        else if(samples[i].find("zjets") < samples[i].size()) samplename = "Z+Jets";
        else if(samples[i].find("other") < samples[i].size()) samplename = "Other";
        bghists.push_back(HnuPlots::HistStruct(samplename, hist));
    }
}

void HnuPlots::mcSystCalc(int forceBins, std::map<std::string, std::vector<double> >* systMap)
{
    using namespace std;
    for(std::vector<HnuPlots::HistStruct >::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ++ihbg)
    {
        TH1 *hmod = NULL, *hnom = NULL, *hmodnorm = NULL, *hnomnorm = NULL;
        hnom = (TH1*)ihbg->hist->Clone();
        if(ihbg->normhist) hnomnorm = (TH1*)ihbg->normhist->Clone();
        ihbg++;
        if(ihbg == bghists.end()) break;

        hmod = (TH1*)ihbg->hist->Clone();
        if(ihbg->normhist) hmodnorm = (TH1*)ihbg->normhist->Clone();

        //renorm the mod hist if needed
        if(hmodnorm && hnomnorm)
        {
            //std::cout << "In range: " << ihbg->normll << " - " << ihbg->normul << "\nScaling by: " << hnomnorm->Integral(hnomnorm->FindBin(ihbg->normll), hnomnorm->FindBin(ihbg->normul)) << " / " << hmodnorm->Integral(hmodnorm->FindBin(ihbg->normll), hmodnorm->FindBin(ihbg->normul)) << " = " << hnomnorm->Integral(hnomnorm->FindBin(ihbg->normll), hnomnorm->FindBin(ihbg->normul)) / hmodnorm->Integral(hmodnorm->FindBin(ihbg->normll), hmodnorm->FindBin(ihbg->normul)) << std::endl;
            hmod->Scale(hnomnorm->Integral(hnomnorm->FindBin(ihbg->normll), hnomnorm->FindBin(ihbg->normul)) / hmodnorm->Integral(hmodnorm->FindBin(ihbg->normll), hmodnorm->FindBin(ihbg->normul)));
        }

        //double bins[] = {0.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 4000.0};
        hmod = hmod->Rebin(sizeof(bins) / sizeof(double) - 1, "hmrb", bins);
        hnom = hnom->Rebin(sizeof(bins) / sizeof(double) - 1, "hnrb", bins);

        std::vector<double> hbins;
        double error = fabs(hmod->Integral() / hnom->Integral() - 1), a = 0.0, e = 0.0, p = 0.0;
        int ilast = hnom->GetNbinsX(), ifirst = hnom->FindBin(bins[0] + 1);
        switch(forceBins)
        {
            case 0:
                hbins.push_back(bins[sizeof(bins) / sizeof(double) - 1]);
                for(int i = ilast; i >= ifirst; i--)
                {
                    a = hnom->IntegralAndError(i, ilast, e);
                    p = hmod->Integral(i, ilast);
                    error = fabs(a / p - 1);
                    //kludge to make 1000 always a systematic boundary
                    if((hnom->GetBinLowEdge(i) > 1000.0 - 0.1) && (hnom->GetBinLowEdge(i) < 1000.0 + 0.1)) hbins.push_back(hnom->GetBinLowEdge(i));
                    else if(e / a < error)
                    {
                        hbins.push_back(hnom->GetBinLowEdge(i));
                        ilast = i - 1;
                    }
                    else if(i == ifirst)
                    {
                        if(e / a < error) hbins.push_back(hnom->GetBinLowEdge(i));
                        else
                        {
                            hbins.push_back(hnom->GetBinLowEdge(i));
                        }
                    }
                    //std::cout << e / a << "\t" << error << std::endl;
                }
                break;
            case 1:
                //for(int i = 1; i <= hnom->GetNbinsX() + 1; i++) hbins.push_back(hnom->GetBinLowEdge(i));
                hbins.push_back(bins[0]);
                hbins.push_back(1000.0);
                hbins.push_back(bins[sizeof(bins) / sizeof(double) - 1]);
                break;
            case 2:
                for(int i = ifirst; i <= ilast + 1; i++)
                {
                    hbins.push_back(hnom->GetBinLowEdge(i));
                }
        }
        std::sort(hbins.begin(), hbins.end());

        TH1* hmp = hmod->Rebin(Int_t(hbins.size() - 1), "newhist", &hbins.front());
        TH1* hnp = hnom->Rebin(Int_t(hbins.size() - 1), "newhist2", &hbins.front());
        //cout << *isys << " uncertainty: " << (hmod->Integral() / hnom->Integral() - 1)*100 << endl;
        hmp->Divide(hnp);
        TH1* hist = (TH1*)hmp->Clone("thehisto");
        hist->Fit("pol0", "QN");

        if(systMap == NULL)
        {
            printf("%s", ihbg->label.c_str());
            for(int i = hnom->FindBin(bins[0] + 1); i < hnom->FindBin(bins[sizeof(bins) / sizeof(double) - 1] + 1); i++)
            {
                printf(",%f", hist->GetBinContent(hist->FindBin(hnom->GetBinLowEdge(i) + 1)));
            }
            printf("\n");
        }
        else
        {
            for(int i = hnom->FindBin(bins[0] + 1); i < hnom->FindBin(bins[sizeof(bins) / sizeof(double) - 1] + 1); i++)
            {
                (*systMap)[ihbg->label].push_back(hist->GetBinContent(hist->FindBin(hnom->GetBinLowEdge(i) + 1)));
            }
        }
    }
}

void HnuPlots::autoSetHistogramAxisTitle(int mode)
{
    std::string histValues(datahist.hist->GetName());
    std::string histValName = histValues.substr(0, histValues.find(";"));
    std::vector<std::string> histQs;
    // read variable names to plot
    for(size_t pos = 0, npos = 0; npos != size_t(-1);pos = npos + 1)
    {
        npos = histValName.find(':', pos + 1);
        histQs.push_back(histValName.substr(pos, npos - pos));
    }

    for(unsigned int i = 0; i < histQs.size(); i++)
    {
        std::string name(histQs[i]);
        std::string* axislabel;
        if(i == 0)      axislabel = &xaxislabel;
        else if(i == 1) axislabel = &yaxislabel;
        else break;
        switch(mode)
        {
            case 0:
            case 4:
            case 6:
            case 7:
                if(name.find("mWR") < name.size()) *axislabel = "M_{#mu#mujj} [TeV]";
                else if(name.find("mWR_1b") < name.size()) *axislabel = "M_{#mu#mubj} [TeV]";
                else if(name.find("mWR_2b") < name.size()) *axislabel = "M_{#mu#mubb} [TeV]";
                else if(name.find("mLL") < name.size()) *axislabel = "M_{#mu#mu} [TeV]";
                else if(name.find("mLL_1b") < name.size()) *axislabel = "M_{#mu#mu} (1 b-tag) [GeV]";
                else if(name.find("mLL_2b") < name.size()) *axislabel = "M_{#mu#mu} (2 b-tag) [GeV]";
                else if(name.find("mLLZoom") < name.size()) *axislabel = "M_{#mu#mu} [GeV]";
                else if(name.find("mLLNorm") < name.size()) *axislabel = "M_{#mu#mu} [TeV]";
                //else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{#mu_{#lower[-0.2]{1}}}} [GeV]";
                //else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{#mu_{#lower[-0.2]{2}}}} [GeV]";
                else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{1}} [GeV]";
                else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{2}} [GeV]";
                else if(name.find("ptL1") < name.size()) *axislabel = "p_{T}(#mu_{1}) [GeV]";
                else if(name.find("ptL2") < name.size()) *axislabel = "p_{T}(#mu_{2}) [GeV]";
                else if(name.find("ptJ1") < name.size()) *axislabel = "p_{T}(j_{1}) [GeV]";
                else if(name.find("ptJ2") < name.size()) *axislabel = "p_{T}(j_{2}) [GeV]";
                else if(name.find("mJJ") < name.size()) *axislabel = "M_{jj} [GeV]";
                else if(name.find("pL1") < name.size()) *axislabel = "p(#mu_{1}) [GeV]";
                else if(name.find("pL2") < name.size()) *axislabel = "p(#mu_{2}) [GeV]";
                else if(name.find("dPhiL") < name.size()) *axislabel = "#Delta#phi(#mu_{1}, #mu_{2}) [GeV]";
                else if(name.find("dEtaL") < name.size()) *axislabel = "#Delta#eta(#mu_{1}, #mu_{2}) [GeV]";
                break;
            case 1:
            case 5:
            case 8:
            case 9:
                if(name.find("mWR_1b") < name.size()) *axislabel = "M_{eebj} [TeV]";
                else if(name.find("mWR_2b") < name.size()) *axislabel = "M_{eebb} [TeV]";
                else if(name.find("mWR") < name.size()) *axislabel = "M_{eejj} [TeV]";
                else if(name.find("mLL") < name.size()) *axislabel = "M_{ee} [TeV]";
                else if(name.find("mLL_1b") < name.size()) *axislabel = "M_{ee} (1 b-tag) [GeV]";
                else if(name.find("mLL_2b") < name.size()) *axislabel = "M_{ee} (2 b-tag) [GeV]";
                else if(name.find("mLLZoom") < name.size()) *axislabel = "M_{ee} [GeV]";
                else if(name.find("mLLNorm") < name.size()) *axislabel = "M_{ee} [TeV]";
                //else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{e_{#lower[-0.2]{1}}}} [GeV]";
                //else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{e_{#lower[-0.2]{2}}}} [GeV]";
                else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{1}} [GeV]";
                else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{2}} [GeV]";
                else if(name.find("mOuR1") < name.size()) *axislabel = "M_{eej1} [GeV]";
                else if(name.find("mOuR2") < name.size()) *axislabel = "M_{eej2} [GeV]";
                else if(name.find("ptL1") < name.size()) *axislabel = "p_{T}(e_{1}) [GeV]";
                else if(name.find("ptL2") < name.size()) *axislabel = "p_{T}(e_{2}) [GeV]";
                else if(name.find("etaL1") < name.size()) *axislabel = "#eta(e_{1})";
                else if(name.find("etaL2") < name.size()) *axislabel = "#eta(e_{2})";
                else if(name.find("phiL1") < name.size()) *axislabel = "#phi(e_{1})";
                else if(name.find("phiL2") < name.size()) *axislabel = "#phi(e_{2})";
                else if(name.find("mJJ") < name.size()) *axislabel = "M_{jj} [GeV]";
                //else if(name.find("mLQmin") < name.size()) *axislabel = "min M_{LQ} [GeV]";
                else if(name.find("st") < name.size()) *axislabel = "S_{T} [GeV]";
                else if(name.find("pL1") < name.size()) *axislabel = "E(e_{1}) [GeV]";
                else if(name.find("pL2") < name.size()) *axislabel = "E(e_{2}) [GeV]";
                else if(name.find("rhL1") < name.size()) *axislabel = "seed crystal E(e_{1}) [GeV]";
                else if(name.find("rhL2") < name.size()) *axislabel = "seed crystal E(e_{2}) [GeV]";
                else if(name.find("seL1") < name.size()) *axislabel = "seed cluster E(e_{1}) [GeV]";
                else if(name.find("seL2") < name.size()) *axislabel = "seed cluster E(e_{2}) [GeV]";
                else if(name.find("rhoScL1") < name.size()) *axislabel = "seed crystal E(e_{1}) / seed cluster E(e_{1})";
                else if(name.find("rhoScL2") < name.size()) *axislabel = "seed crystal E(e_{2}) / seed cluster E(e_{2})";
                else if(name.find("dPhiL") < name.size()) *axislabel = "#Delta#phi(e_{1}, e_{2}) [GeV]";
                else if(name.find("dEtaL") < name.size()) *axislabel = "#Delta#eta(e_{1}, e_{2}) [GeV]";
                else if(name.find("jmult") < name.size()) *axislabel = "jet multiplicity";
                else if(name.find("bmult") < name.size()) *axislabel = "b jet multiplicity";
                else if(name.find("mLQmin") < name.size()) *axislabel = "M_{ej}^{min} [GeV]";
                else if(name.find("mLQmax") < name.size()) *axislabel = "M_{ej}^{max} [GeV]";
                else if(name.find("mLQavg") < name.size()) *axislabel = "M_{ej}^{avg} [GeV]";
                break;
            case 2:
            case 3:
                if(name.find("mWR") < name.size()) *axislabel = "M_{e#mujj} [TeV]";
                else if(name.find("mLL") < name.size()) *axislabel = "M_{e#mu} [GeV]";
                else if(name.find("mLLZoom") < name.size()) *axislabel = "M_{e#mu} [GeV]";
                else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{#mu}} [GeV]";
                else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{#lower[-0.2]{2}}} [GeV]";
                break;
        }
        if(name.find("mJJ") < name.size()) *axislabel = "M_{jj} [GeV]";
        else if(name.find("n_vertex") < name.size()) *axislabel = "N primary vertex";
        else if(name.find("ptJ1") < name.size()) *axislabel = "p_{T}(j_{1}) [GeV]";
        else if(name.find("ptJ2") < name.size()) *axislabel = "p_{T}(j_{2}) [GeV]";
        else if(name.find("etaJ1") < name.size()) *axislabel = "#eta(j_{1})";
        else if(name.find("etaJ2") < name.size()) *axislabel = "#eta(j_{2})";
        else if(name.find("phiJ1") < name.size()) *axislabel = "#phi(j_{1})";
        else if(name.find("phiJ2") < name.size()) *axislabel = "#phi(j_{2})";
        else if(name.find("met") < name.size()) *axislabel = "Missing E_{t} [GeV]";
        else if(name.find("dEtaJ") < name.size()) *axislabel = "#Delta#eta(j_{1}, j_{2}) [GeV]";
        else if(name.find("dPhiJ") < name.size()) *axislabel = "#Delta#phi(j_{1}, j_{2}) [GeV]";
        if(!axislabel->size()) *axislabel = "Sorry, no approperiate label found";
    }
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

void HnuPlots::setCompPlot(bool cp)
{
    plotSMoData = cp;
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

void HnuPlots::setYRange(double min, double max)
{
    ymin = min;
    ymax = max;
}

void HnuPlots::setSavePlots(bool sp)
{
    saveplots = sp;
}

//cutlevel strings
const std::string cutlevels[] = {
    "cut0_none",
    "cut1_LLJJpt",
    "cut2_TrigMatches",
    "cut3_Vertex",
    "cut4_L1HighPt",
    "cut5_diLmass",
    "cut6_mWRmass",
    "cutX_LLpt",
    "cut5a_loDiLmass",
    "cut4a_L1HighPt_1b",
    "cut4b_L1HighPt_2b",
    "cut4c_ZPeak",
    "LQ1Cuts",
    "LQ1Cuts2",
    "LQ1Cuts3",
    "LQ1Cuts4",
    "LQ1Cuts5",
    ""
};

const std::string cutlevelsTop[] = {
    "cut0_none",
    "cut1_LLJJpt",
    "cut2_TrigMatches",
    "cut3_Vertex",
    "cut4_Mu1HighPt",
    "cut5_diLmass",
    "cut6_mWRmass",
    "cutX_LLpt",
    "cut5a_loDiLmass",
    "cut4a_Mu1HighPt_1b",
    "cut4b_Mu1HighPt_2b",
    ""
};

//Integrated lumi
const static double lumi2012mm = /*5295 + 7002 + 7264*/19700, lumi2012ee = /*5280 + 7031 + 7274*/19700;
//MC xsecs - https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
const static double xsecttbar = 23.64, xsecZJ = 3503.71, xsecZZ = 17.721, xsecWZ = 33.59, xsecWW = 55.475, xsectW = 11.1, xsectbarW = 11.1, xsecWJ = 36257.2, xsecZ0J = 2950.0, xsecZ1J = 561.0, xsecZ2J = 181.0, xsecZ3J = 51.1, xsecZ4J = 23.04;
//MC total events
const static double Nttbar = 4246444, NZJ = 28807863, NZZ = 9739908, NWZ = 10000283, NWW = 10000431, NtW = 497658, NtbarW = 493460, NZ0J = 28807863, NZ1J = 23745248, NZ2J = 21521261, NZ3J = 10602630, NZ4J = 5499858;
//muon k factors
const static double k_mm_ddtop = /*0.620166*/0.631595,                           k_mm_Zscale = /*1.00858*//*1.02005*/1.02701*0.999979, k_mm_NNLOZ = 1.2036, k_top = 1.13159;
//electron k factors
const static double k_ee_ddtop = /*0.518549*/0.524452 * lumi2012ee / lumi2012mm, k_ee_Zscale = /*0.963943*//*0.939217 0.973471 1.05259*/ 1.00043*0.99879, k_ee_NNLOZ = 1.1893;

//data files
const std::string data_ee("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_2/Elec-Run2012ABCD-22Jan2013-v1.root");
const std::string data_mm("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_2/Mu-Run2012ABCD-22Jan2013-v1.root");
const std::string data_em("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_2/Mu-Run2012ABCD-22Jan2013-v1.root");

//const std::string data_ee("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_10/run2012ABCD-electron.root");
//const std::string data_mm("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_10/run2012ABCD-muon.root");
//const std::string data_em("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_10/run2012ABCD-muon.root");


//mc file
const std::string mc_tt(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_TTJets_FullLeptMGDecays_8TeV-madgraph.root");
const std::string mc_ZJ(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_START53_V7A_skim.root");
const std::string mc_ZZ(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_ZZ_TuneZ2star_8TeV_pythia6_tauola_START53_V7A-v1.root");
const std::string mc_WZ(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WZ_TuneZ2star_8TeV_pythia6_tauola_START53_V7A-v1.root");
const std::string mc_WW(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WW_TuneZ2star_8TeV_pythia6_tauola_START53_V7A-v1.root");
const std::string mc_tW(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_START53_V7A-v1.root");
const std::string mc_tbarW("/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_START53_V7A-v1.root");
const std::string mc_Z0J(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DY0JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_START53_V7A.root");
const std::string mc_Z1J(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_START53_V7A.root");
const std::string mc_Z2J(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root");
const std::string mc_Z3J(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root");
const std::string mc_Z4J(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_START53_V7A-v1_partial.root");
const std::string mc_SZJ(  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/sherpaDY_full.root");
const std::string mc_WJ(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/prelimWJets.root");

//double lumi2011AMu24 = 216, lumi2011AMu40 = 2056.0, lumi2011B = 2719.0;  //pixel only lumi
//double lumi2011AMu24 = 216.2, lumi2011AMu40 = 1956.7, lumi2011B = 2510.5;  //HF lumi

HeavyNuTree::HNuSlopeFitInfo ll, ul;

void plotMCVar(int cutlevel, std::string plot, int rebin = 5, std::string xaxis = "M_{W_{R}} [GeV]", bool rescale = false, bool log = true)
{

    using namespace std;

    //background legend label, TFile
    std::vector<std::vector<HnuPlots::FileStruct> > bg, sig;

    ll.cutlevel = cutlevel;
    ul.cutlevel = 1000;


    vector<HnuPlots::FileStruct> bgm1, bgm2, bgm3;

    //bgm1.push_back(HnuPlots::FileStruct("Gen",    "/local/cms/user/pastika/heavyNuShape/HeavyNu_accept_3000_187.root",    "hNuGen2012/cut5_diLmass/m4obj", 1.0, 1.0, 1.0,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    //bgm2.push_back(HnuPlots::FileStruct("Reco",    "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-187_TuneZ2star_8TeV-pythia6-tauola.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));

    bgm1.push_back(HnuPlots::FileStruct("P1", data_em,    "hNuEMu/"    + cutlevels[cutlevel] + "/" + "ptL1;", 1.0, 1.0, 1.0, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    bgm2.push_back(HnuPlots::FileStruct("P2", data_em,    "hNuEMu/"    + cutlevels[cutlevel] + "/" + "ptL1;SS>0.5", 1.0, 1.0, 1.0, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    
    bg.push_back(bgm1);
    bg.push_back(bgm2);
    //bg.push_back(bgm3);

    //data
    HnuPlots::FileStruct data("Data", "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-187_TuneZ2star_8TeV-pythia6-tauola.root", "hNuMu40/" + cutlevels[cutlevel] + "/" + plot);

    HnuPlots hps(data, bg, sig, 0.0);
    hps.setXAxisTitle(xaxis.c_str());
    hps.setYAxisTitle("Events");
    hps.setLog(log);
    hps.setRebin(rebin);
    hps.plotMCComp(rescale);
    hps.integrals(600,4001);
}

void setBgandData(int mode, HnuPlots::FileStruct& data, std::vector<std::vector<HnuPlots::FileStruct> >& bg, double& lumi, int cutlevel = 5, std::string plot = "mWR", bool lt = false, bool hft = false, bool isZJ = false)
{
    char fdata[256];

    std::vector<HnuPlots::FileStruct> bgTT, bgZJ, bgZ1J, bgZ2J, bgZ3J, bgZ4J, bgOther;

    ll.cutlevel = cutlevel;
    ll.mlljj    = 0.0;
    ll.mll      = 0.0;
    ll.l1pt     = 0.0;
    ll.l1eta    = -10.0;
    ll.l1phi    = -10.0;
    ll.l2pt     = 0.0;
    ll.l2eta    = -10.0;
    ll.l2phi    = -10.0;
    ll.j1pt     = 0.0;
    ll.j1eta    = -10.0;
    ll.j1phi    = -10.0;
    ll.j2pt     = 0.0;
    ll.j2eta    = -10.0;
    ll.j2phi    = -10.0;
    ul.cutlevel = 100;
    ul.mlljj    = 8000.0;
    ul.mll      = 8000.0;
    ul.l1pt     = 40000.0;
    ul.l1eta    = 10.0;
    ul.l1phi    = 10.0;
    ul.l2pt     = 40000.0;
    ul.l2eta    = 10.0;
    ul.l2phi    = 10.0;
    ul.j1pt     = 40000.0;
    ul.j1eta    = 10.0;
    ul.j1phi    = 10.0;
    ul.j2pt     = 40000.0;
    ul.j2eta    = 10.0;
    ul.j2phi    = 10.0;

    double kf = 1.0;

    //background
    switch(mode)
    {
        case 0:  //muon plots
            //bgZJ.push_back(HnuPlots::FileStruct("DD Z+Jets",   data_mm,  "hNu/"        + cutlevels[11]    + "/" + plot,     1.0,     1.0,      3.44921630331921220e-02,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 4000, 100, &ll, &ul));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa", mc_SZJ,    "hNuMu40/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    0.719638*k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    k_mm_Zscale / NZJ,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 4000, 100, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,   k_mm_ddtop,                                     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    - k_mm_ddtop / NZJ,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 4000, 100, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuEMu/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    -k_mm_ddtop*0.916553*k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     - k_mm_ddtop / NtW,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  - k_mm_ddtop / NtbarW,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     - k_mm_ddtop / NZZ,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     - k_mm_ddtop / NWZ,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     - k_mm_ddtop / NWW,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.24820,     -k_mm_ddtop / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.08888 ,     -k_mm_ddtop / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW,                                   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    0.5 / NtW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 0.5 / NtbarW,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

//            bgZ1J.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_CWR/WW2Jets_EW6_TuneZ2star_8TeV-phantom-tauola.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.0993,     1.0 / 496500,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.24820,     1.0 / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.08888 ,     1.0 / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            //bgTT.push_back(HnuPlots::FileStruct(   "Data Fit",   "mudatafit_gt_800.root",   "mWR", 1.0, 1.0, 1.0,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false, false, 0, 0, -1, &ll, &ul));

            if(isZJ)
            {
                bg.push_back(bgZJ);
                bg.push_back(bgTT);
                bg.push_back(bgOther);
            }
            else
            {
                bg.push_back(bgTT);
                bg.push_back(bgZJ);
                bg.push_back(bgOther);
            }
            //bg.push_back(bgZ1J);
            //bg.push_back(bgQCD);

        case 7: //Limit input plots
            sprintf(fdata, "%s", data_mm.c_str());
            lumi += lumi2012mm;
            data.histpath = "hNu/" + cutlevels[cutlevel] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;
            break;
        case 1:  // electron plots
            //bgZJ.push_back(HnuPlots::FileStruct(   "DD Z+Jets",   data_ee,  "hNuE/"       + cutlevels[11]          + "/" + plot,     1.0,     1.0,      3.64878416472316086e-02,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa", mc_SZJ,    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    0.688348*k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets",   mc_Z0J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_Zscale * k_ee_NNLOZ / NZ0J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   k_ee_Zscale * k_ee_NNLOZ / NZ1J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   k_ee_Zscale * k_ee_NNLOZ / NZ2J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   k_ee_Zscale * k_ee_NNLOZ / NZ3J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   k_ee_Zscale * k_ee_NNLOZ / NZ4J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecttbar, k_top / Nttbar,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,        1.0,     1.0,   k_ee_ddtop,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuEMu/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    - 0.89016*k_ee_Zscale*k_ee_ddtop / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    - k_ee_ddtop / NZJ,    "", 0.0, 0.0 , true, 1, true, 0.0, 0.0, lt));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - k_ee_ddtop / NtW,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - k_ee_ddtop / NtbarW,                           "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - k_ee_ddtop / NZZ,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - k_ee_ddtop / NWZ,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.24820,     -k_mm_ddtop / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.08888 ,     -k_mm_ddtop / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,     0.5 / NtW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW,  0.5 / NtbarW,                                   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

//            bgZ1J.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_CWR/WW2Jets_EW6_TuneZ2star_8TeV-phantom-tauola.root",    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.0993,     1.0 / 496500,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.24820,     1.0 / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.08888 ,     1.0 / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            //bgTT.push_back(HnuPlots::FileStruct(   "Data Fit",   "edatafit_gt_800_lt_1800.root",   "mWR", 1.0, 1.0, 1.0,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false, false, 0, 0, -1, &ll, &ul));

            if(isZJ)
            {
                bg.push_back(bgZJ);
                bg.push_back(bgTT);
                bg.push_back(bgOther);
            }
            else
            {
                bg.push_back(bgTT);
                bg.push_back(bgZJ);
                bg.push_back(bgOther);
            }
            //bg.push_back(bgZ1J);
            //bg.push_back(bgQCD);

            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuE/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));

        case 8:  //Limit input plots
            sprintf(fdata, "%s", data_ee.c_str());
            lumi += lumi2012ee;
            data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;
            break;
        case 2:  // emu plots
            kf = 1/k_mm_ddtop;
            bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, kf * k_top / Nttbar, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bg.push_back(bgTT);
        case 3:
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets",   mc_Z0J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   k_mm_Zscale * k_mm_NNLOZ / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   k_mm_Zscale * k_mm_NNLOZ / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   k_mm_Zscale * k_mm_NNLOZ / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   k_mm_Zscale * k_mm_NNLOZ / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   k_mm_Zscale * k_mm_NNLOZ / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    1.0 / NtW,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 1.0 / NtbarW, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    1.0 / NZZ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    1.0 / NWZ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    1.0 / NWW,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            //bgZ1J.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_CWR/WW2Jets_EW6_TuneZ2star_8TeV-phantom-tauola.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.0993,     1.0 / 496500,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.24820,     1.0 / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.08888 ,     1.0 / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            sprintf(fdata, "%s", data_em.c_str());
            lumi += lumi2012mm;
            data.histpath = "hNuEMu/" + cutlevels[cutlevel] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;

            bg.push_back(bgZJ);
            bg.push_back(bgOther);
            bg.push_back(bgZ1J);
            break;
        case 4:
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J /* 3.44921630331921220e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J /* 3.44921630331921220e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J /* 3.44921630331921220e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J /* 3.44921630331921220e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J /* 3.44921630331921220e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,   k_mm_ddtop,                                     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     - k_mm_ddtop / NtW,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  - k_mm_ddtop / NtbarW,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     - k_mm_ddtop / NZZ,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     - k_mm_ddtop / NWZ,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     - k_mm_ddtop / NWW,                             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW,                                   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            bg.push_back(bgZJ);
            bg.push_back(bgTT);
            bg.push_back(bgOther);

            sprintf(fdata, "%s", data_mm.c_str());
            lumi += lumi2012mm;
            data.histpath = "hNu/" + cutlevels[11] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;
            break;
        case 5:
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets",   mc_Z0J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,    k_ee_Zscale * k_ee_NNLOZ / NZ0J /* 3.64878416472316086e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,    k_ee_Zscale * k_ee_NNLOZ / NZ1J /* 3.64878416472316086e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,    k_ee_Zscale * k_ee_NNLOZ / NZ2J /* 3.64878416472316086e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,    k_ee_Zscale * k_ee_NNLOZ / NZ3J /* 3.64878416472316086e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,    k_ee_Zscale * k_ee_NNLOZ / NZ4J /* 3.64878416472316086e-02*/,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));//bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecttbar, k_top / Nttbar,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot,        1.0,     1.0,   k_ee_ddtop,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - k_ee_ddtop / NtW,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - k_ee_ddtop / NtbarW,                           "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - k_ee_ddtop / NZZ,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - k_ee_ddtop / NWZ,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"  + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            bg.push_back(bgZJ);
            bg.push_back(bgTT);
            bg.push_back(bgOther);

            sprintf(fdata, "%s", data_ee.c_str());
            lumi += lumi2012ee;
            data.histpath = "hNuE/" + cutlevels[11] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;
            break;
        case 6: //Z plots
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    k_mm_Zscale / NZJ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "DY#lower[-0.20]{+}Jets Madgraph 0#lower[-0.20]{+}1 Jet" ,   mc_Z0J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ0J/xsecZJ, k_mm_NNLOZ * k_mm_Zscale / NZ0J / (k_ee_Zscale / NZJ) / 0.719638, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 1 Jet"   ,   mc_Z1J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ1J/xsecZJ, k_mm_NNLOZ * k_mm_Zscale / NZ1J / (k_ee_Zscale / NZJ) / 0.719638, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ2J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 2 Jets" ,   mc_Z2J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ2J/xsecZJ, k_mm_NNLOZ * k_mm_Zscale / NZ2J / (k_ee_Zscale / NZJ) / 0.719638, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ3J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 3 Jets" ,   mc_Z3J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ3J/xsecZJ, k_mm_NNLOZ * k_mm_Zscale / NZ3J / (k_ee_Zscale / NZJ) / 0.719638, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ4J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 4+ Jets",   mc_Z4J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ4J/xsecZJ, k_mm_NNLOZ * k_mm_Zscale / NZ4J / (k_ee_Zscale / NZJ) / 0.719638, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuMu40/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
            sprintf(fdata, "%s", mc_SZJ.c_str());
            lumi += lumi2012mm;
            data.label = "Z+Jets Sherpa";
            data.histpath = "hNuMu40/" + cutlevels[cutlevel] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;

            bg.push_back(bgZJ);
            //bg.push_back(bgZ1J);
            bg.push_back(bgZ2J);
            bg.push_back(bgZ3J);
            bg.push_back(bgZ4J);
            //bg.push_back(bgQCD);
            break;
        case 9: //Z plots
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    k_mm_Zscale / NZJ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "DY+Jets Madgraph 0+1 Jet" ,   mc_Z0J,   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ0J/xsecZJ, k_ee_NNLOZ * k_ee_Zscale / NZ0J / (k_ee_Zscale / NZJ) / 0.688348, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 1 Jet"   ,   mc_Z1J,   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ1J/xsecZJ, k_ee_NNLOZ * k_ee_Zscale / NZ1J / (k_ee_Zscale / NZJ) / 0.688348, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ2J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 2 Jets" ,   mc_Z2J,   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ2J/xsecZJ, k_ee_NNLOZ * k_ee_Zscale / NZ2J / (k_ee_Zscale / NZJ) / 0.688348, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ3J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 3 Jets" ,   mc_Z3J,   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ3J/xsecZJ, k_ee_NNLOZ * k_ee_Zscale / NZ3J / (k_ee_Zscale / NZJ) / 0.688348, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZ4J.push_back(HnuPlots::FileStruct(   "Z+Jets Madgraph 4+ Jets",   mc_Z4J,   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, xsecZ4J/xsecZJ, k_ee_NNLOZ * k_ee_Zscale / NZ4J / (k_ee_Zscale / NZJ) / 0.688348, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuMu40/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
            sprintf(fdata, "%s", mc_SZJ.c_str());
            lumi += lumi2012mm;
            data.label = "Z+Jets Sherpa";
            data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
            data.loadtuple = lt;
            data.histFromTuple = hft;
            data.thll = 0;
            data.thul = 0;
            data.thb = -1;
            data.tpll = &ll;
            data.tpul = &ul;

            bg.push_back(bgZJ);
            //bg.push_back(bgZ1J);
            bg.push_back(bgZ2J);
            bg.push_back(bgZ3J);
            bg.push_back(bgZ4J);
            //bg.push_back(bgQCD);
            break;
    }

    //bg.push_back(bgOther2);

    //data
    data.label = "Data";
    data.file = fdata;
}

void makeCutString(const int cutlevel, std::string plotname, std::string& retStr)
{

    char cutString[256];
    // record base cut level
    sprintf(cutString, "_cut%d", cutlevel);

    bool inv = false;
    if(plotname.find('!') < plotname.size())
    {
        inv = true;
        plotname.erase(plotname.find('!'), 1);
    }

    size_t cutStart = plotname.find(";");
    if(cutStart != size_t(-1) && inv) sprintf(cutString, "%s_not", cutString);
    // read cuts to implament
    for(size_t pos = cutStart + 1, npos = 0; npos != size_t(-1);pos = npos + 1)
    {
        npos = plotname.find(';', pos + 1);
        std::string tmp2 = plotname.substr(pos, npos - pos);
        size_t sepPos = 0;
        std::string cutType;
        if     ((sepPos = tmp2.find('>')) != size_t(-1)) cutType = "gt";
        else if((sepPos = tmp2.find('<')) != size_t(-1)) cutType = "lt";
        else continue;
        std::string t1 = tmp2.substr(0, sepPos), t2 = tmp2.substr(sepPos + 1, (size_t(-1)));
        t1.erase(remove(t1.begin(),t1.end(),' '),t1.end());
        t2.erase(remove(t2.begin(),t2.end(),' '),t2.end());
        sprintf(cutString, "%s_%s_%s_%s", cutString, t1.c_str(), cutType.c_str(), t2.c_str());
    }

    retStr = std::string(cutString);
}

void plot2012(int mode = 0, int cutlevel = 5, std::string plot = "mWR", int rebin = 5, bool log = true, double xmin = 0.0, double xmax = 3500.0, bool autoY = true)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    // if asked get hist from tuple, or if required
    bool hft = true;
    if(plot.find(';') != size_t(-1)) hft = true;
    else if(plot.find(':') != size_t(-1)) hft = true;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig, vsig2, vsig3, vsig4, vsig5, vsig6;
    setBgandData(mode, data, bg, lumi, cutlevel, plot, true, hft);

    std::cout << "Lumi:" << lumi << std::endl;

    std::string histograms = "", normhist = "";
    int signormbin = 0;

    switch(mode)
    {
        case 0:
            histograms = "hNuMu40/" + cutlevels[cutlevel] + "/" + plot;
            normhist = "hNuMu40/mc_type";
            signormbin = 3;
            break;
        case 1:
            histograms = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
            normhist = "hNuE/mc_type";
            signormbin = 2;
            break;
        case 2:
            histograms = "hNuEMu/" + cutlevels[cutlevel] + "/" + plot;
            normhist = "hNuEMu/mc_type";
            signormbin = 4;
            break;
    }

    //signal
    if(mode <= 1)
    {
        //Nominal signal points -- Sean comment these out before you plot more signal
        //if(rebin > 0)
        //{
        //    vsig.push_back(HnuPlots::FileStruct("#lower[0.31]{#splitline{M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.5 TeV}{M_{N} = M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}/2}}",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 0.002286, 1.140, normhist, 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //}
        //else
        //{
        //    vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.5 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 0.002286, 1.140, normhist, 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //    vsig2.push_back(HnuPlots::FileStruct("#lower[0.31]{#splitline{M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.5 TeV unbinned}{M_{N} = M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}/2}}",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 0.002286, 1.140, normhist, 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul, true));
        //}
        //sig.push_back(vsig);
        //if(rebin <= 0) sig.push_back(vsig2);
        
        //sample gen signal point -- Sean add individual signal points here
        //Format  (modify stared fields)     label*           filepath*                                                                                          tupple folder / plotname  lumi  xsec*   kfactor/Nevts*  the rest is a magic incantation that should not be changed
        vsig.push_back( HnuPlots::FileStruct("test signal"  ,  "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_8/src/HeavyNu/AnalysisModules/HeavyNu_accept_1000_25.root", "hNuGen2012/" + plot, lumi, 0.002286, 1.140/1000, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //vsig2.push_back(HnuPlots::FileStruct("test signal 2",  "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_8/src/HeavyNu/AnalysisModules/HeavyNu_accept_1000_250.root", "hNuGen2012/" + plot, lumi, 0.002286, 1.140/1000, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        
        //Then add the individual signal points to the list of signal points
        sig.push_back(vsig);
        //sig.push_back(vsig2);
    }
    else if(!hft && mode == 2)
    {
        ////vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 1.1 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 0.013339, 1.214 * 0.5 * 0.75, normhist, 0.0, 0.0, true, signormbin));
        vsig2.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 1.0 TeV",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",    "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.667875, 1.340 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        vsig3.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 1.5 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",    "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.082688, 1.293 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        vsig4.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 2.0 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.013339, 1.214 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        vsig5.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 2.5 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.002286, 1.140 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        vsig6.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 3.0 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.000393, 1.151 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));

        ////sig.push_back(vsig);
        //sig.push_back(vsig2);
        //sig.push_back(vsig3);
        //sig.push_back(vsig4);
        //sig.push_back(vsig5);
        //sig.push_back(vsig6);
    }

    HnuPlots hps(data, bg, sig, lumi);
    std::string clstring;
    makeCutString(cutlevel, plot, clstring);
    switch(mode)
    {
        case 7:
            rebin = -1;
            xmax = 4000;
        case 0:
            if(rebin > 0) hps.setFormLabel("hNu_mm_2012" + clstring);
            else          hps.setFormLabel("hNu_mm_ls_2012" + clstring);
            hps.setSavePlots(true);
            if(!plot.compare("mWR") || !plot.compare("mWR;"))
            {
                if(cutlevel == 5) hps.setYRange(0.06, 3000);
                hps.loadSystFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/systematicsdb_mu_2012.csv", "/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/ratesdb.csv", (mode == 7));
                hps.mcBgShape();
            }
            break;
        case 8:
            rebin = -1;
            xmax = 4000;
        case 1:
            if(rebin > 0) hps.setFormLabel("hNu_ee_2012" + clstring);
            else          hps.setFormLabel("hNu_ee_ls_2012" + clstring);
            hps.setSavePlots(true);
            if(!plot.compare("mWR") || !plot.compare("mWR;"))
            {
                if(cutlevel == 5) hps.setYRange(0.06, 3000);
                hps.loadSystFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/systematicsdb_elec_2012.csv", "/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/ratesdb_elec.csv", (mode == 8));
                hps.mcBgShape();
            }
            break;
        case 2:
            hps.setFormLabel("hNu_em_2012" + clstring);
            hps.setSavePlots(true);
            break;
        case 3:
            hps.setFormLabel("hNu_em2_2012" + clstring);
            hps.setSavePlots(true);
            hps.setCompPlot(false);
            break;
        case 4:
            hps.setFormLabel("ddZ_mm_2012" + clstring);
            hps.setSavePlots(true);
            break;
        case 5:
            hps.setFormLabel("ddZ_ee_2012" + clstring);
            hps.setSavePlots(true);
            break;
        case 6:
            hps.setFormLabel("sherpaZvsMadgraph_mm_2012" + clstring);
            hps.setSavePlots(true);
            break;
        case 9:
            hps.setFormLabel("sherpaZvsMadgraph_ee_2012" + clstring);
            hps.setSavePlots(true);
            break;
        default:
            hps.setSavePlots(false);
            break;
    }
    if(autoY) hps.setYAxisTitle("please auto set the axis");
    else hps.setYAxisTitle("Events");
    hps.autoSetHistogramAxisTitle(mode);
    hps.setRebin(rebin);
    hps.setLog(log);
    hps.setCompPlot(true);
    hps.setXRange(xmin, xmax);
    hps.plot();
    hps.integrals(600, 60000);
}

void plotRatios(int mode = 0, int cutlevel = 5, std::string plot = "mWR", int rebin = 5)
{
    std::vector<std::vector<HnuPlots::FileStruct> > bg, sig;
    HnuPlots::FileStruct data;
    char fdata[128];

    std::vector<HnuPlots::FileStruct> bg1, bg2, bg3, bg4;

    switch(mode)
    {
        case 0:
            bg1.push_back(HnuPlots::FileStruct(   "DD t#bar{t}",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,   k_mm_ddtop));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     - k_mm_ddtop / NtW));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  - k_mm_ddtop / NtbarW));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     - k_mm_ddtop / NZZ));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     - k_mm_ddtop / NWZ));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     - k_mm_ddtop / NWW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW));

            bg2.push_back(HnuPlots::FileStruct(   "MC t#bar{t}",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J));
            bg2.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW));

            bg3.push_back(HnuPlots::FileStruct("DD t#bar{t} fitted",    "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_6_patch1/src/HeavyNu/Tools/muonfittedDDtop.root",    "fitted", 1.0, 1.0 , 1.0));

            bg4.push_back(HnuPlots::FileStruct("MC t#bar{t} fitted",    "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_6_patch1/src/HeavyNu/Tools/muonfittedMCtop.root",    "fitted", 1.0, 1.0 , 1.0));

            bg.push_back(bg1);
            bg.push_back(bg2);
            bg.push_back(bg3);
            bg.push_back(bg4);

            sprintf(fdata, "%s", data_mm.c_str());
            //lumi += lumi2012mm;
            data.histpath = "hNu/" + cutlevels[cutlevel] + "/" + plot;
            data.label = "Data";
            data.file = fdata;
            break;
        case 1:
            bg1.push_back(HnuPlots::FileStruct(   "DD t#bar{t}",   mc_Z0J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,    k_ee_Zscale * k_ee_NNLOZ / NZ0J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,    k_ee_Zscale * k_ee_NNLOZ / NZ1J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,    k_ee_Zscale * k_ee_NNLOZ / NZ2J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,    k_ee_Zscale * k_ee_NNLOZ / NZ3J));
            bg1.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,    k_ee_Zscale * k_ee_NNLOZ / NZ4J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,   k_ee_ddtop));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,    - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,    - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,    - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,    - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,    - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,     - k_ee_ddtop / NtW));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW,  - k_ee_ddtop / NtbarW));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,     - k_ee_ddtop / NZZ));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,     - k_ee_ddtop / NWZ));
            bg1.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,     - k_ee_ddtop / NWW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,     1.0 / NtW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW,  1.0 / NtbarW));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,     1.0 / NZZ));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,     1.0 / NWZ));
            bg1.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,     1.0 / NWW));

            bg2.push_back(HnuPlots::FileStruct(   "MC t#bar{t}",   mc_Z0J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,    k_ee_Zscale * k_ee_NNLOZ / NZ0J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,    k_ee_Zscale * k_ee_NNLOZ / NZ1J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,    k_ee_Zscale * k_ee_NNLOZ / NZ2J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,    k_ee_Zscale * k_ee_NNLOZ / NZ3J));
            bg2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,    k_ee_Zscale * k_ee_NNLOZ / NZ4J));
            bg2.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecttbar, k_top / Nttbar));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,     1.0 / NtW));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW,  1.0 / NtbarW));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,     1.0 / NZZ));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,     1.0 / NWZ));
            bg2.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,     1.0 / NWW));

            bg3.push_back(HnuPlots::FileStruct("DD t#bar{t} fitted",    "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_6_patch1/src/HeavyNu/Tools/elecfittedDDtop.root",    "fitted", 1.0, 1.0 , 1.0));

            bg4.push_back(HnuPlots::FileStruct("MC t#bar{t} fitted",    "/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_6_patch1/src/HeavyNu/Tools/elecfittedMCtop.root",    "fitted", 1.0, 1.0 , 1.0));

            bg.push_back(bg1);
            bg.push_back(bg2);
            bg.push_back(bg3);
            bg.push_back(bg4);

            sprintf(fdata, "%s", data_ee.c_str());
            //lumi += lumi2012mm;
            data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
            data.label = "Data";
            data.file = fdata;
            break;
        default:
            double tlumi = 19700;
            setBgandData(mode, data, bg, tlumi);
            break;
    }

    HnuPlots hps(data, bg, sig, 19700);
    hps.setYAxisTitle("please auto set the axis");
    hps.setYAxisTitle("Data / Background");
    hps.autoSetHistogramAxisTitle(mode);
    hps.setRebin(rebin);
    hps.plotRatios();
}

void plotMCFits(int mode = 0, int cutlevel = 5, bool log = true)
{
    using namespace std;

    char plot[] = "mWR";
    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    setBgandData(mode, data, bg, lumi, cutlevel, plot);

    std::string bgestfname;
    //if(is2011A) bgestfname += "2011A.txt";
    //else  bgestfname += "2011B.txt";

    HnuPlots hps(data, bg, sig, lumi);
    hps.autoSetHistogramAxisTitle(mode);
    hps.setRebin(5);
    hps.setLog(log);
    switch(mode)
    {

        case 0:
            hps.setFormLabel("bgFits_mm_2012");
            bgestfname = "bgest_mm.txt";
            break;
        case 1:
            hps.setFormLabel("bgFits_ee_2012");
            bgestfname = "bgest_ee.txt";
            break;
        case 2:
            hps.setFormLabel("bgFits_em_2012");
            bgestfname = "bgest_em.txt";
            break;
    }
    hps.setXRange(0, 2500);
    hps.plotMCShape(bgestfname);
}

void plotMCShapes(int mode = 0, int cutlevel = 5)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    setBgandData(mode, data, bg, lumi, cutlevel, "mWR", true);

    HnuPlots hps(data, bg, sig, lumi);
    switch(mode)
    {

        case 0:
            hps.mcBgShape(cutlevel, "_mm");
            break;
        case 1:
            hps.mcBgShape(cutlevel, "_ee");
            break;
    }

}

void plotDDZJNorm(bool isMuon = true, int cutlevel = 4, bool log = true)//, std::string sample = "")
{
    using namespace std;

    char plot[] = "mWR", fdata[256];
    double lumi = 0.0;
    string datahistname;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther, bgTT;
    if(isMuon)
    {
        bgZJ.push_back(HnuPlots::FileStruct("DD Z+Jets",   data_mm,  "hNu/"        + cutlevels[11]    + "/" + plot,     1.0,     1.0,      1.0,                 ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsecZJ,    - 1.0 / NZJ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_tW,    "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsectW,    - 1.0 / NtW,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_tbarW, "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsectbarW, - 1.0 / NtbarW, ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZZ,    "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsecZZ,    - 1.0 / NZZ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_WZ,    "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsecWZ,    - 1.0 / NWZ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_WW,    "hNuMu40/"    + cutlevels[11]    + "/" + plot, lumi2012mm, xsecWW,    - 1.0 / NWW,    ""));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,   k_mm_ddtop,                        "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    - k_mm_ddtop / NZJ,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     - k_mm_ddtop / NtW,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  - k_mm_ddtop / NtbarW,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     - k_mm_ddtop / NZZ,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     - k_mm_ddtop / NWZ,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     - k_mm_ddtop / NWW,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tt,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar,                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW,                         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ,                         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ,                         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW,                         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("QCD",      data_mm,  "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, 1.0,     1.0,           1.0       ,                     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    -1.0 / NZJ,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectW,    -1.0 / NtW,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectbarW, -1.0 / NtbarW,                   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZZ,    -1.0 / NZZ,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0 / NWZ,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWW,    -1.0 / NWW,                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0,          "hNuMu1QCD/cutlevel", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    data_mm,  "hNuFakeMuGoodEwgtMu/" + cutlevelsTop[cutlevel] + "/" + plot, 1.0,     1.0, -k_mm_ddtop * 0.11,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        sprintf(fdata, "%s", data_mm.c_str());
        lumi += lumi2012mm;
        datahistname = "hNu/" + cutlevels[cutlevel] + "/" + plot;
    }
    else
    {
        bgZJ.push_back(HnuPlots::FileStruct(   "DD Z+Jets",   data_ee,  "hNuE/"       + cutlevels[11]          + "/" + plot,     1.0,     1.0,      1.0,                 ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsecZJ,    - 1.0 / NZJ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_tW,    "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsectW,    - 1.0 / NtW,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_tbarW, "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsectbarW, - 1.0 / NtbarW, ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZZ,    "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsecZZ,    - 1.0 / NZZ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_WZ,    "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsecWZ,    - 1.0 / NWZ,    ""));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_WW,    "hNuE/"       + cutlevels[11]          + "/" + plot, lumi2012mm, xsecWW,    - 1.0 / NWW,    ""));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,        1.0,     1.0,   k_ee_ddtop,                         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevelsTop[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    - k_ee_ddtop / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - k_ee_ddtop / NtW,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - k_ee_ddtop / NtbarW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - k_ee_ddtop / NZZ,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - k_ee_ddtop / NWZ,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,                          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

        sprintf(fdata, "%s", data_ee.c_str());
        lumi += lumi2012ee;
        datahistname = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
    }
    bg.push_back(bgZJ);
    bg.push_back(bgTT);
    bg.push_back(bgOther);

    //data
    HnuPlots::FileStruct data("Data", fdata, datahistname);

    HnuPlots hps(data, bg, sig, lumi);
    hps.autoSetHistogramAxisTitle(!isMuon);
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(5);
    hps.setLog(log);
    std::string flabel = "ddZJnorm_";
    std::string clstring;
    makeCutString(cutlevel, "", clstring);
    if(isMuon) hps.setFormLabel(flabel + "_mm_2012" + clstring);
    else hps.setFormLabel(flabel + "_ee_2012" + clstring);
    hps.setSavePlots(true);
    hps.setXRange(0.0, 2500.0);
    hps.scaleByShape(600, 2500, 1);
    //hps.plot();
    //hps.plotNorm(120, 200);
}

void plotDDTTNorm(bool isMuon = true, int nb = 1, int cutlevel = 4, bool log = true)//, std::string sample = "")
{
    using namespace std;

    char plot[16] = "mLL", fdata[256];
    double lumi = 0.0;
    string datahistname;

    if(nb == 1) sprintf(plot, "mLL_1b");
    else if(nb == 2) sprintf(plot, "mLL_2b");

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther, bgTT;

    //data
    HnuPlots::FileStruct data;

    if(isMuon)
    {
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    k_mm_Zscale / NZJ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   k_mm_Zscale * k_mm_NNLOZ / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   k_mm_Zscale * k_mm_NNLOZ / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   k_mm_Zscale * k_mm_NNLOZ / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   k_mm_Zscale * k_mm_NNLOZ / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   k_mm_Zscale * k_mm_NNLOZ / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,  1.0,                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    - k_mm_ddtop / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   - k_mm_Zscale * k_mm_NNLOZ / NZ0J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   - k_mm_Zscale * k_mm_NNLOZ / NZ1J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   - k_mm_Zscale * k_mm_NNLOZ / NZ2J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   - k_mm_Zscale * k_mm_NNLOZ / NZ3J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   - k_mm_Zscale * k_mm_NNLOZ / NZ4J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    - 1.0 / NtW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, - 1.0 / NtbarW,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    - 1.0 / NZZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    - 1.0 / NWZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    - 1.0 / NWW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tt,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    1.0 / NtW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 1.0 / NtbarW,           "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    1.0 / NZZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    1.0 / NWZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    1.0 / NWW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuMu40/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", data_mm.c_str());
        lumi += lumi2012mm;
        data.histpath = "hNu/" + cutlevels[cutlevel] + "/" + plot;

        bg.push_back(bgTT);
        bg.push_back(bgZJ);
        bg.push_back(bgOther);
        //bg.push_back(bgQCD);
    }
    else
    {
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_Zscale * k_ee_NNLOZ / NZ0J,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   k_ee_Zscale * k_ee_NNLOZ / NZ1J,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   k_ee_Zscale * k_ee_NNLOZ / NZ2J,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   k_ee_Zscale * k_ee_NNLOZ / NZ3J,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   k_ee_Zscale * k_ee_NNLOZ / NZ4J,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,       1.0,     1.0,    1.0,                     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    - k_ee_ddtop / NZJ,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_Zscale * k_ee_NNLOZ / NZ0J,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_Zscale * k_ee_NNLOZ / NZ1J,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_Zscale * k_ee_NNLOZ / NZ2J,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_Zscale * k_ee_NNLOZ / NZ3J,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_Zscale * k_ee_NNLOZ / NZ4J,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - 1.0 / NtW,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - 1.0 / NtbarW,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - 1.0 / NZZ,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - 1.0 / NWZ,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - 1.0 / NWW,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuE/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));

        sprintf(fdata, "%s", data_ee.c_str());
        lumi += lumi2012ee;
        data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;

        bg.push_back(bgTT);
        bg.push_back(bgZJ);
        bg.push_back(bgOther);
        //bg.push_back(bgQCD);
    }

    data.label = "Data";
    data.file = fdata;

    HnuPlots hps(data, bg, sig, lumi);
    hps.autoSetHistogramAxisTitle(!isMuon);
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(1);
    hps.setLog(log);
    std::string flabel = "ddtopnorm_";
    std::string flend = "";
    if(nb == 1) flend = "_1b";
    else if(nb == 2) flend = "_2b";
    std::string clstring;
    makeCutString(cutlevel, "", clstring);
    if(isMuon) hps.setFormLabel(flabel + "mm_2012" + flend + clstring);
    else       hps.setFormLabel(flabel + "ee_2012" + flend + clstring);
    hps.setSavePlots(true);
    hps.setXRange(60.0, 500.0);
    hps.scaleByShape(120, 200, 2);
    //hps.plotNorm(120, 200);
}

void plotTTBarNorm(int cutlevel = 5, bool log = true)
{

    using namespace std;

    char plot[] = "mWR";
    double lumi = 0.0;
    std::string datahistname;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    //data
    HnuPlots::FileStruct data;

    setBgandData(2, data, bg, lumi, cutlevel, plot, true);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("M_{e#mujj} [GeV]");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(4);
    hps.setLog(log);
    hps.setFormLabel("ttnorm_mm_2012_MC");
    hps.setSavePlots(false);
    //hps.plotNorm(20.0, 6000.0);
    hps.scaleByShape(20.0, 6000.0, 1);
}

void plotAllNorm(int mode = 0, int cutlevel = 5, bool log = true)
{

    using namespace std;

    char plot[] = "mWR";
    double lumi = 0.0;
    std::string datahistname;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    //data
    HnuPlots::FileStruct data;

    setBgandData(mode, data, bg, lumi, cutlevel, plot, true);

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("M_{e#mujj} [GeV]");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(5);
    hps.setLog(log);
    hps.setFormLabel("ttnorm_mm_2012_MC");
    hps.setSavePlots(false);
    //hps.plotNorm(20.0, 6000.0);
    hps.scaleByShape(20.0, 6000.0, 3);
}

void plotTTBarDDEffBasedNorm(int cutlevel = 4, const double xmin = 120.0, const double xmax = 200.0)//, std::string sample = "")
{
    using namespace std;

    char plot[] = "mLL", fdata[256], fdata2[256];
    double lumi = 0.0, lumi2 = 0.0;
    std::string datahistname, datahistname2;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, bg2, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther, bgZJ2, bgOther2;

    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J,   ""));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J,   ""));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J,   ""));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J,   ""));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J,   ""));
    bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW,                         ""));
    bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW,                      ""));
    bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ,                         ""));
    bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ,                         ""));
    bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW,                         ""));
    sprintf(fdata, "%s", data_mm.c_str());
    lumi += lumi2012mm;
    datahistname = "hNu/" + cutlevels[cutlevel] + "/" + plot;

    bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_Zscale * k_mm_NNLOZ / NZ0J,     ""));
    bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   k_ee_Zscale * k_mm_NNLOZ / NZ1J,     ""));
    bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   k_ee_Zscale * k_mm_NNLOZ / NZ2J,     ""));
    bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   k_ee_Zscale * k_mm_NNLOZ / NZ3J,     ""));
    bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   k_ee_Zscale * k_mm_NNLOZ / NZ4J,     ""));
    bgOther2.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,                           ""));
    bgOther2.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,                        ""));
    bgOther2.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,                           ""));
    bgOther2.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,                           ""));
    bgOther2.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"   + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,                           ""));
    sprintf(fdata2, "%s", data_ee.c_str());
    lumi2 += lumi2012ee;
    datahistname2 = "hNuE/" + cutlevels[cutlevel] + "/" + plot;

    bg.push_back(bgZJ);
    bg.push_back(bgOther);
    bg2.push_back(bgZJ2);
    bg2.push_back(bgOther2);

    //data

    double int_mm = 0, int_ee = 0, err_mm = 0, err_ee = 0;

    if(true)
    {
        HnuPlots::FileStruct data("Data", fdata, datahistname);
        HnuPlots hps(data, bg, sig, lumi);
        hps.integrals(xmin, xmax, &int_mm, &err_mm);
    }
    if(true)
    {

        HnuPlots::FileStruct data2("Data", fdata2, datahistname2);
        HnuPlots hps2(data2, bg2, sig, lumi2);
        hps2.integrals(xmin, xmax, &int_ee, &err_ee);
    }

    //double mu_eff = 0.901;
    double mu_eff = 0.871660;

    int_ee *= lumi / lumi2;

    std::cout << err_ee << "\t" << err_mm << std::endl;
    std::cout << "Muon dd top scale factor: " << 0.5 * sqrt(int_mm / int_ee) * (2 - mu_eff) << " +/- "           << (1.0 / 4) * sqrt(int_mm / int_ee) * sqrt(pow(err_mm / int_mm, 2) + pow(err_ee / int_ee, 2)) * (2 - mu_eff) << std::endl;
    std::cout << "Elec dd top scale factor: " << 0.5 * sqrt(int_ee / int_mm) * lumi2 / lumi / mu_eff  << " +/- " << (1.0 / 4) * sqrt(int_ee / int_mm) * sqrt(pow(err_mm / int_mm, 2) + pow(err_ee / int_ee, 2)) * lumi2 / lumi / mu_eff << std::endl;
}

void plotTTBarMCNorm(bool isMuon = true, int cutlevel = 5, int nb = 0, bool log = true)//, std::string sample = "")
{
    using namespace std;

    char plot[16], fdata[256];
    double lumi = 0.0;

    switch(nb)
    {
        case 0:
        default:
            sprintf(plot, "mWR");
            break;
        case 1:
            sprintf(plot, "mWR_1b");
            break;
        case 2:
            sprintf(plot, "mWR_2b");
            break;
    }

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgTT;

    bgTT.push_back(HnuPlots::FileStruct("t#bar{t} e#mu", mc_tt, "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, ""));
    lumi += lumi2012mm;

    bg.push_back(bgTT);

    //data
    sprintf(fdata, "%s", mc_tt.c_str());
    HnuPlots::FileStruct data("", fdata, "");

    if(isMuon)
    {
        data.histpath = "hNuMu40/" + cutlevels[cutlevel] + "/" + plot;
        data.label = "t#bar{t} #mu#mu";
    }
    else
    {
        data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;
        data.label = "t#bar{t} ee";
    }

    HnuPlots hps(data, bg, sig, lumi);
    hps.setXAxisTitle("M_{e#mujj} [GeV]");
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(2);
    hps.setLog(log);
    //if(is2011A && !is2011B) hps.setFormLabel("ttnorm_2011A");
    //else if(!is2011A && is2011B) hps.setFormLabel("ttnorm_2011B");
    std::string clstring;
    makeCutString(cutlevel, "", clstring);
    if(isMuon) hps.setFormLabel("ttMCnorm_mm_2012" + clstring);
    else       hps.setFormLabel("ttMCnorm_ee_2012" + clstring);
    hps.setSavePlots(true);
    hps.plotNorm(40.0, 6000.0);
}

void plotZJNorm(int mode = 0, int cutlevel = 4, bool log = true)//, std::string sample = "")
{
    using namespace std;

    char plot[] = "mLLNorm";//, fdata[256];
    double lumi = 0.0;
    string datahistname;
    //bool lt = true, hft = true;

    //data
    HnuPlots::FileStruct data;

    ll.cutlevel = cutlevel;
    ll.mlljj    = 0.0;
    ll.mll      = 0.0;
    ll.l1pt     = 0.0;
    ll.l1eta    = -10.0;
    ll.l1phi    = -10.0;
    ll.l2pt     = 0.0;
    ll.l2eta    = -10.0;
    ll.l2phi    = -10.0;
    ll.j1pt     = 0.0;
    ll.j1eta    = -10.0;
    ll.j1phi    = -10.0;
    ll.j2pt     = 0.0;
    ll.j2eta    = -10.0;
    ll.j2phi    = -10.0;
    ul.cutlevel = 100;
    ul.mlljj    = 8000.0;
    ul.mll      = 8000.0;
    ul.l1pt     = 40000.0;
    ul.l1eta    = 10.0;
    ul.l1phi    = 10.0;
    ul.l2pt     = 40000.0;
    ul.l2eta    = 10.0;
    ul.l2phi    = 10.0;
    ul.j1pt     = 40000.0;
    ul.j1eta    = 10.0;
    ul.j1phi    = 10.0;
    ul.j2pt     = 40000.0;
    ul.j2eta    = 10.0;
    ul.j2phi    = 10.0;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther, bgTT;
    /*if(isMuon)
    {
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa", mc_SZJ,    "hNuMu40/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J, k_mm_NNLOZ  / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J, k_mm_NNLOZ  / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J, k_mm_NNLOZ  / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J, k_mm_NNLOZ  / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J, k_mm_NNLOZ  / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,  k_mm_ddtop,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    - k_mm_ddtop / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    - k_mm_ddtop / NtW,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, - k_mm_ddtop / NtbarW,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    - k_mm_ddtop / NZZ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    - k_mm_ddtop / NWZ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    - k_mm_ddtop / NWW,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.24820,     -k_mm_ddtop / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.08888 ,     -k_mm_ddtop / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tt,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    1.0 / NtW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 1.0 / NtbarW,           "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    1.0 / NZZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    1.0 / NWZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    1.0 / NWW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.24820,     1.0 / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, 0.08888 ,     1.0 / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        //bgOther.push_back(HnuPlots::FileStruct("QCD",      data_mm,  "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, 1.0,     1.0,           1.0       ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    -1.0 / NZJ,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectW,    -1.0 / NtW,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectbarW, -1.0 / NtbarW,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZZ,    -1.0 / NZZ,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0 / NWZ,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWW,    -1.0 / NWW,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0, "hNuMu1QCD/cutlevel", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    data_mm,  "hNuFakeMuGoodEwgtMu/" + cutlevelsTop[cutlevel] + "/" + plot, 1.0,     1.0, -k_mm_ddtop * 0.11,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

        
        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuMu40/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", data_mm.c_str());
        lumi += lumi2012mm;
        data.histpath = "hNu/" + cutlevels[cutlevel] + "/" + plot;

        bg.push_back(bgZJ);
        bg.push_back(bgTT);
        bg.push_back(bgOther);

    }
    else
    {
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa", mc_SZJ,    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_NNLOZ  / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   k_ee_NNLOZ  / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   k_ee_NNLOZ  / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   k_ee_NNLOZ  / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   k_ee_NNLOZ  / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, 1.0 / Nttbar,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,        1.0,     1.0,   k_ee_ddtop,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevelsTop[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    - k_ee_ddtop / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - k_ee_ddtop / NtW,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - k_ee_ddtop / NtbarW,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - k_ee_ddtop / NZZ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - k_ee_ddtop / NWZ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.24820,     -k_mm_ddtop / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        bgTT.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuEMu/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.08888 ,     -k_mm_ddtop / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WpWpqq_8TeV.root",    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.24820,     1.0 / 99985,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        bgOther.push_back(HnuPlots::FileStruct("EWK WW 2j",    "/local/cms/user/pastika/heavyNuAnalysis_2012/WWqq/WmWmqq_8TeV.root",    "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, 0.08888 ,     1.0 / 96392,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
        //bgOther.push_back(HnuPlots::FileStruct("QCD",      data_ee,  "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, 1.0,     1.0,           1.0      ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZJ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    -1.0 / NZJ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsectW,    -1.0 / NtW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsectbarW, -1.0 / NtbarW,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZZ,    -1.0 / NZZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWZ,    -1.0 / NWZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWW,    -1.0 / NWW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WJ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWZ,    -1.0, "hNuE1QCD/cutlevel", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    data_mm,  "hNuGoodMuFakeEwgtE/" + cutlevelsTop[cutlevel] + "/" + plot, 1.0,     1.0,  -k_ee_ddtop * 0.03,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuE/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
  
        sprintf(fdata, "%s", data_ee.c_str());
        lumi += lumi2012ee;
        data.histpath = "hNuE/" + cutlevels[cutlevel] + "/" + plot;

        bg.push_back(bgZJ);
        bg.push_back(bgTT);
        bg.push_back(bgOther);
    }

    //data
    data.loadtuple = lt;
    data.histFromTuple = hft;
    data.thll = 0;
    data.thul = 0;
    data.thb = -1;
    data.tpll = &ll;
    data.tpul = &ul;

    data.label = "Data";
    data.file = fdata;*/
    
    setBgandData(mode, data, bg, lumi, cutlevel, plot, true, true, true);

    HnuPlots hps(data, bg, sig, lumi);
    hps.autoSetHistogramAxisTitle(mode);
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(10);
    hps.setLog(log);
    hps.setCompPlot(false);
    std::string clstring;
    makeCutString(cutlevel, "", clstring);
    if(mode == 0)      hps.setFormLabel("zjnorm_mm_2012" + clstring);
    else if(mode == 1) hps.setFormLabel("zjnorm_ee_2012" + clstring);
    hps.setSavePlots(true);
    hps.setXRange(0.06, 0.5);
    hps.plotNorm(0.06, 0.5);
}

void plotCutFlow(int mode = 0)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig;
    setBgandData(mode, data, bg, lumi, 0, "cutlevel", true, true);

    //signal
    std::string histograms = "", normhist = "";
    int signormbin = 0;

    switch(mode)
    {

        case 0:
            histograms = "hNuMu40/cutlevel";
            normhist = "hNuMu40/mc_type";
            signormbin = 3;
            break;
        case 1:
            histograms = "hNuE/cutlevel";
            normhist = "hNuE/mc_type";
            signormbin = 2;
            break;
        case 2:
            histograms = "hNuEMu/cutlevel";
            normhist = "hNuEMu/mc_type";
            break;
    }

    //signal
    //vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.0 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 0.013339, 1.214, normhist, 0.0, 0.0, true, signormbin));
    //vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.4 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2400_MNu-1200_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 0.003225, 1.164, normhist, 0.0, 0.0, true, signormbin));//, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
    vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.5 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 0.002286, 1.140, normhist, 0.0, 0.0, true, signormbin));//s, true, 0.0, 0.0, true, hft, 0, 0, -1, &ll, &ul));
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
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root",                   "hNuQCDMu24/" + cutlevel + "/" + plot,  216.2, 3160,  1.49 * 0.02269,  "hNuMu24/cutlevel"));
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011A_QCD_jan6.root",                   "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 3160,  1.49 * 0.02269,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_QCD_dec27.root",           "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  16.17, 1.06,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011A_QCD_dec27.root",           "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 16.17, 1.06,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root",    "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root",    "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root", "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.3,  1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011A_QCD_dec27.root", "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  43, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  18.2, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu24/" + cutlevel + "/" + plot, 216.2,  5.9, 1.0,  "hNuMu24/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011A_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-dec23.root");
        lumi += 216.2 + 1956.7;
    }
    else
    {
        bgZJ.push_back(   HnuPlots::FileStruct("Z+Jets",   "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_2011B_QCD_jan6.root",                    "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 3160, 1.48 * 0.02269 * 1.16,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6_2011B_QCD_dec27.root",              "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 16.17, 0.95,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa_2011B_QCD_jan9.root",                            "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 31314, 1.48 * 0.02269 * 1.16,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_QCD_dec27.root",    "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.3, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_2011B_QCD_dec27.root", "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.3,  1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 43, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 18.2, 1.0,  "hNuMu40/cutlevel"));
        bgOther.push_back(HnuPlots::FileStruct("Other", "/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola_2011B_QCD_dec27.root",                "hNuQCDMu40/" + cutlevel + "/" + plot, 2510.5 + 216.2 + 1956.7, 5.9, 1.0,  "hNuMu40/cutlevel"));
        sprintf(fdata, "%s", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root");
        lumi += 2510.5 + 216.2 + 1956.7;
    }
    //bg.push_back(bgTT);
    bg.push_back(bgZJ);
    bg.push_back(bgOther);

    //data
    //HnuPlots::FileStruct data("Data", "data-run2011a-run2011b-nov11.root", "hNu/" + cutlevels[cutlevel] + "/" + plot);
    HnuPlots::FileStruct data("Data", fdata, "hNuQCD/" + cutlevel + "/" + plot);

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

void plotSigEff(int mode)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data("Data", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root", "hNu/cutlevel");

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig10, vsig11, vsig12, vsig13, vsig14, vsig15, vsig16, vsig17, vsig18, vsig19, vsig20, vsig21, vsig22, vsig23, vsig24, vsig25, vsig26, vsig27, vsig28, vsig29, vsig30, vsig31, vsig32, vsig33, vsig34, vsig35, vsig36, vsig37, vsig38, vsig39, vsig40, vsig41, vsig42, vsig43;
    vector<HnuPlots::FileStruct> vsig44, vsig45, vsig46, vsig47, vsig48, vsig49, vsig50, vsig51, vsig52, vsig53, vsig54, vsig55, vsig56, vsig57, vsig58, vsig59, vsig60, vsig61, vsig62, vsig63;
    //setBgandData(true, true data, bg, lumi, 9, "cutlevel");

    std::string histograms = "", normhist = "", label = "";
    int signormbin = 0;

    switch(mode)
    {

        case 0:
            histograms = "hNuMu40/cut5_diLmass/mWR";
            normhist = "hNuMu40/mc_type";
            label = "";
            signormbin = 3;
            lumi = lumi2012mm;
            break;
        case 1:
            histograms = "hNuE/cut5_diLmass/mWR";
            normhist = "hNuE/mc_type";
            label = "";
            signormbin = 2;
            lumi = lumi2012ee;
            break;
        case 2:
            histograms = "hNuEMu/cut5_diLmass/mWR";
            normhist = "hNuEMu/mc_type";
            label = "";
            lumi = lumi2012mm;
            break;
    }

    //signal  the xsecs and k-factors here are wrong because they do not matter for eff!
    vsig10.push_back(HnuPlots::FileStruct("signal_700_350" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-700_MNu-350_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig11.push_back(HnuPlots::FileStruct("signal_800_400" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-800_MNu-400_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig12.push_back(HnuPlots::FileStruct("signal_900_450" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-900_MNu-450_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig13.push_back(HnuPlots::FileStruct("signal_1000_500" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig14.push_back(HnuPlots::FileStruct("signal_1100_550" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig15.push_back(HnuPlots::FileStruct("signal_1200_600" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1200_MNu-600_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig16.push_back(HnuPlots::FileStruct("signal_1300_650" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1300_MNu-650_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig17.push_back(HnuPlots::FileStruct("signal_1400_700" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1400_MNu-700_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig18.push_back(HnuPlots::FileStruct("signal_1500_750" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig19.push_back(HnuPlots::FileStruct("signal_1600_800" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1600_MNu-800_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig20.push_back(HnuPlots::FileStruct("signal_1700_850" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1700_MNu-850_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig21.push_back(HnuPlots::FileStruct("signal_1800_900" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1800_MNu-900_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig22.push_back(HnuPlots::FileStruct("signal_1900_950" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1900_MNu-950_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig23.push_back(HnuPlots::FileStruct("signal_2000_1000" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig24.push_back(HnuPlots::FileStruct("signal_2100_1050" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2100_MNu-1050_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig25.push_back(HnuPlots::FileStruct("signal_2200_1100" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2200_MNu-1100_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig26.push_back(HnuPlots::FileStruct("signal_2300_1150" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2300_MNu-1150_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig27.push_back(HnuPlots::FileStruct("signal_2400_1200" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2400_MNu-1200_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig28.push_back(HnuPlots::FileStruct("signal_2500_1250" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig29.push_back(HnuPlots::FileStruct("signal_2600_1300" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2600_MNu-1300_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig30.push_back(HnuPlots::FileStruct("signal_2700_1350" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2700_MNu-1350_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig31.push_back(HnuPlots::FileStruct("signal_2800_1400" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2800_MNu-1400_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig32.push_back(HnuPlots::FileStruct("signal_2900_1450" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig33.push_back(HnuPlots::FileStruct("signal_3000_1500" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig34.push_back(HnuPlots::FileStruct("signal_3100_1550" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3100_MNu-1550_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig35.push_back(HnuPlots::FileStruct("signal_3200_1600" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3200_MNu-1600_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig36.push_back(HnuPlots::FileStruct("signal_3300_1650" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3300_MNu-1650_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig37.push_back(HnuPlots::FileStruct("signal_3400_1700" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3400_MNu-1700_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig38.push_back(HnuPlots::FileStruct("signal_3500_1750" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3500_MNu-1750_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig39.push_back(HnuPlots::FileStruct("signal_3600_1800" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3600_MNu-1800_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig40.push_back(HnuPlots::FileStruct("signal_3700_1850" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3700_MNu-1850_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig41.push_back(HnuPlots::FileStruct("signal_3800_1900" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3800_MNu-1900_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig42.push_back(HnuPlots::FileStruct("signal_3900_1950" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3900_MNu-1950_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig43.push_back(HnuPlots::FileStruct("signal_4000_2000" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-4000_MNu-2000_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));

    vsig44.push_back(HnuPlots::FileStruct("signal_1000_62"   + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1000_MNu-62_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig45.push_back(HnuPlots::FileStruct("signal_1000_125"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1000_MNu-125_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig46.push_back(HnuPlots::FileStruct("signal_1000_187"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1000_MNu-187_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig47.push_back(HnuPlots::FileStruct("signal_1000_250"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1000_MNu-250_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig48.push_back(HnuPlots::FileStruct("signal_1000_833"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1000_MNu-833_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig49.push_back(HnuPlots::FileStruct("signal_1500_93"   + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1500_MNu-93_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig50.push_back(HnuPlots::FileStruct("signal_1500_187"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1500_MNu-187_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig51.push_back(HnuPlots::FileStruct("signal_1500_281"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1500_MNu-281_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig52.push_back(HnuPlots::FileStruct("signal_1500_375"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1500_MNu-375_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig53.push_back(HnuPlots::FileStruct("signal_1500_1250" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-1500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig54.push_back(HnuPlots::FileStruct("signal_2000_125"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-2000_MNu-125_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig55.push_back(HnuPlots::FileStruct("signal_2000_250"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-2000_MNu-250_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig56.push_back(HnuPlots::FileStruct("signal_2000_375"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-2000_MNu-375_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig57.push_back(HnuPlots::FileStruct("signal_2000_500"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-2000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig58.push_back(HnuPlots::FileStruct("signal_2000_1666" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-2000_MNu-1666_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig59.push_back(HnuPlots::FileStruct("signal_3000_187"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-187_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig60.push_back(HnuPlots::FileStruct("signal_3000_375"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-375_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig61.push_back(HnuPlots::FileStruct("signal_3000_562"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-562_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig62.push_back(HnuPlots::FileStruct("signal_3000_750"  + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));
    vsig63.push_back(HnuPlots::FileStruct("signal_3000_2500" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MWR-3000_MNu-2500_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi, 1.0, 1.0 / (lumi), normhist, 0.0, 0.0, true, signormbin));

    sig.push_back(vsig10);
    sig.push_back(vsig11);
    sig.push_back(vsig12);
    sig.push_back(vsig13);
    sig.push_back(vsig14);
    sig.push_back(vsig15);
    sig.push_back(vsig16);
    sig.push_back(vsig17);
    sig.push_back(vsig18);
    sig.push_back(vsig19);
    sig.push_back(vsig20);
    sig.push_back(vsig21);
    sig.push_back(vsig22);
    sig.push_back(vsig23);
    sig.push_back(vsig24);
    sig.push_back(vsig25);
    sig.push_back(vsig26);
    sig.push_back(vsig27);
    sig.push_back(vsig28);
    sig.push_back(vsig29);
    sig.push_back(vsig30);
    sig.push_back(vsig31);
    sig.push_back(vsig32);
    sig.push_back(vsig33);
    sig.push_back(vsig34);
    sig.push_back(vsig35);
    sig.push_back(vsig36);
    sig.push_back(vsig37);
    sig.push_back(vsig38);
    sig.push_back(vsig39);
    sig.push_back(vsig40);
    sig.push_back(vsig41);
    sig.push_back(vsig42);
    sig.push_back(vsig43);
    sig.push_back(vsig44);
    sig.push_back(vsig45);
    sig.push_back(vsig46);
    sig.push_back(vsig47);
    sig.push_back(vsig48);
    sig.push_back(vsig49);
    sig.push_back(vsig50);
    sig.push_back(vsig51);
    sig.push_back(vsig52);
    sig.push_back(vsig53);
    sig.push_back(vsig54);
    sig.push_back(vsig55);
    sig.push_back(vsig56);
    sig.push_back(vsig57);
    sig.push_back(vsig58);
    sig.push_back(vsig59);
    sig.push_back(vsig60);
    sig.push_back(vsig61);
    sig.push_back(vsig62);
    sig.push_back(vsig63);

    HnuPlots hps(data, bg, sig, lumi);
    hps.sigEff();
    hps.sigRMS();
    hps.sigStatErr();
}

void plotSigMatch(int mode = 0)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data("Data", "/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root", "hNu/cutlevel");

    std::string histograms = "", label = "";
    int signormbin = 0;

    switch(mode)
    {

        case 0:
            histograms = "hNuMu40/cut5_diLmass/nuLMatchedJets";
            label = "_mm";
            signormbin = 3;
            break;
        case 1:
            histograms = "hNuE/cut5_diLmass/nuLMatchedJets";
            label = "_ee";
            signormbin = 2;
            break;
        case 2:
            histograms = "hNuTopMu40/cut5_diLmass/nuLMatchedJets";
            label = "_em";
            break;
    }


    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig10, vsig11, vsig12, vsig13, vsig14, vsig15, vsig16, vsig17, vsig18, vsig19, vsig20, vsig21, vsig22, vsig23, vsig24, vsig25, vsig26, vsig27, vsig28, vsig29, vsig30, vsig31, vsig32, vsig33;
    //setBgandData(true, true data, bg, lumi, 9, "cutlevel");

    //signal  the xsecs and k-factors here are wrong because they do not matter for eff!
    vsig10.push_back(HnuPlots::FileStruct("signal_700_350" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-700_MNu-350_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig11.push_back(HnuPlots::FileStruct("signal_800_400" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-800_MNu-400_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig12.push_back(HnuPlots::FileStruct("signal_900_450" + label,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-900_MNu-450_TuneZ2star_8TeV-pythia6-tauola.root",   histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig13.push_back(HnuPlots::FileStruct("signal_1000_500" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig14.push_back(HnuPlots::FileStruct("signal_1100_550" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig15.push_back(HnuPlots::FileStruct("signal_1200_600" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1200_MNu-600_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig16.push_back(HnuPlots::FileStruct("signal_1300_650" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1300_MNu-650_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig17.push_back(HnuPlots::FileStruct("signal_1400_700" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1400_MNu-700_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig18.push_back(HnuPlots::FileStruct("signal_1500_750" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig19.push_back(HnuPlots::FileStruct("signal_1600_800" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1600_MNu-800_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig20.push_back(HnuPlots::FileStruct("signal_1700_850" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1700_MNu-850_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig21.push_back(HnuPlots::FileStruct("signal_1800_900" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1800_MNu-900_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig22.push_back(HnuPlots::FileStruct("signal_1900_950" + label,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-1900_MNu-950_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig23.push_back(HnuPlots::FileStruct("signal_2000_1000" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig24.push_back(HnuPlots::FileStruct("signal_2100_1050" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2100_MNu-1050_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig25.push_back(HnuPlots::FileStruct("signal_2200_1100" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2200_MNu-1100_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig26.push_back(HnuPlots::FileStruct("signal_2300_1150" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2300_MNu-1150_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig27.push_back(HnuPlots::FileStruct("signal_2400_1200" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2400_MNu-1200_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig28.push_back(HnuPlots::FileStruct("signal_2500_1250" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig29.push_back(HnuPlots::FileStruct("signal_2600_1300" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2600_MNu-1300_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig30.push_back(HnuPlots::FileStruct("signal_2700_1350" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2700_MNu-1350_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig31.push_back(HnuPlots::FileStruct("signal_2800_1400" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2800_MNu-1400_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig32.push_back(HnuPlots::FileStruct("signal_2900_1450" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));
    vsig33.push_back(HnuPlots::FileStruct("signal_3000_1500" + label, "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/heavynu_2012Bg_heavyNuAnalysis_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root", histograms, lumi2012mm, 1.0, 1.0 / (lumi2012mm), "hNuMu40/mc_type", 0.0, 0.0, true, signormbin));

    sig.push_back(vsig10);
    sig.push_back(vsig11);
    sig.push_back(vsig12);
    sig.push_back(vsig13);
    sig.push_back(vsig14);
    sig.push_back(vsig15);
    sig.push_back(vsig16);
    sig.push_back(vsig17);
    sig.push_back(vsig18);
    sig.push_back(vsig19);
    sig.push_back(vsig20);
    sig.push_back(vsig21);
    sig.push_back(vsig22);
    sig.push_back(vsig23);
    sig.push_back(vsig24);
    sig.push_back(vsig25);
    sig.push_back(vsig26);
    sig.push_back(vsig27);
    sig.push_back(vsig28);
    sig.push_back(vsig29);
    sig.push_back(vsig30);
    sig.push_back(vsig31);
    sig.push_back(vsig32);
    sig.push_back(vsig33);

    HnuPlots hps(data, bg, sig, lumi);
    hps.sigMatch();
}

void mcSystCalcSetBg(int mode, std::vector<std::vector<HnuPlots::FileStruct> >& bg, std::string sample, std::string uncert, int cutlevel, std::string plot)
{
    std::vector<HnuPlots::FileStruct> bgZJ, bgOther, bgZJ2, bgOther2;
    std::vector<HnuPlots::FileStruct> vsig10, vsig11, vsig12, vsig13, vsig14, vsig15, vsig16, vsig17, vsig18, vsig19, vsig20, vsig21, vsig22, vsig23, vsig24, vsig25, vsig26, vsig27, vsig28, vsig29, vsig30, vsig31, vsig32, vsig33, vsig34, vsig35, vsig36, vsig37, vsig38, vsig39, vsig40, vsig41, vsig42, vsig43;
    std::vector<HnuPlots::FileStruct> vsig102, vsig112, vsig122, vsig132, vsig142, vsig152, vsig162, vsig172, vsig182, vsig192, vsig202, vsig212, vsig222, vsig232, vsig242, vsig252, vsig262, vsig272, vsig282, vsig292, vsig302, vsig312, vsig322, vsig332, vsig342, vsig352, vsig362, vsig372, vsig382, vsig392, vsig402, vsig412, vsig422, vsig432;

    std::string amode, tmode;
    double lumi2012 = 0.0;
    switch(mode)
    {
        case 0:
            amode = "Mu40";
            lumi2012 = lumi2012mm;
            tmode = "";
            break;
        case 1:
            amode = "E";
            lumi2012 = lumi2012ee;
            tmode = "";
            break;
        case 2:
            amode = "EMu";
            lumi2012 = lumi2012mm;
            tmode = "";
            break;
    }

    //background
    bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample,   mc_Z0J, "hNu" + amode +          "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample,   mc_Z1J, "hNu" + amode +          "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample,   mc_Z2J, "hNu" + amode +          "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample,   mc_Z3J, "hNu" + amode +          "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample,   mc_Z4J, "hNu" + amode +          "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));

    bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample,   mc_Z0J, "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample,   mc_Z1J, "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample,   mc_Z2J, "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample,   mc_Z3J, "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample,   mc_Z4J, "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));

    //bgZJ.push_back(HnuPlots::FileStruct(   "zjets" + tmode + "," + sample, mc_ZJ,         "hNu" + amode + "/" +          cutlevels[5] + "/" + plot, lumi2012, xsecZJ,    k_mm_Zscale / NZJ, "hNu" + amode + "/" +          cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    //bgZJ2.push_back(HnuPlots::FileStruct(  "zjets" + tmode + "," + sample, mc_ZJ,         "hNu" + amode + uncert + "/" + cutlevels[5] + "/" + plot, lumi2012, xsecZJ,    k_mm_Zscale / NZJ, "hNu" + amode + uncert + "/" + cutlevels[4] + "/mLLZoom", 0.0, 0.0, true, 1, false, 60.0, 120.0));
    bgOther.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_tW,     "hNu" + amode + "/" +          cutlevels[cutlevel] + "/" + plot, lumi2012, xsectW,    1.0 / NtW,    "hNu" + amode + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_tbarW,  "hNu" + amode + "/" +          cutlevels[cutlevel] + "/" + plot, lumi2012, xsectbarW, 1.0 / NtbarW, "hNu" + amode + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_ZZ,     "hNu" + amode + "/" +          cutlevels[cutlevel] + "/" + plot, lumi2012, xsecZZ,    1.0 / NZZ,    "hNu" + amode + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_WZ,     "hNu" + amode + "/" +          cutlevels[cutlevel] + "/" + plot, lumi2012, xsecWZ,    1.0 / NWZ,    "hNu" + amode + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_WW,     "hNu" + amode + "/" +          cutlevels[cutlevel] + "/" + plot, lumi2012, xsecWW,    1.0 / NWW,    "hNu" + amode + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther2.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_tW,    "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, lumi2012, xsectW,    1.0 / NtW,    "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther2.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_tbarW, "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, lumi2012, xsectbarW, 1.0 / NtbarW, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther2.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_ZZ,    "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, lumi2012, xsecZZ,    1.0 / NZZ,    "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther2.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_WZ,    "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, lumi2012, xsecWZ,    1.0 / NWZ,    "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    bgOther2.push_back(HnuPlots::FileStruct("other" + tmode + "," + sample,  mc_WW,    "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, lumi2012, xsecWW,    1.0 / NWW,    "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    //signal  the xsecs and k-factors here are wrong because they do not matter for systematics!
    vsig10.push_back(HnuPlots::FileStruct("signal_700_350" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-700_MNu-350_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig11.push_back(HnuPlots::FileStruct("signal_800_400" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-800_MNu-400_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig12.push_back(HnuPlots::FileStruct("signal_900_450" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-900_MNu-450_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig13.push_back(HnuPlots::FileStruct("signal_1000_500" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig14.push_back(HnuPlots::FileStruct("signal_1100_550" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig15.push_back(HnuPlots::FileStruct("signal_1200_600" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1200_MNu-600_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig16.push_back(HnuPlots::FileStruct("signal_1300_650" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1300_MNu-650_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig17.push_back(HnuPlots::FileStruct("signal_1400_700" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1400_MNu-700_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig18.push_back(HnuPlots::FileStruct("signal_1500_750" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig19.push_back(HnuPlots::FileStruct("signal_1600_800" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1600_MNu-800_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig20.push_back(HnuPlots::FileStruct("signal_1700_850" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1700_MNu-850_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig21.push_back(HnuPlots::FileStruct("signal_1800_900" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1800_MNu-900_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig22.push_back(HnuPlots::FileStruct("signal_1900_950" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1900_MNu-950_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig23.push_back(HnuPlots::FileStruct("signal_2000_1000" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig24.push_back(HnuPlots::FileStruct("signal_2100_1050" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2100_MNu-1050_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig25.push_back(HnuPlots::FileStruct("signal_2200_1100" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2200_MNu-1100_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig26.push_back(HnuPlots::FileStruct("signal_2300_1150" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2300_MNu-1150_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig27.push_back(HnuPlots::FileStruct("signal_2400_1200" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2400_MNu-1200_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig28.push_back(HnuPlots::FileStruct("signal_2500_1250" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig29.push_back(HnuPlots::FileStruct("signal_2600_1300" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2600_MNu-1300_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig30.push_back(HnuPlots::FileStruct("signal_2700_1350" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2700_MNu-1350_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig31.push_back(HnuPlots::FileStruct("signal_2800_1400" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2800_MNu-1400_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig32.push_back(HnuPlots::FileStruct("signal_2900_1450" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig33.push_back(HnuPlots::FileStruct("signal_3000_1500" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig34.push_back(HnuPlots::FileStruct("signal_3100_1550" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3100_MNu-1550_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig35.push_back(HnuPlots::FileStruct("signal_3200_1600" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3200_MNu-1600_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig36.push_back(HnuPlots::FileStruct("signal_3300_1650" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3300_MNu-1650_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig37.push_back(HnuPlots::FileStruct("signal_3400_1700" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3400_MNu-1700_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig38.push_back(HnuPlots::FileStruct("signal_3500_1750" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3500_MNu-1750_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig39.push_back(HnuPlots::FileStruct("signal_3600_1800" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3600_MNu-1800_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig40.push_back(HnuPlots::FileStruct("signal_3700_1850" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3700_MNu-1850_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig41.push_back(HnuPlots::FileStruct("signal_3800_1900" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3800_MNu-1900_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig42.push_back(HnuPlots::FileStruct("signal_3900_1950" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3900_MNu-1950_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig43.push_back(HnuPlots::FileStruct("signal_4000_2000" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-4000_MNu-2000_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode +           "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode +          "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));

    //signal  the xsecs and k-factors here are wrong because they do not matter for systematics!
    vsig102.push_back(HnuPlots::FileStruct("signal_700_350" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-700_MNu-350_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig112.push_back(HnuPlots::FileStruct("signal_800_400" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-800_MNu-400_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig122.push_back(HnuPlots::FileStruct("signal_900_450" + tmode + "," + sample,   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-900_MNu-450_TuneZ2star_8TeV-pythia6-tauola.root",   "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig132.push_back(HnuPlots::FileStruct("signal_1000_500" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig142.push_back(HnuPlots::FileStruct("signal_1100_550" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig152.push_back(HnuPlots::FileStruct("signal_1200_600" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1200_MNu-600_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig162.push_back(HnuPlots::FileStruct("signal_1300_650" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1300_MNu-650_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig172.push_back(HnuPlots::FileStruct("signal_1400_700" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1400_MNu-700_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig182.push_back(HnuPlots::FileStruct("signal_1500_750" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig192.push_back(HnuPlots::FileStruct("signal_1600_800" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1600_MNu-800_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig202.push_back(HnuPlots::FileStruct("signal_1700_850" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1700_MNu-850_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig212.push_back(HnuPlots::FileStruct("signal_1800_900" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1800_MNu-900_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig222.push_back(HnuPlots::FileStruct("signal_1900_950" + tmode + "," + sample,  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1900_MNu-950_TuneZ2star_8TeV-pythia6-tauola.root",  "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig232.push_back(HnuPlots::FileStruct("signal_2000_1000" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig242.push_back(HnuPlots::FileStruct("signal_2100_1050" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2100_MNu-1050_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig252.push_back(HnuPlots::FileStruct("signal_2200_1100" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2200_MNu-1100_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig262.push_back(HnuPlots::FileStruct("signal_2300_1150" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2300_MNu-1150_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig272.push_back(HnuPlots::FileStruct("signal_2400_1200" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2400_MNu-1200_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig282.push_back(HnuPlots::FileStruct("signal_2500_1250" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig292.push_back(HnuPlots::FileStruct("signal_2600_1300" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2600_MNu-1300_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig302.push_back(HnuPlots::FileStruct("signal_2700_1350" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2700_MNu-1350_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig312.push_back(HnuPlots::FileStruct("signal_2800_1400" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2800_MNu-1400_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig322.push_back(HnuPlots::FileStruct("signal_2900_1450" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig332.push_back(HnuPlots::FileStruct("signal_3000_1500" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig342.push_back(HnuPlots::FileStruct("signal_3100_1550" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3100_MNu-1550_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig352.push_back(HnuPlots::FileStruct("signal_3200_1600" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3200_MNu-1600_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig362.push_back(HnuPlots::FileStruct("signal_3300_1650" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3300_MNu-1650_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig372.push_back(HnuPlots::FileStruct("signal_3400_1700" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3400_MNu-1700_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig382.push_back(HnuPlots::FileStruct("signal_3500_1750" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3500_MNu-1750_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig392.push_back(HnuPlots::FileStruct("signal_3600_1800" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3600_MNu-1800_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig402.push_back(HnuPlots::FileStruct("signal_3700_1850" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3700_MNu-1850_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig412.push_back(HnuPlots::FileStruct("signal_3800_1900" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3800_MNu-1900_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig422.push_back(HnuPlots::FileStruct("signal_3900_1950" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3900_MNu-1950_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));
    vsig432.push_back(HnuPlots::FileStruct("signal_4000_2000" + tmode + "," + sample, "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_rerecoData/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-4000_MNu-2000_TuneZ2star_8TeV-pythia6-tauola.root", "hNu" + amode + uncert + "/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "hNu" + amode + uncert + "/cut0_none/n_vertex", 0.0, 0.0, true, 1, false, 0.0, 1000.0));

    bg.push_back(bgZJ);
    bg.push_back(bgZJ2);
    bg.push_back(bgOther);
    bg.push_back(bgOther2);
    if(mode == 0 || mode == 1)
    {

        bg.push_back(vsig10);
        bg.push_back(vsig102);
        bg.push_back(vsig11);
        bg.push_back(vsig112);
        bg.push_back(vsig12);
        bg.push_back(vsig122);
        bg.push_back(vsig13);
        bg.push_back(vsig132);
        bg.push_back(vsig14);
        bg.push_back(vsig142);
        bg.push_back(vsig15);
        bg.push_back(vsig152);
        bg.push_back(vsig16);
        bg.push_back(vsig162);
        bg.push_back(vsig17);
        bg.push_back(vsig172);
        bg.push_back(vsig18);
        bg.push_back(vsig182);
        bg.push_back(vsig19);
        bg.push_back(vsig192);
        bg.push_back(vsig20);
        bg.push_back(vsig202);
        bg.push_back(vsig21);
        bg.push_back(vsig212);
        bg.push_back(vsig22);
        bg.push_back(vsig222);
        bg.push_back(vsig23);
        bg.push_back(vsig232);
        bg.push_back(vsig24);
        bg.push_back(vsig242);
        bg.push_back(vsig25);
        bg.push_back(vsig252);
        bg.push_back(vsig26);
        bg.push_back(vsig262);
        bg.push_back(vsig27);
        bg.push_back(vsig272);
        bg.push_back(vsig28);
        bg.push_back(vsig282);
        bg.push_back(vsig29);
        bg.push_back(vsig292);
        bg.push_back(vsig30);
        bg.push_back(vsig302);
        bg.push_back(vsig31);
        bg.push_back(vsig312);
        bg.push_back(vsig32);
        bg.push_back(vsig322);
        bg.push_back(vsig33);
        bg.push_back(vsig332);
        bg.push_back(vsig34);
        bg.push_back(vsig342);
        bg.push_back(vsig35);
        bg.push_back(vsig352);
        bg.push_back(vsig36);
        bg.push_back(vsig362);
        bg.push_back(vsig37);
        bg.push_back(vsig372);
        bg.push_back(vsig38);
        bg.push_back(vsig382);
        bg.push_back(vsig39);
        bg.push_back(vsig392);
        bg.push_back(vsig40);
        bg.push_back(vsig402);
        bg.push_back(vsig41);
        bg.push_back(vsig412);
        bg.push_back(vsig42);
        bg.push_back(vsig422);
        bg.push_back(vsig43);
        bg.push_back(vsig432);
    }
}

void plotMCSystCalc(int mode = 0)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile

    vector<pair<string, int> > uncerts;

    switch(mode)
    {
        case 0:
            uncerts.push_back(make_pair("jesHi", 1));
            uncerts.push_back(make_pair("jesLo", 1));
            uncerts.push_back(make_pair("jerHi", 1));
            uncerts.push_back(make_pair("jerLo", 1));
            //uncerts.push_back("mer");
            uncerts.push_back(make_pair("midHi", 1));
            uncerts.push_back(make_pair("midLo", 1));
            uncerts.push_back(make_pair("trigHi", 1));
            uncerts.push_back(make_pair("trigLo", 1));
            uncerts.push_back(make_pair("puHi", 1));
            uncerts.push_back(make_pair("puLo", 1));
            break;
        case 1:
            uncerts.push_back(make_pair("jesHi", 1));
            uncerts.push_back(make_pair("jesLo", 1));
            uncerts.push_back(make_pair("jerHi", 1));
            uncerts.push_back(make_pair("jerLo", 1));
            uncerts.push_back(make_pair("escale", 2));
            uncerts.push_back(make_pair("idHi", 1));
            uncerts.push_back(make_pair("idLo", 1));
            uncerts.push_back(make_pair("er", 2));
            uncerts.push_back(make_pair("puHi", 1));
            uncerts.push_back(make_pair("puLo", 1));
            break;
    }

    map<string, vector<double> > systMap;
    for(vector<pair<string, int> >::const_iterator i = uncerts.begin(); i != uncerts.end(); ++i)
    {
        vector<vector<HnuPlots::FileStruct> > bg, sig;
        mcSystCalcSetBg(mode, bg, i->first, i->first, 5, "mWR");
        HnuPlots hps(data, bg, sig, lumi);
        hps.mcSystCalc(i->second, &systMap);
    }

    //combime the highs and the lows
    for(map<string, vector<double> >::iterator highit = systMap.begin(); highit != systMap.end(); )
    {
        size_t reploc = highit->first.find("Hi", highit->first.size() - 2);
        if(reploc != string::npos)
        {
            string key = highit->first;
            key.erase(reploc, string::npos);
            map<string, vector<double> >::iterator lowit = systMap.find(key + "Lo");
            for(unsigned int i = 0; i < highit->second.size() && i < lowit->second.size(); i++)
            {
                double high = fabs(highit->second[i] - 1);
                double low = fabs(lowit->second[i] - 1);
                systMap[key].push_back(1 + (high + low) / 2);
            }
            systMap.erase(lowit);
            systMap.erase(highit);
            highit = systMap.begin();
        }
        else ++highit;
    }

    for(map<string, vector<double> >::const_iterator mit = systMap.begin(); mit != systMap.end(); ++mit)
    {
        printf("%s", mit->first.c_str());
        for(vector<double>::const_iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
        {

            printf(",%f", *vit);
        }
        printf("\n");
    }
}

void plotall()
{
    plot2012(0, 4, "mWR;");
    plot2012(0, 5, "mWR;");
    plot2012(0, 6, "mWR;");
    plot2012(0, 4, "mLL;");
    plot2012(0, 5, "mLL;");
    plot2012(0, 6, "mLL;");
    plot2012(0, 4, "mJJ;");
    plot2012(0, 5, "mJJ;");
    plot2012(0, 6, "mJJ;");
    plot2012(0, 4, "mNuR1;");
    plot2012(0, 5, "mNuR1;");
    plot2012(0, 6, "mNuR1;");
    plot2012(0, 4, "mNuR2;");
    plot2012(0, 5, "mNuR2;");
    plot2012(0, 6, "mNuR2;");
    plot2012(0, 4, "ptL1;");
    plot2012(0, 5, "ptL1;");
    plot2012(0, 6, "ptL1;");
    plot2012(0, 4, "ptL2;");
    plot2012(0, 5, "ptL2;");
    plot2012(0, 6, "ptL2;");
    plot2012(0, 4, "ptJ1;");
    plot2012(0, 5, "ptJ1;");
    plot2012(0, 6, "ptJ1;");
    plot2012(0, 4, "ptJ2;");
    plot2012(0, 5, "ptJ2;");
    plot2012(0, 6, "ptJ2;");
    plot2012(7);
    plot2012(4);
    plot2012(0, 4, "mWR;mLL>120;mLL<200");
    plotMCFits(0);

    plot2012(1, 4, "mWR;");
    plot2012(1, 5, "mWR;");
    plot2012(1, 6, "mWR;");
    plot2012(1, 4, "mLL;");
    plot2012(1, 5, "mLL;");
    plot2012(1, 6, "mLL;");
    plot2012(1, 4, "mJJ;");
    plot2012(1, 5, "mJJ;");
    plot2012(1, 6, "mJJ;");
    plot2012(1, 4, "mNuR1;");
    plot2012(1, 5, "mNuR1;");
    plot2012(1, 6, "mNuR1;");
    plot2012(1, 4, "mNuR2;");
    plot2012(1, 5, "mNuR2;");
    plot2012(1, 6, "mNuR2;");
    plot2012(1, 4, "ptL1;");
    plot2012(1, 5, "ptL1;");
    plot2012(1, 6, "ptL1;");
    plot2012(1, 4, "ptL2;");
    plot2012(1, 5, "ptL2;");
    plot2012(1, 6, "ptL2;");
    plot2012(1, 4, "ptJ1;");
    plot2012(1, 5, "ptJ1;");
    plot2012(1, 6, "ptJ1;");
    plot2012(1, 4, "ptJ2;");
    plot2012(1, 5, "ptJ2;");
    plot2012(1, 6, "ptJ2;");
    plot2012(1, 4, "etaL1;", 5, true, -3, 3);
    plot2012(1, 5, "etaL1;", 5, true, -3, 3);
    plot2012(1, 6, "etaL1;", 5, true, -3, 3);
    plot2012(1, 4, "etaL2;", 5, true, -3, 3);
    plot2012(1, 5, "etaL2;", 5, true, -3, 3);
    plot2012(1, 6, "etaL2;", 5, true, -3, 3);
    plot2012(1, 4, "etaJ1;", 5, true, -3, 3);
    plot2012(1, 5, "etaJ1;", 5, true, -3, 3);
    plot2012(1, 6, "etaJ1;", 5, true, -3, 3);
    plot2012(1, 4, "etaJ2;", 5, true, -3, 3);
    plot2012(1, 5, "etaJ2;", 5, true, -3, 3);
    plot2012(1, 6, "etaJ2;", 5, true, -3, 3);
    plot2012(1, 4, "phiL1;", 2, true, -4, 4);
    plot2012(1, 5, "phiL1;", 2, true, -4, 4);
    plot2012(1, 6, "phiL1;", 2, true, -4, 4);
    plot2012(1, 4, "phiL2;", 2, true, -4, 4);
    plot2012(1, 5, "phiL2;", 2, true, -4, 4);
    plot2012(1, 6, "phiL2;", 2, true, -4, 4);
    plot2012(1, 4, "phiJ1;", 2, true, -4, 4);
    plot2012(1, 5, "phiJ1;", 2, true, -4, 4);
    plot2012(1, 6, "phiJ1;", 2, true, -4, 4);
    plot2012(1, 4, "phiJ2;", 2, true, -4, 4);
    plot2012(1, 5, "phiJ2;", 2, true, -4, 4);
    plot2012(1, 6, "phiJ2;", 2, true, -4, 4);
    plot2012(8);
    plot2012(5);
    plot2012(1, 4, "mWR;mLL>120;mLL<200");
    plot2012(1, 4, "mLL;mWR>900;mWR<1200");
    plot2012(1, 4, "mLL;!mWR>900;mWR<1200");
    plot2012(1, 5, "dPhiL;mWR>900;mWR<1200", 3, true, -4, 4);
    plot2012(1, 5, "dEtaL;mWR>900;mWR<1200", 3, true, -4, 4);
    plot2012(1, 5, "dPhiJ;mWR>900;mWR<1200", 3, true, -4, 4);
    plot2012(1, 5, "dEtaJ;mWR>900;mWR<1200", 3, true, -4, 4);
    plot2012(1, 4, "pL1;");
    plot2012(1, 5, "mWR;pL1>600;pL1<900");
    plot2012(1, 5, "mWR;!pL1>600;pL1<900");
    plot2012(1, 4, "mLL;pL1>600;pL1<900");
    plot2012(1, 4, "mLL;!pL1>600;pL1<900");
    plot2012(1, 5, "n_vertex;");
    plot2012(1, 5, "n_vertex;mWR>900;mWR<1200");
    plot2012(1, 5, "n_vertex;mWR>1800;mWR<2200");
    plot2012(1, 5, "mWR;n_vertex<15");
    plot2012(1, 5, "mWR;!n_vertex<15");
    plotMCFits(1);
    plot2012(1, 4, "mLL;mWR>1800;mWR<2200");
    plot2012(1, 5, "mNuR1;mWR>1800;mWR<2200");
    plot2012(1, 5, "mNuR2;mWR>1800;mWR<2200");
    plot2012(1, 5, "ptL1;mWR>1800;mWR<2200");
    plot2012(1, 5, "ptL2;mWR>1800;mWR<2200");
    plot2012(1, 5, "ptJ1;mWR>1800;mWR<2200");
    plot2012(1, 5, "ptJ2;mWR>1800;mWR<2200");
    plot2012(1, 5, "etaL1;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "etaL2;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "etaJ1;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "etaJ2;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "phiL1;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "phiL2;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "phiJ1;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "phiJ2;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "mJJ;mWR>1800;mWR<2200");
    plot2012(1, 5, "mOuR1;mWR>1800;mWR<2200");
    plot2012(1, 5, "mOuR2;mWR>1800;mWR<2200");
    plot2012(1, 5, "pL1;mWR>1800;mWR<2200");
    plot2012(1, 5, "pL2;mWR>1800;mWR<2200");
    plot2012(1, 5, "dPhiL;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "dEtaL;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "dPhiJ;mWR>1800;mWR<2200", 3, true, -4, 4);
    plot2012(1, 5, "dEtaJ;mWR>1800;mWR<2200", 3, true, -4, 4);


    plot2012(2, 5, "mWR");
    plot2012(3, 5, "mWR;");

    printf("\n\nZJ norm\nMuon\n");
    plotZJNorm(true);
    printf("\n");
    printf("Electron\n");
    plotZJNorm(false);
    printf("\n");

    printf("\nTT 1b\nMuon\n");
    plotDDTTNorm(true, 1);
    printf("\n");
    printf("Electron\n");
    plotDDTTNorm(false, 1);
    printf("\n");

    printf("\nTT 2b\nMuon\n");
    plotDDTTNorm(true, 2);
    printf("\n");
    printf("Electron\n");
    plotDDTTNorm(false, 2);
    printf("\n");

    printf("\nTT eff\n");
    plotTTBarDDEffBasedNorm();
    printf("\n");

    printf("\nTT MC\nMuon\n");
    plotTTBarMCNorm(true, 5, 0, true);
    printf("\n");
    printf("Electron\n");

    plotTTBarMCNorm(false, 5, 0, true);
    printf("\n");

    printf("\nMuon Cut Flow\n");
    plotCutFlow(0);
    printf("\n");
    printf("\nElectron Cut Flow\n");
    plotCutFlow(1);
    printf("\n");
}

int main()
{
    plot2012(0, 5);
    //plotall();
}
