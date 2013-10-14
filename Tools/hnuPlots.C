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
    kCyan + 2,
    kOrange + 1,
    kGreen - 3,
    kYellow + 4,
    kMagenta - 1,
    kRed,
    kBlue,
    kGreen
};
const int NCOLORS = sizeof(colors) / sizeof(int);

const int shcolors[] = {
    kRed,
    kBlue,
    kGreen,
    kYellow,
    kMagenta,
    kOrange
};
const int NSHCOLORS = sizeof(shcolors) / sizeof(int);

const int hatchs[] = {
    3004,
    3005,
    3002,
    3013,
    3006,
    3007,
    3017,
    3018
};
const int NHATCHS = sizeof(hatchs) / sizeof(int);

double bins[] = {600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2200.0, 4000.0};
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

        FileStruct();
        FileStruct(std::string l, std::string f, std::string h);
        FileStruct(std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh = "", double cl = 0.0, double ch = 0.0, bool px = true, int nb = 1, bool un = true, double nll = 0.0, double nul = 0.0, bool loadtuple = false, bool lhft = false, double ll = 0.0, double ul = 1.0, double bn = 1.0, HeavyNuTree::HNuSlopeFitInfo* tll = 0, HeavyNuTree::HNuSlopeFitInfo* tul = 0);
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

        HistStruct();
        HistStruct(std::string l, TH1* h, TH1* nh = NULL, double nll = 0.0, double nul = 0.0);
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
    void setSavePlots(bool sp);
    void autoSetHistogramAxisTitle(int mode = 0);

private:
    std::vector<HistStruct> bghists;
    std::vector<HistStruct> sighists;
    HistStruct datahist;

    std::vector<std::vector<float> > systematics;
    std::vector<float> shapeerr;

    int rebin, nhist;
    double iLumi, xmin, xmax;
    std::string xaxislabel;
    std::string yaxislabel;
    std::string formlabel;
    bool autosort, islog, saveplots, plotSMoData;
    double sigscale;

    TH1* project(TH2* h2d, double cl, double ch, bool porjx = true);
    int projcount;

    TH1* histFromTuple(std::string label, std::string histpath, double thll, double thul, double nb, std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >& bgtvec, HeavyNuTree::HNuSlopeFitInfo *ll = NULL, HeavyNuTree::HNuSlopeFitInfo *ul = NULL);
    bool runFilter(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iE);
    void histFromDataCard(std::map<std::pair<std::string, std::string>, std::vector<float> >& uncerts);

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
                    bgtvec.push_back(std::make_pair(hnt->event_, scale));
                    ibgtvec.push_back(std::make_pair(hnt->event_, scale));
                }
                while(hnt->GetNextEvent());
                delete [] hnt;
            }

            // gethistogram
            if(ibgf->loadtuple && ibgf->histFromTuple)
            {
                h = histFromTuple(ibgf->label, ibgf->histpath.substr(ibgf->histpath.rfind("/") + 1, ibgf->histpath.size()), ibgf->thll, ibgf->thul, ibgf->thb, ibgtvec, ibgf->tpll, ibgf->tpul);
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

            TH1 * h;
            if(isigf->loadtuple && isigf->histFromTuple)
            {
                h = histFromTuple(isigf->label, isigf->histpath.substr(isigf->histpath.rfind("/") + 1, isigf->histpath.size()), isigf->thll, isigf->thul, isigf->thb, sigtvec, isigf->tpll, isigf->tpul);
            }
            else h = (TH1*)file->Get(isigf->histpath.c_str());
            //if(h) h = (TH1*)h->Clone();
            if(h && hn)
            {
                std::cout << hn->GetBinContent(isigf->normbin) << std::endl;
                h->Scale(isigf->intLumi * isigf->cs * isigf->kfactor / hn->GetBinContent(isigf->normbin));
                if(first)
                {
                    sighists.push_back(HnuPlots::HistStruct(isigf->label, (TH1*)h->Clone()));
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
            h = histFromTuple(fdata.label, fdata.histpath.substr(fdata.histpath.rfind("/") + 1, fdata.histpath.size()), fdata.thll, fdata.thul, fdata.thb, dtvec, fdata.tpll, fdata.tpul);
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

HnuPlots::HistStruct::HistStruct(std::string l, TH1* h, TH1* nh, double nll, double nul)
{

    label = l;
    hist = h;
    normhist = nh;
    normll = nll;
    normul = nul;
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

TH1* HnuPlots::histFromTuple(std::string label, std::string histValues, double thll, double thul, double nb, std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >& bgtvec, HeavyNuTree::HNuSlopeFitInfo *ll, HeavyNuTree::HNuSlopeFitInfo *ul)
{    
    std::vector<std::string> histQs;
    for(size_t pos = 0;;) 
    {
        size_t npos = histValues.find(':', pos + 1);
        histQs.push_back(histValues.substr(pos, npos - pos));
        if(npos == size_t(-1)) break;
        pos = npos + 1;
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
            if(histV.compare("mWR") == 0) vlim.push_back(Limits(0.0, 4000.0, 100));
            else if(histV.compare("mLL") == 0 || histV.compare("mJJ") == 0) vlim.push_back(Limits(0.0, 2000.0, 100));
            else if(histV.compare("mNuR1") == 0 || histV.compare("mNuR2") == 0) vlim.push_back(Limits(0.0, 3000.0, 150));
            else if(histV.compare("mOuR1") == 0 || histV.compare("mOuR2") == 0) vlim.push_back(Limits(0.0, 2000.0, 100));
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
        }
    }

    TH1 *hist = 0;
    
    if(vlim.size() == 1)       hist = new TH1D(hname, hname, vlim[0].nb, vlim[0].thll, vlim[0].thul);
    else if (vlim.size() == 2) hist = new TH2D(hname, htitle, vlim[0].nb, vlim[0].thll, vlim[0].thul, vlim[1].nb, vlim[1].thll, vlim[1].thul);
    else 
    {
        printf("!!!Too many histogram dimmensions!!!\n");
        return 0;
    }

    for(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iT = bgtvec.begin(); iT != bgtvec.end(); ++iT)
    {
        std::vector<double> values;
        
        if(ll)
        {
            if((iT->first.cutlevel < ll->cutlevel)) continue;
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
            if((iT->first.cutlevel > ul->cutlevel)) continue;
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
        if(iT->first.weight > 1000 || iT->first.weight < 0.001) continue;

        double pL1 = iT->first.l1pt * cosh(iT->first.l1eta);
        double pL2 = iT->first.l2pt * cosh(iT->first.l2eta);
        double pLmin = std::min(pL1, pL2);
        double pLmax = std::max(pL1, pL2);

        double gR9L1 = iT->first.sE1 - iT->first.rhE1;
        double gR9L2 = iT->first.sE2 - iT->first.rhE2;
        double rhoScL1 = iT->first.rhE1 / iT->first.sE1;
        double rhoScL2 = iT->first.rhE2 / iT->first.sE2;

        //TLorentzVector L1, L2;
        //L1.SetPtEtaPhiM(iT->first.l1pt, iT->first.l1eta, iT->first.l1phi, 0);
        //L2.SetPtEtaPhiM(iT->first.l2pt, iT->first.l2eta, iT->first.l2phi, 0);
        //if(fabs(ROOT::Math::VectorUtil::DeltaPhi(L1, L2)) > 3) continue;

        //if(std::max(iT->first.sE1, iT->first.sE2) < 600 || std::max(iT->first.sE1, iT->first.sE2) > 900) continue;
        //if(pLmax < 500 || pLmax > 900) continue;
        //if(iT->first.mlljj < 900 || iT->first.mlljj > 1200) continue;
        if(iT->first.mlljj < 1800 || iT->first.mlljj > 2200) continue;
        //if(iT->first.mlljj > 800 && iT->first.mlljj < 1200) continue;
        //if(iT->first.mll < 150) continue;
        //if(pL2 > 700) continue;
        //if(iT->first.mlljj < 1000 || iT->first.mlljj > 1200) continue;
        //if(iT->first.mll < 500) continue;
        //if(iT->first.rhE1 / iT->first.sE1 < 0.1) continue;
        //if(fabs(iT->first.l1eta) < 1.5 && fabs(iT->first.l2eta) < 1.5) continue; //NOT EB-EB
        //if(fabs(iT->first.l1eta) > 1.5 && fabs(iT->first.l2eta) > 1.5) continue; //NOT EE-EE
        //if(fabs(iT->first.l1eta) > 1.5 || fabs(iT->first.l2eta) > 1.5) continue; //EB-EB
        //if(fabs(iT->first.l1eta) < 1.5 || fabs(iT->first.l2eta) < 1.5) continue; //EE-EE
        //if(iT->first.n_primaryVertex < 15) continue;
        //if(!runFilter(iT)) continue;
        

        int nbjet = 0;
        if(iT->first.j1B + iT->first.j2B < nbjet) continue;

        for(std::vector<std::string>::const_iterator ihlabel = histQs.begin(); ihlabel != histQs.end(); ++ihlabel)
        {
            std::string histpath = *ihlabel;
            if(histpath.compare("mWR") == 0)        values.push_back(iT->first.mlljj);
            else if(histpath.compare("mLL") == 0)   values.push_back(iT->first.mll  );
            else if(histpath.compare("mNuR1") == 0)
            {
                TLorentzVector J1, J2, L;
                J1.SetPtEtaPhiM(iT->first.j1pt, iT->first.j1eta, iT->first.j1phi, 0);
                J2.SetPtEtaPhiM(iT->first.j2pt, iT->first.j2eta, iT->first.j2phi, 0);
                L.SetPtEtaPhiM(iT->first.l1pt, iT->first.l1eta, iT->first.l1phi, 0);
                values.push_back((J1 + J2 + L).M());
            }
            else if(histpath.compare("mNuR2") == 0)
            {
                TLorentzVector J1, J2, L;
                J1.SetPtEtaPhiM(iT->first.j1pt, iT->first.j1eta, iT->first.j1phi, 0);
                J2.SetPtEtaPhiM(iT->first.j2pt, iT->first.j2eta, iT->first.j2phi, 0);
                L.SetPtEtaPhiM(iT->first.l2pt, iT->first.l2eta, iT->first.l2phi, 0);
                values.push_back((J1 + J2 + L).M());
            }
            else if(histpath.compare("mOuR1") == 0)
            {
                TLorentzVector L1, L2, J;
                L1.SetPtEtaPhiM(iT->first.l1pt, iT->first.l1eta, iT->first.l1phi, 0);
                L2.SetPtEtaPhiM(iT->first.l2pt, iT->first.l2eta, iT->first.l2phi, 0);
                J.SetPtEtaPhiM(iT->first.j1pt, iT->first.j1eta, iT->first.j1phi, 0);
                values.push_back((L1 + L2 + J).M());
            }
            else if(histpath.compare("mOuR2") == 0)
            {
                TLorentzVector L1, L2, J;
                L1.SetPtEtaPhiM(iT->first.l1pt, iT->first.l1eta, iT->first.l1phi, 0);
                L2.SetPtEtaPhiM(iT->first.l2pt, iT->first.l2eta, iT->first.l2phi, 0);
                J.SetPtEtaPhiM(iT->first.j2pt, iT->first.j2eta, iT->first.j2phi, 0);
                values.push_back((L1 + L2 + J).M());
            }
            else if(histpath.compare("ptL1") == 0)  values.push_back(iT->first.l1pt );
            //else if(histpath.compare("etaL1") == 0) values.push_back((pL1 > pL2)?iT->first.l1eta:iT->first.l2eta);
            //else if(histpath.compare("phiL1") == 0) values.push_back((pL1 > pL2)?iT->first.l1phi:iT->first.l2phi);
            else if(histpath.compare("etaL1") == 0) values.push_back(iT->first.l1eta);
            else if(histpath.compare("phiL1") == 0) values.push_back(iT->first.l1phi);
            else if(histpath.compare("ptL2") == 0)  values.push_back(iT->first.l2pt );
            //else if(histpath.compare("etaL2") == 0) values.push_back((pL1 < pL2)?iT->first.l1eta:iT->first.l2eta);
            //else if(histpath.compare("phiL2") == 0) values.push_back((pL1 < pL2)?iT->first.l1phi:iT->first.l2phi);
            else if(histpath.compare("etaL2") == 0) values.push_back(iT->first.l2eta);
            else if(histpath.compare("phiL2") == 0) values.push_back(iT->first.l2phi);
            else if(histpath.compare("ptJ1") == 0)  values.push_back(iT->first.j1pt );
            else if(histpath.compare("etaJ1") == 0) values.push_back(iT->first.j1eta);
            else if(histpath.compare("phiJ1") == 0) values.push_back(iT->first.j1phi);
            else if(histpath.compare("ptJ2") == 0)  values.push_back(iT->first.j2pt );
            else if(histpath.compare("etaJ2") == 0) values.push_back(iT->first.j2eta);
            else if(histpath.compare("phiJ2") == 0) values.push_back(iT->first.j2phi);
            else if(histpath.compare("pL1") == 0)   values.push_back(pLmax);
            else if(histpath.compare("pL2") == 0)   values.push_back(pLmin);
            else if(histpath.compare("mJJ") == 0)
            {
                TLorentzVector J1, J2;
                J1.SetPtEtaPhiM(iT->first.j1pt, iT->first.j1eta, iT->first.j1phi, 0);
                J2.SetPtEtaPhiM(iT->first.j2pt, iT->first.j2eta, iT->first.j2phi, 0);
                values.push_back((J1 + J2).M() );
            }
            else if(histpath.compare("n_vertex") == 0) values.push_back(iT->first.n_primaryVertex);
            else if(histpath.compare("seL1") == 0)     values.push_back(std::max(iT->first.sE1, iT->first.sE2));
            else if(histpath.compare("seL2") == 0)     values.push_back(std::min(iT->first.sE1, iT->first.sE2));
            else if(histpath.compare("rhL1") == 0)     values.push_back(std::max(iT->first.rhE1, iT->first.rhE2));
            else if(histpath.compare("rhL2") == 0)     values.push_back(std::min(iT->first.rhE1, iT->first.rhE2));
            else if(histpath.compare("gR9L1") == 0)    values.push_back((pL1 > pL2)?gR9L1:gR9L2);
            else if(histpath.compare("gR9L2") == 0)    values.push_back((pL1 < pL2)?gR9L1:gR9L2);
            else if(histpath.compare("rhoScL1") == 0)  values.push_back((pL1 > pL2)?rhoScL1:rhoScL2);
            else if(histpath.compare("rhoScL2") == 0)  values.push_back((pL1 < pL2)?rhoScL1:rhoScL2);
            else if(histpath.compare("met") == 0)      values.push_back(iT->first.met);
            else if(histpath.compare("dEtaL") == 0)    values.push_back(iT->first.l1eta - iT->first.l2eta);
            else if(histpath.compare("dEtaJ") == 0)    values.push_back(iT->first.j1eta - iT->first.j2eta);
            else if(histpath.compare("dPhiL") == 0)    
            {
                TLorentzVector L1, L2;
                L1.SetPtEtaPhiM(iT->first.l1pt, iT->first.l1eta, iT->first.l1phi, 0);
                L2.SetPtEtaPhiM(iT->first.l2pt, iT->first.l2eta, iT->first.l2phi, 0);
                values.push_back(ROOT::Math::VectorUtil::DeltaPhi(L1, L2));
            }
            else if(histpath.compare("dPhiJ") == 0)    
            {
                TLorentzVector J1, J2;
                J1.SetPtEtaPhiM(iT->first.j1pt, iT->first.j1eta, iT->first.j1phi, 0);
                J2.SetPtEtaPhiM(iT->first.j2pt, iT->first.j2eta, iT->first.j2phi, 0);
                values.push_back(ROOT::Math::VectorUtil::DeltaPhi(J1, J2));
            }
        }
        
        if(values.size() == 1) hist->Fill(values[0], iT->first.weight);
        else if(values.size() == 2) ((TH2*)hist)->Fill(values[0], values[1], iT->first.weight);
        
        //printf("%d:%d:%d\n", iT->first.run, iT->first.ls, iT->first.event);
        
    }

    return hist;
}

bool HnuPlots::runFilter(std::vector<std::pair<HeavyNuTree::HNuSlopeFitInfo, double> >::const_iterator iE)
{
    char pkey[128];
    sprintf(pkey,"%d:%d:%d", iE->first.run, iE->first.ls, iE->first.event);
    std::string key(pkey);
    
    std::set<std::string> events;
    
    events.insert("190705:125:138707430");
    events.insert("196203:111:91807338");
    events.insert("196453:10:9358199");
    events.insert("196453:883:761216455");
    events.insert("194912:323:545723590");
    events.insert("196199:520:569357512");
    events.insert("194825:154:166070019");
    events.insert("195013:402:601005884");
    events.insert("195390:24:19520685");
    events.insert("194076:734:719190420");
    events.insert("196453:1643:1438165030");
    events.insert("191264:157:164680827");
    events.insert("193541:287:174824785");
    events.insert("191247:107:167104820");
    events.insert("193621:221:197334826");
    events.insert("195552:1213:1413997178");
    events.insert("196197:523:777269889");
    events.insert("196250:39:73142375");
    events.insert("195530:390:554937459");
    events.insert("196438:598:492420742");
    events.insert("195163:201:272542209");
    events.insert("195013:359:536196304");
    events.insert("194076:261:280918454");
    events.insert("196452:630:875969154");
    events.insert("194464:142:208604846");
    events.insert("194455:55:32556189");
    events.insert("194210:58:58280648");
    events.insert("195950:709:663851870");
    events.insert("191277:53:51859398");
    events.insert("194429:160:155488974");
    events.insert("194224:153:245691275");
    events.insert("195776:304:249648204");
    events.insert("196199:51:62636583");
    events.insert("196364:890:812956406");
    events.insert("196452:190:226390737");
    events.insert("195390:249:197412270");
    events.insert("195251:144:254613585");
    events.insert("195396:105:118958995");
    events.insert("194480:270:234952983");
    events.insert("194428:235:246021234");
    events.insert("194050:1801:1561086221");
    events.insert("194643:223:283029706");
    events.insert("195930:303:262921310");
    events.insert("194643:73:96160284");
    events.insert("194533:101:116820658");
    events.insert("196334:167:233947563");
    events.insert("194224:63:105088730");
    events.insert("198212:251:140506263");
    events.insert("198208:122:95438269");
    events.insert("198230:376:395750377");
    events.insert("198269:88:145535146");
    events.insert("198271:609:707755328");
    events.insert("198271:574:668661571");
    events.insert("198487:383:448070154");
    events.insert("201191:153:194774190");
    events.insert("201191:1380:1738901057");
    events.insert("201191:855:1177128029");
    events.insert("201191:329:506815057");
    events.insert("202478:586:638731130");
    events.insert("202504:1342:1528576884");
    events.insert("202504:597:777603766");
    events.insert("202973:1066:943524588");
    events.insert("202973:32:35901862");
    events.insert("203002:551:728595239");
    events.insert("203002:1110:1345254735");
    events.insert("203002:682:875890936");
    events.insert("201625:766:993192102");
    events.insert("199008:76:46821012");
    events.insert("199336:396:371022277");
    events.insert("199336:113:112492563");
    events.insert("199356:197:194760721");
    events.insert("199429:3:2824720");
    events.insert("199435:559:638334633");
    events.insert("199435:312:361578631");
    events.insert("199435:637:725824161");
    events.insert("199436:557:371094848");
    events.insert("199569:345:395211885");
    events.insert("199569:46:57998469");
    events.insert("201707:707:862384763");
    events.insert("201707:387:508357418");
    events.insert("202013:14:22173400");
    events.insert("202045:430:486277357");
    events.insert("202060:557:698198797");
    events.insert("202060:335:425882011");
    events.insert("202060:238:287005652");
    events.insert("200075:443:531775950");
    events.insert("200091:214:242291765");
    events.insert("200091:586:709855001");
    events.insert("200244:59:96270806");
    events.insert("200245:93:108373628");
    events.insert("200525:609:784057891");
    events.insert("200525:786:974845738");
    events.insert("202087:738:881809665");
    events.insert("202178:901:1069311071");
    events.insert("202272:311:365859285");
    events.insert("202299:524:724043236");
    events.insert("202299:168:195058159");
    events.insert("202314:235:331270229");
    events.insert("200976:128:79776907");
    events.insert("200991:683:855410441");
    events.insert("200991:642:808094824");
    events.insert("200991:220:316237305");
    events.insert("200992:59:53661396");
    events.insert("201196:297:247709162");
    events.insert("201202:151:143333421");
    events.insert("201602:634:851797398");
    events.insert("201613:230:372265502");
    events.insert("199608:614:716998937");
    events.insert("199752:170:214118956");
    events.insert("199754:933:873005082");
    events.insert("199754:115:131824414");
    events.insert("199804:366:425220581");
    events.insert("199804:160:157408527");
    events.insert("199812:275:321510128");
    events.insert("199876:375:424041189");
    events.insert("199961:155:164977314");
    events.insert("199973:88:52408823");
    events.insert("200041:1022:1189939829");
    events.insert("206745:1144:1156506678");
    events.insert("206745:1745:1583890103");
    events.insert("206906:110:137197727");
    events.insert("207233:58:51316882");
    events.insert("207273:714:775754120");
    events.insert("207905:198:259518227");
    events.insert("207920:89:72140784");
    events.insert("208307:506:717647682");
    events.insert("208353:300:356001943");
    events.insert("208391:704:886577304");
    events.insert("208487:621:899301705");
    events.insert("208487:354:568995139");
    events.insert("208541:137:207923730");
    events.insert("206302:33:52193690");
    events.insert("206331:232:318273204");
    events.insert("206448:4:3259361");
    events.insert("206448:374:367483936");
    events.insert("206466:93:171933899");
    events.insert("206484:318:452791668");
    events.insert("206542:244:423881463");
    events.insert("206542:565:854555458");
    events.insert("206574:80:134760174");
    events.insert("206594:236:341472582");
    events.insert("208551:307:508743224");
    events.insert("208686:142:176573370");
    events.insert("204599:248:369145054");
    events.insert("204601:97:137400872");
    events.insert("205158:204:275312091");
    events.insert("205193:849:1087562044");
    events.insert("205193:881:1123797334");
    events.insert("205193:124:132784573");
    events.insert("205217:294:294904978");
    events.insert("205236:275:382572514");
    events.insert("205310:98:84737023");
    events.insert("205344:819:887286778");
    events.insert("205344:1378:1353831950");
    events.insert("205526:261:235416659");
    events.insert("205526:37:35531724");
    events.insert("205617:538:552856136");
    events.insert("205666:219:338960296");
    events.insert("203912:339:400823882");
    events.insert("203912:620:709658290");
    events.insert("203987:755:844756663");
    events.insert("204113:244:353151649");
    events.insert("204238:37:59447365");
    events.insert("204544:358:475186135");
    events.insert("204564:676:739538862");
    events.insert("204564:88:106661836");
    events.insert("204577:518:640366965");
    events.insert("207397:130:218826419");
    events.insert("207487:60:69012218");
    events.insert("207515:525:822748782");
    events.insert("207515:736:1103307498");
    events.insert("207515:144:182196526");
    events.insert("207515:560:872945761");
    events.insert("207884:45:69726318");
    events.insert("207884:109:168129245");
    events.insert("207905:963:1261676675");
    events.insert("205718:357:602911878");
    events.insert("205774:16:25437700");
    events.insert("205781:73:105816168");
    events.insert("206187:601:881513255");
    events.insert("206207:487:720837075");
    events.insert("206210:232:217763299");
    events.insert("206210:272:252752342");
    events.insert("206243:816:1129109305");
    events.insert("206246:504:485730004");


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

HnuPlots::FileStruct::FileStruct(std::string l, std::string f, std::string h, double iL, double c, double kf, std::string nh, double cl, double ch, bool px, int nb, bool un, double nll, double nul, bool lt, bool lhft, double ll, double ul, double bn, HeavyNuTree::HNuSlopeFitInfo* tll, HeavyNuTree::HNuSlopeFitInfo * tul)
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
    gStyle->SetHatchesSpacing(0.7);

    if(rebin > 1)
    {
        datahist.hist->Rebin(rebin);
        for(vector<HnuPlots::HistStruct>::const_iterator ihbg = bghists.begin(); ihbg != bghists.end(); ihbg++)
        {
            ihbg->hist->Rebin(rebin);
        }
        for(vector<HnuPlots::HistStruct>::const_iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
        {
            ihsig->hist->Rebin(rebin);
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

    if(!yaxislabel.compare("please auto set the axis"))
    {
        char temp[128];
        if(xaxislabel.find("GeV") < xaxislabel.size())
        {
            if(rebin >= 0) sprintf(temp, "Events / %.0f GeV", datahist.hist->GetBinWidth(1));
            else sprintf(temp, "Events / GeV");
            yaxislabel = temp;
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
        fontScale = 8.0 / 9;
    }
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.06);
    //c1->SetMargin(0.15, 0.1, 0.1, 0.1);

    //TLegend *leg = new TLegend(0.52, 0.67, 0.94, 0.91);
    TLegend *leg = new TLegend(0.52, 0.67, 0.90, 0.91);
    leg->SetFillStyle(0); //Color(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    float dataintegral = 0.0;
    if(rebin >= 0) dataintegral = datahist.hist->Integral(0, datahist.hist->GetNbinsX() + 1);
    else dataintegral = datahist.hist->Integral(1, datahist.hist->GetNbinsX(), "width");
    char datahllabel[128];
    sprintf(datahllabel, "%s (%.0f)", datahist.label.c_str(), dataintegral);
    leg->AddEntry(datahist.hist, datahllabel);

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
        else integral = ihbg->hist->Integral(1, ihbg->hist->GetNbinsX(), "width");
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)", ihbg->label.c_str(), floor(integral + 0.5));
        leg->AddEntry(ihbg->hist, hllabel);
    }
    //double sigMaxMin = datahist.hist->GetMaximum();
    for(vector<HnuPlots::HistStruct>::const_iterator ihsig = sighists.begin(); ihsig != sighists.end(); ihsig++)
    {
        float integral = 0.0;
        if(rebin >= 0) integral = ihsig->hist->Integral(0, ihsig->hist->GetNbinsX() + 1);
        else integral = ihsig->hist->Integral(1, ihsig->hist->GetNbinsX(), "width");
        char hllabel[128];
        sprintf(hllabel, "%s (%.0f)", ihsig->label.c_str(), floor(integral + 0.5));
        leg->AddEntry(ihsig->hist, hllabel);
        //sigMaxMin = std::min(ihsig->hist->GetMaximum(), sigMaxMin);
    }

    TH1 *dummy = new TH1F("dummy", "dummy", 1000, datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
    if(xmin != xmax) dummy->GetXaxis()->SetRangeUser(xmin, xmax);
    //dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    if(islog)
    {
        dummy->GetYaxis()->SetRangeUser(std::max(0.0001, 0.2 * std::min(hbg->GetMaximum(), 0.3 * datahist.hist->GetMinimum(0.0001))), std::max(hbg->GetMaximum(), datahist.hist->GetMaximum())*4);
        gPad->SetLogy(1);
    }
    else
    {
        dummy->GetYaxis()->SetRangeUser(0.001, std::max(hbg->GetMaximum(), datahist.hist->GetMaximum())*1.2);
    }
    dummy->GetYaxis()->SetTitle(yaxislabel.c_str());
    dummy->GetYaxis()->SetTitleOffset(1.05 / fontScale);
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->SetStats(0);
    if(plotSMoData) dummy->GetXaxis()->SetTitle(0);
    else            dummy->GetXaxis()->SetTitle(xaxislabel.c_str());
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    TLatex mark;
    mark.SetTextSize(0.04 * 7 / 6.5 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    datahist.hist->SetMarkerColor(kBlack);
    datahist.hist->SetMarkerStyle(20);
    datahist.hist->SetLineWidth(2.0);

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
        isig->hist->Draw("hist same");
    }
    datahist.hist->Draw("same pe");
    fixOverlay();
    leg->Draw("same");
    mark.DrawLatex(0.15, 0.95, "CMS Preliminary");
    mark.DrawLatex(0.68, 0.95, lumistamp);
    fixOverlay();

    if(plotSMoData)
    {
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

        double d2ymax = max(2.2, min(4.5, chdata->GetMaximum(25.0)*1.2));
        TH1 *dummy2 = new TH1F("dummy2", "dummy2", 1000, datahist.hist->GetBinLowEdge(1), datahist.hist->GetBinLowEdge(datahist.hist->GetNbinsX()) + datahist.hist->GetBinWidth(datahist.hist->GetNbinsX()));
        dummy2->GetXaxis()->SetTitle(xaxislabel.c_str());
        dummy2->GetXaxis()->SetTitleOffset(1.05);
        dummy2->GetYaxis()->SetRangeUser(0, d2ymax);
        dummy2->GetYaxis()->SetNdivisions(3, 5, 0);
        dummy2->GetYaxis()->SetTitle("Data/SM");
        dummy2->GetYaxis()->SetTitleOffset(0.32);
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
        fline->SetLineColor(kRed);

        dummy2->Draw();

        //if(chsig) chsig->Draw("hist same");

        TExec * setex2 = new TExec("setex2", "gStyle->SetErrorX(0.5)");
        setex2->Draw();

        TH1 **tgs = new TH1*[systematics.size()];
        int itg = 0;
        const int ebColors[] = {kYellow, kGreen}, NEBCOLORS = sizeof(ebColors) / sizeof(int);
        for(std::vector< std::vector<float> >::const_iterator sit = systematics.begin(); sit != systematics.end(); ++sit)
        {
            char hname[128];
            sprintf(hname, "hsyst_%d", itg);
            tgs[itg] = new TH1F(hname, hname, sizeof(bins) / sizeof(double) - 1, bins);
            for(int i = 1; i <= tgs[itg]->GetNbinsX(); i++)
            {
                tgs[itg]->SetBinContent(i, 1.0);
                if(i - 1 < (int)sit->size()) tgs[itg]->SetBinError(i, fabs(sit->at(i - 1)));
                else if(sit->size() > 0) tgs[itg]->SetBinError(i, fabs(sit->back()));
            }
            tgs[itg]->SetFillColor(ebColors[itg % NEBCOLORS]);
            tgs[itg]->SetMarkerStyle(0);
            tgs[itg]->Draw("E2 same");
            itg++;
        }

        if(systematics.size())
        {
            //TLine *sysStartLine = new TLine(600.0, std::max(-0.9, std::min(0.0, 1.0 - 1.2 * tgs[0]->GetBinError(tgs[0]->FindBin(3500)))), 600.0, 600.0);
            TLine *sysStartLine = new TLine(600.0, 0, 600.0, 600.0);
            sysStartLine->SetLineColor(kBlack);
            sysStartLine->SetLineStyle(2);
            sysStartLine->Draw();
        }

        TExec *setex = new TExec("setex", "gStyle->SetErrorX(0.0)");
        setex->Draw();

        fline->Draw("same");
        chdata->Draw("same");

        TLine *tl = new TLine();
        tl->SetLineColor(kBlack);
        tl->SetLineWidth(2);
        for(int i = 1; i <= chdata->GetNbinsX(); i++)
        {
            if(chdata->GetBinCenter(i) < xmin || chdata->GetBinCenter(i) > xmax) continue;
            if(chbg->GetBinContent(i) > 0.0001 && (d2ymax > chdata->GetBinContent(i) + chdata->GetBinError(i)))   tl->DrawLine(chdata->GetBinCenter(i), std::max(0.0, std::min(d2ymax, chdata->GetBinContent(i) + chdata->GetBinError(i))), chdata->GetBinCenter(i), std::max(0.0, chdata->GetBinContent(i) - chdata->GetBinError(i)));
            else if((chbg->GetBinContent(i) < 0.0001) && (datahist.hist->GetBinContent(i) > 0)) tl->DrawLine(chdata->GetBinCenter(i), 0.0, chdata->GetBinCenter(i), d2ymax);
        }
        fixOverlay();
    }

    if(saveplots)
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

    //if(bghists.size() >= 2)
    //{
    //    for(int ibin = bghists[0].hist->FindBin(600); ibin <= bghists[0].hist->GetNbinsX(); ibin++)
    //    {
    //        printf("%f,", bghists[0].hist->GetBinContent(ibin) / bghists[1].hist->GetBinContent(ibin));
    //    }
    //    printf("\n");
    //    for(int ibin = bghists[0].hist->FindBin(600); ibin <= bghists[0].hist->GetNbinsX(); ibin++)
    //    {
    //        double b1 = bghists[0].hist->GetBinContent(ibin);
    //        double b2 = bghists[1].hist->GetBinContent(ibin);
    //        double e1 = bghists[0].hist->GetBinError(ibin);
    //        double e2 = bghists[1].hist->GetBinError(ibin);
    //        
    //        printf("%f,", (b1 / b2) * sqrt(e1*e1/(b1*b1) + e2*e2/(b2*b2)));
    //    }
    //    printf("\n");
    //}
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
        printf("%s,%s", i->label.c_str(), "sigeff");
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
    string systsd[] = {"gamerr", "norm", "pdf", "fact", "ren", "lumi", "jes", "jer", "mer", "muonid", "trig", "pu", "id", "escale"};
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
    std::vector<std::string> histQs;
    for(size_t pos = 0;;) 
    {
        size_t npos = histValues.find(':', pos + 1);
        histQs.push_back(histValues.substr(pos, npos - pos));
        if(npos == size_t(-1)) break;
        pos = npos + 1;
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
            case 7:
                if(name.find("mWR") < name.size()) *axislabel = "M_{#mu#mujj} [GeV]";
                else if(name.find("mWR_1b") < name.size()) *axislabel = "M_{#mu#mubj} [GeV]";
                else if(name.find("mWR_2b") < name.size()) *axislabel = "M_{#mu#mubb} [GeV]";
                else if(name.find("mLL") < name.size()) *axislabel = "M_{#mu#mu} [GeV]";
                else if(name.find("mLL_1b") < name.size()) *axislabel = "M_{#mu#mu} (1 b-tag) [GeV]";
                else if(name.find("mLL_2b") < name.size()) *axislabel = "M_{#mu#mu} (2 b-tag) [GeV]";
                else if(name.find("mLLZoom") < name.size()) *axislabel = "M_{#mu#mu} [GeV]";
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
                if(name.find("mWR_1b") < name.size()) *axislabel = "M_{eebj} [GeV]";
                else if(name.find("mWR_2b") < name.size()) *axislabel = "M_{eebb} [GeV]";
                else if(name.find("mWR") < name.size()) *axislabel = "M_{eejj} [GeV]";
                else if(name.find("mLL") < name.size()) *axislabel = "M_{ee} [GeV]";
                else if(name.find("mLL_1b") < name.size()) *axislabel = "M_{ee} (1 b-tag) [GeV]";
                else if(name.find("mLL_2b") < name.size()) *axislabel = "M_{ee} (2 b-tag) [GeV]";
                else if(name.find("mLLZoom") < name.size()) *axislabel = "M_{ee} [GeV]";
                //else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{e_{#lower[-0.2]{1}}}} [GeV]";
                //else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{e_{#lower[-0.2]{2}}}} [GeV]";
                else if(name.find("mNuR1") < name.size()) *axislabel = "M_{N_{1}} [GeV]";
                else if(name.find("mNuR2") < name.size()) *axislabel = "M_{N_{2}} [GeV]";
                else if(name.find("ptL1") < name.size()) *axislabel = "p_{T}(e_{1}) [GeV]";
                else if(name.find("ptL2") < name.size()) *axislabel = "p_{T}(e_{2}) [GeV]";
                else if(name.find("etaL1") < name.size()) *axislabel = "#eta(e_{1})";
                else if(name.find("etaL2") < name.size()) *axislabel = "#eta(e_{2})";
                else if(name.find("phiL1") < name.size()) *axislabel = "#phi(e_{1})";
                else if(name.find("phiL2") < name.size()) *axislabel = "#phi(e_{2})";
                else if(name.find("mJJ") < name.size()) *axislabel = "M_{jj} [GeV]";
                else if(name.find("mLQmin") < name.size()) *axislabel = "min M_{LQ} [GeV]";
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
                break;
            case 2:
            case 3:
                if(name.find("mWR") < name.size()) *axislabel = "M_{e#mujj} [GeV]";
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
const static double k_mm_ddtop = /*0.620166*/0.631595,                           k_mm_Zscale = /*1.00858*//*1.02005*/1.027, k_mm_NNLOZ = 1.2036, k_top = 1.13159;
//electron k factors
const static double k_ee_ddtop = /*0.518549*/0.524452 * lumi2012ee / lumi2012mm, k_ee_Zscale = /*0.963943*//*0.939217 0.973471 1.05259*/ 1.00043, k_ee_NNLOZ = 1.1893;

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
const std::string mc_WJ(   "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12_approval_Jun23/prelimWJets.root");

//double lumi2011AMu24 = 216, lumi2011AMu40 = 2056.0, lumi2011B = 2719.0;  //pixel only lumi
//double lumi2011AMu24 = 216.2, lumi2011AMu40 = 1956.7, lumi2011B = 2510.5;  //HF lumi

HeavyNuTree::HNuSlopeFitInfo ll, ul;

void plotMCVar(int cutlevel, std::string plot, int rebin = 5, std::string xaxis = "M_{W_{R}} [GeV]", bool rescale = false, bool log = true)
{

    using namespace std;

    //background legend label, TFile
    std::vector<std::vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgZJ2;
    //bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets madgraph", mc_Z0J, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J, k_mm_NNLOZ / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets",        mc_Z1J, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J, k_mm_NNLOZ / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets",        mc_Z2J, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J, k_mm_NNLOZ / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets",        mc_Z3J, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J, k_mm_NNLOZ / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets",        mc_Z4J, "hNuMu40/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J, k_mm_NNLOZ / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //
    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

    //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    1.0 / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuMu40/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //
    bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets madgraph", mc_Z0J, "hNuE/" + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J, k_ee_NNLOZ / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets", mc_Z1J, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J, k_ee_NNLOZ / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets", mc_Z2J, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J, k_ee_NNLOZ / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets", mc_Z3J, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J, k_ee_NNLOZ / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    bgZJ2.push_back(HnuPlots::FileStruct("Z+Jets", mc_Z4J, "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J, k_ee_NNLOZ / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

    //bgZJ .push_back(HnuPlots::FileStruct("0.8 < mWR < 1.4 TeV", "/local/cms/user/dahmes/forJoe/V03-00-12/bump1TeV/electron-run2012ABCD-V03-00-12.root", "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct("full data", "/local/cms/user/dahmes/forJoe/V03-00-12/electron-run2012ABCD-V03-00-12.root",   "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

    ll.cutlevel = cutlevel;
    ul.cutlevel = 1000;
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, true, 0, 4000, 100, &ll, &ul));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, true, 0, 4000, 100, &ll, &ul));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, true, 0, 4000, 100, &ll, &ul));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, true, 0, 4000, 100, &ll, &ul));
    bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, true, true, 0, 4000, 100, &ll, &ul));

    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, false));

    //bgZJ .push_back(HnuPlots::FileStruct(   "M_{W_{R}} = 2.5 TeV",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_7/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", "hNuE/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "escale",        "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_7/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root", "hNuEescale/" + cutlevels[cutlevel] + "/" + plot, 1.0, 1.0, 1.0,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

    //bgZJ .push_back(HnuPlots::FileStruct(   "M_{#mu#mujj} rms",   "file.root", "rmsm", 1.0, 1.0, 1.0,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
    //bgZJ2.push_back(HnuPlots::FileStruct(   "M_{eejj} rms",       "file.root", "rmse", 1.0, 1.0, 1.0,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));

    bg.push_back(bgZJ2);
    //bg.push_back(bgZJ);
    //bg.push_back(bg3);

    //data
    HnuPlots::FileStruct data("Data", "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12/heavyun_2012Data_2012ABC_mu.root", "hNu/" + cutlevels[cutlevel] + "/" + plot);

    HnuPlots hps(data, bg, sig, 0.0);
    hps.setXAxisTitle(xaxis.c_str());
    hps.setYAxisTitle("Events");
    hps.setLog(log);
    hps.setRebin(rebin);
    hps.plotMCComp(rescale);
}

void setBgandData(int mode, HnuPlots::FileStruct& data, std::vector<std::vector<HnuPlots::FileStruct> >& bg, double& lumi, int cutlevel = 5, std::string plot = "mWR", bool lt = false)
{
    char fdata[128];

    std::vector<HnuPlots::FileStruct> bgTT, bgZJ, bgZ1J, bgZ2J, bgZ3J, bgZ4J, bgOther;

    bool hft = false;
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

    //background
    switch(mode)
    {

        case 0:  //muon plots
            //bgZJ.push_back(HnuPlots::FileStruct("DD Z+Jets",   data_mm,  "hNu/"        + cutlevels[11]    + "/" + plot,     1.0,     1.0,      3.42977226076537911e-02,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 4000, 100, &ll, &ul));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuMu40/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    0.916553*k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    k_mm_Zscale / NZJ,                  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 4000, 100, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
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
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,     1.0 / NtW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW,  1.0 / NtbarW,                                   "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,     1.0 / NZZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,     1.0 / NWZ,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,     1.0 / NWW,                                      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));

            bg.push_back(bgTT);
            bg.push_back(bgZJ);
            bg.push_back(bgOther);
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
            //bgZJ.push_back(HnuPlots::FileStruct(   "DD Z+Jets",   data_ee,  "hNuE/"       + cutlevels[11]          + "/" + plot,     1.0,     1.0,      3.63902344861221985e-02,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    0.89016*k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_Zscale * k_ee_NNLOZ / NZ0J,                 "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
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
            bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,                              "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,                                    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,                                       "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            
            bg.push_back(bgTT);
            bg.push_back(bgZJ);
            bg.push_back(bgOther);
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
            bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bg.push_back(bgTT);
        case 3:
            //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   k_mm_Zscale * k_mm_NNLOZ / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   k_mm_Zscale * k_mm_NNLOZ / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   k_mm_Zscale * k_mm_NNLOZ / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   k_mm_Zscale * k_mm_NNLOZ / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   k_mm_Zscale * k_mm_NNLOZ / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    1.0 / NtW,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 1.0 / NtbarW, "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    1.0 / NZZ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    1.0 / NWZ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuEMu/" + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    1.0 / NWW,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
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
            break;
        case 4:
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,    k_mm_Zscale * k_mm_NNLOZ / NZ0J / 3.44921630331921220e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,    k_mm_Zscale * k_mm_NNLOZ / NZ1J / 3.44921630331921220e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,    k_mm_Zscale * k_mm_NNLOZ / NZ2J / 3.44921630331921220e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,    k_mm_Zscale * k_mm_NNLOZ / NZ3J / 3.44921630331921220e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,    k_mm_Zscale * k_mm_NNLOZ / NZ4J / 3.44921630331921220e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bg.push_back(bgZJ);
            
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
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,    k_ee_Zscale * k_ee_NNLOZ / NZ0J / 3.64878416472316086e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,    k_ee_Zscale * k_ee_NNLOZ / NZ1J / 3.64878416472316086e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,    k_ee_Zscale * k_ee_NNLOZ / NZ2J / 3.64878416472316086e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,    k_ee_Zscale * k_ee_NNLOZ / NZ3J / 3.64878416472316086e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"    + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,    k_ee_Zscale * k_ee_NNLOZ / NZ4J / 3.64878416472316086e-02,                "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt, hft, 0, 0, -1, &ll, &ul));
            bg.push_back(bgZJ);
            
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
            bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   k_mm_Zscale / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZ1J.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   k_mm_Zscale / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZ2J.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   k_mm_Zscale / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZ3J.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   k_mm_Zscale / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
            bgZ4J.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   k_mm_Zscale / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));

            //bgTT.push_back(HnuPlots::FileStruct("t#bar{t}", "/local/cms/user/pastika/heavyNuAnalysis_2012/Summer12/heavynu_2011Bg_mumu_test_heavyNuAnalysis_TTBar_Skim.root", "hNuMu40/" + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, 225.197,  69620.0 / 6736135 * 1.5, "hNuMu40/cutlevel"));
            sprintf(fdata, "%s", mc_ZJ.c_str());
            lumi += lumi2012mm;
            data.histpath = "hNuMu40/" + cutlevels[cutlevel] + "/" + plot;

            bg.push_back(bgZJ);
            bg.push_back(bgZ1J);
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

void plot2012(int mode = 0, int cutlevel = 5, std::string plot = "mWR", int rebin = 5, bool log = true, double xmin = 0.0, double xmax = 3500.0, bool autoY = true)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig, vsig2, vsig3, vsig4, vsig5, vsig6;
    setBgandData(mode, data, bg, lumi, cutlevel, plot, true);

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
    //vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.0 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall_11/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 0.013339, 1.214, normhist, 0.0, 0.0, true, signormbin));
    //vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 1.1 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1100_MNu-550_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 0.013339, 1.214 * 0.5 * 0.75, normhist, 0.0, 0.0, true, signormbin));
    //vsig2.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 1.0 TeV",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1000_MNu-500_TuneZ2star_8TeV-pythia6-tauola.root",    "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.667875, 1.340 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin));
    //vsig3.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 1.5 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-1500_MNu-750_TuneZ2star_8TeV-pythia6-tauola.root",    "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.082688, 1.293 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin));
    //vsig4.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 2.0 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.013339, 1.214 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin));
    //vsig5.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 2.5 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2500_MNu-1250_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.002286, 1.140 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin));
    //vsig6.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}}(N_{#tau}) = 3.0 TeV ",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-3000_MNu-1500_TuneZ2star_8TeV-pythia6-tauola.root",   "hTauX/" + cutlevels[cutlevel] + "/" + plot, lumi, 0.000393, 1.151 * 0.062, "hTauX/mc_type", 0.0, 0.0, true, signormbin));
    //if(mode <= 1) sig.push_back(vsig);
    //else if(mode == 2)
    //{
    //    sig.push_back(vsig2);
    //    sig.push_back(vsig3);
    //    sig.push_back(vsig4);
    //    sig.push_back(vsig5);
    //    sig.push_back(vsig6);
    //}

    HnuPlots hps(data, bg, sig, lumi);
    switch(mode)
    {
        case 7:
            rebin = -1;
            xmax = 4000;
        case 0:
            if(rebin > 0) hps.setFormLabel("hNu_mm_2012");
            else          hps.setFormLabel("hNu_mm_ls_2012");
            hps.setSavePlots(true);
            if(!plot.compare("mWR"))
            {
                hps.loadSystFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/systematicsdb_mu_2012.csv", "/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/ratesdb.csv", (mode == 7));
                hps.mcBgShape();
            }
            break;
        case 8:
            rebin = -1;
            xmax = 4000;
        case 1:
            if(rebin > 0) hps.setFormLabel("hNu_ee_2012");
            else          hps.setFormLabel("hNu_ee_ls_2012");
            hps.setSavePlots(true);
            if(!plot.compare("mWR"))
            {

                hps.loadSystFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/systematicsdb_elec_2012.csv", "/home/ugrad/pastika/cms/HeavyNu/CMSSW_6_1_1/src/HeavyNu/Limits/ctool/ratesdb_elec.csv", (mode == 8));
                hps.mcBgShape();
            }
            break;
        case 2:
            hps.setFormLabel("hNu_em_2012");
            hps.setSavePlots(true);
            break;
        case 3:
            hps.setFormLabel("hNu_em2_2012");
            hps.setSavePlots(true);
            hps.setCompPlot(false);
            break;
        case 4:
        case 5:
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
    if(isMuon) hps.setFormLabel(flabel + "_mm_2012");

    else hps.setFormLabel(flabel + "_ee_2012");
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
        //bgOther.push_back(HnuPlots::FileStruct("QCD",      data_mm,  "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, 1.0,     1.0,           1.0       ,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    -1.0 / NZJ,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectW,    -1.0 / NtW,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsectbarW, -1.0 / NtbarW,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZZ,    -1.0 / NZZ,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0 / NWZ,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWW,    -1.0 / NWW,          "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WJ,    "hNuMu1QCD/"  + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecWZ,    -1.0, "hNuMu1QCD/cutlevel", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    data_mm,  "hNuFakeMuGoodEwgtMu/" + cutlevelsTop[cutlevel] + "/" + plot, 1.0,     1.0, -k_mm_ddtop * 0.11,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));

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
        //bgOther.push_back(HnuPlots::FileStruct("QCD",      data_ee,  "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, 1.0,     1.0,           1.0      ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZJ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    -1.0 / NZJ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsectW,    -1.0 / NtW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsectbarW, -1.0 / NtbarW,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZZ,    -1.0 / NZZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWZ,    -1.0 / NWZ,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWW,    -1.0 / NWW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WJ,    "hNuE1QCD/"   + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecWZ,    -1.0, "hNuE1QCD/cutlevel", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));
        //bgOther.push_back(HnuPlots::FileStruct("Other",    data_mm,  "hNuGoodMuFakeEwgtE/" + cutlevelsTop[cutlevel] + "/" + plot, 1.0,     1.0,  -k_ee_ddtop * 0.03,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0, lt));

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
    if(isMuon) hps.setFormLabel(flabel + "mm_2012" + flend);
    else       hps.setFormLabel(flabel + "ee_2012" + flend);
    hps.setSavePlots(true);
    hps.setXRange(60.0, 500.0);
    hps.scaleByShape(60, 200, 2);
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
    if(isMuon) hps.setFormLabel("ttMCnorm_mm_2012");
    else hps.setFormLabel("ttMCnorm_ee_2012");
    hps.setSavePlots(true);
    hps.plotNorm(40.0, 6000.0);
}

void plotZJNorm(bool isMuon = true, int cutlevel = 4, bool log = true)//, std::string sample = "")
{
    using namespace std;

    char plot[] = "mLLZoom", fdata[256];
    double lumi = 0.0;
    string datahistname;

    //data
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> bgZJ, bgOther, bgTT;
    if(isMuon)
    {
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    1.0 / NZJ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuMu40/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012mm, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J, k_mm_NNLOZ  / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J, k_mm_NNLOZ  / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J, k_mm_NNLOZ  / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J, k_mm_NNLOZ  / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J, k_mm_NNLOZ  / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,         1.0,     1.0,  k_mm_ddtop,             "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZJ,    - k_mm_ddtop / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ0J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ0J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ1J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ1J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ2J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ2J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ3J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ3J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZ4J,   - k_mm_ddtop * k_mm_Zscale * k_mm_NNLOZ / NZ4J,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    - k_mm_ddtop / NtW,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, - k_mm_ddtop / NtbarW,  "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    - k_mm_ddtop / NZZ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    - k_mm_ddtop / NWZ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    - k_mm_ddtop / NWW,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tt,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, k_top / Nttbar,         "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectW,    1.0 / NtW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsectbarW, 1.0 / NtbarW,           "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecZZ,    1.0 / NZZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWZ,    1.0 / NWZ,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuMu40/"    + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecWW,    1.0 / NWW,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
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
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets Sherpa",   "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_5/heavynu_2012Bg_DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa_START53_V7C-v2.root",    "hNuE/"       + cutlevels[cutlevel]    + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_ZJ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    k_ee_Zscale / NZJ,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z0J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   k_ee_NNLOZ  / NZ0J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z1J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   k_ee_NNLOZ  / NZ1J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z2J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   k_ee_NNLOZ  / NZ2J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z3J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   k_ee_NNLOZ  / NZ3J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgZJ.push_back(HnuPlots::FileStruct(   "Z+Jets",   mc_Z4J,   "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   k_ee_NNLOZ  / NZ4J, "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct("t#bar{t} (MC)", mc_tt,  "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012mm, xsecttbar, 1.0 / Nttbar,       "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", data_em,  "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot,        1.0,     1.0,   k_ee_ddtop,              "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        //bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZJ,    "hNuEMu/"     + cutlevelsTop[cutlevel] + "/" + plot, lumi2012ee, xsecZJ,    - k_ee_ddtop / NZJ,    "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z0J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ0J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ0J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z1J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ1J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ1J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z2J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ2J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ2J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z3J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ3J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ3J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_Z4J,   "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZ4J,   - k_ee_ddtop * k_ee_Zscale * k_ee_NNLOZ / NZ4J,     "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    - k_ee_ddtop / NtW,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_tbarW, "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, - k_ee_ddtop / NtbarW,   "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_ZZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    - k_ee_ddtop / NZZ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WZ,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    - k_ee_ddtop / NWZ,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgTT.push_back(HnuPlots::FileStruct(   "t#bar{t}", mc_WW,    "hNuEMu/"     + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    - k_ee_ddtop / NWW,      "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectW,    1.0 / NtW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_tbarW, "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsectbarW, 1.0 / NtbarW,            "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_ZZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecZZ,    1.0 / NZZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WZ,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWZ,    1.0 / NWZ,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
        bgOther.push_back(HnuPlots::FileStruct("Other",    mc_WW,    "hNuE/"       + cutlevels[cutlevel] + "/" + plot, lumi2012ee, xsecWW,    1.0 / NWW,               "", 0.0, 0.0, true, 1, true, 0.0, 0.0));
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
    data.label = "Data";
    data.file = fdata;

    HnuPlots hps(data, bg, sig, lumi);
    hps.autoSetHistogramAxisTitle(!isMuon);
    hps.setYAxisTitle("please auto set the axis");
    hps.setRebin(20);
    hps.setLog(log);
    if(isMuon) hps.setFormLabel("zjnorm_mm_2012");

    else hps.setFormLabel("zjnorm_ee_2012");
    hps.setSavePlots(true);
    hps.setXRange(60.0, 500.0);
    hps.plotNorm(60.0, 120.0);
}

void plotCutFlow(int mode = 0)
{
    using namespace std;

    double lumi = 0.0;
    HnuPlots::FileStruct data;

    //background legend label, TFile
    vector<vector<HnuPlots::FileStruct> > bg, sig;
    vector<HnuPlots::FileStruct> vsig;
    setBgandData(mode, data, bg, lumi, 17, "cutlevel");

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
    vsig.push_back(HnuPlots::FileStruct("M_{#lower[-0.1]{W_{#lower[-0.2]{R}}}} = 2.0 TeV",  "/local/cms/user/pastika/heavyNuAnalysis_2012/Fall12_4/heavynu_2012Bg_WRToNuLeptonToLLJJ_MW-2000_MNu-1000_TuneZ2star_8TeV-pythia6-tauola.root",  histograms, lumi, 0.013339, 1.214, normhist, 0.0, 0.0, true, signormbin));
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
    //setBgandData(true, true data, bg, lumi, 9, "cutlevel");

    std::string histograms = "", normhist = "", label = "";
    int signormbin = 0;

    switch(mode)
    {

        case 0:
            histograms = "hNuMu40/cut5_diLmass/mWR";
            normhist = "hNuMu40/mc_type";
            label = "_mm";
            signormbin = 3;
            lumi = lumi2012mm;
            break;
        case 1:
            histograms = "hNuE/cut5_diLmass/mWR";
            normhist = "hNuE/mc_type";
            label = "_ee";
            signormbin = 2;
            lumi = lumi2012ee;
            break;
        case 2:
            histograms = "hNuEMu/cut5_diLmass/mWR";
            normhist = "hNuEMu/mc_type";
            label = "_em";
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
    plot2012(0, 4, "mWR");
    plot2012(0, 5, "mWR");
    plot2012(0, 6, "mWR");
    plot2012(0, 4, "mLL");
    plot2012(0, 5, "mLL");
    plot2012(0, 6, "mLL");
    plot2012(0, 4, "mJJ");
    plot2012(0, 5, "mJJ");
    plot2012(0, 6, "mJJ");
    plot2012(0, 4, "mNuR1");
    plot2012(0, 5, "mNuR1");
    plot2012(0, 6, "mNuR1");
    plot2012(0, 4, "mNuR2");
    plot2012(0, 5, "mNuR2");
    plot2012(0, 6, "mNuR2");
    plot2012(0, 4, "ptL1");
    plot2012(0, 5, "ptL1");
    plot2012(0, 6, "ptL1");
    plot2012(0, 4, "ptL2");
    plot2012(0, 5, "ptL2");
    plot2012(0, 6, "ptL2");
    plot2012(0, 4, "ptJ1");
    plot2012(0, 5, "ptJ1");
    plot2012(0, 6, "ptJ1");
    plot2012(0, 4, "ptJ2");
    plot2012(0, 5, "ptJ2");
    plot2012(0, 6, "ptJ2");
    plot2012(0, 5, "mWR", -1, true, 0.0, 4000);
    plot2012(1, 4, "mWR");
    plot2012(1, 5, "mWR");
    plot2012(1, 6, "mWR");
    plot2012(1, 4, "mLL");
    plot2012(1, 5, "mLL");
    plot2012(1, 6, "mLL");
    plot2012(1, 4, "mJJ");
    plot2012(1, 5, "mJJ");
    plot2012(1, 6, "mJJ");
    plot2012(1, 4, "mNuR1");
    plot2012(1, 5, "mNuR1");
    plot2012(1, 6, "mNuR1");
    plot2012(1, 4, "mNuR2");
    plot2012(1, 5, "mNuR2");
    plot2012(1, 6, "mNuR2");
    plot2012(1, 4, "ptL1");
    plot2012(1, 5, "ptL1");
    plot2012(1, 6, "ptL1");
    plot2012(1, 4, "ptL2");
    plot2012(1, 5, "ptL2");
    plot2012(1, 6, "ptL2");
    plot2012(1, 4, "ptJ1");
    plot2012(1, 5, "ptJ1");
    plot2012(1, 6, "ptJ1");
    plot2012(1, 4, "ptJ2");
    plot2012(1, 5, "ptJ2");
    plot2012(1, 6, "ptJ2");
    plot2012(1, 4, "etaL1", 5, true, -3, 3);
    plot2012(1, 5, "etaL1", 5, true, -3, 3);
    plot2012(1, 6, "etaL1", 5, true, -3, 3);
    plot2012(1, 4, "etaL2", 5, true, -3, 3);
    plot2012(1, 5, "etaL2", 5, true, -3, 3);
    plot2012(1, 6, "etaL2", 5, true, -3, 3);
    plot2012(1, 4, "etaJ1", 5, true, -3, 3);
    plot2012(1, 5, "etaJ1", 5, true, -3, 3);
    plot2012(1, 6, "etaJ1", 5, true, -3, 3);
    plot2012(1, 4, "etaJ2", 5, true, -3, 3);
    plot2012(1, 5, "etaJ2", 5, true, -3, 3);
    plot2012(1, 6, "etaJ2", 5, true, -3, 3);
    plot2012(1, 4, "phiL1", 2, true, -4, 4);
    plot2012(1, 5, "phiL1", 2, true, -4, 4);
    plot2012(1, 6, "phiL1", 2, true, -4, 4);
    plot2012(1, 4, "phiL2", 2, true, -4, 4);
    plot2012(1, 5, "phiL2", 2, true, -4, 4);
    plot2012(1, 6, "phiL2", 2, true, -4, 4);
    plot2012(1, 4, "phiJ1", 2, true, -4, 4);
    plot2012(1, 5, "phiJ1", 2, true, -4, 4);
    plot2012(1, 6, "phiJ1", 2, true, -4, 4);
    plot2012(1, 4, "phiJ2", 2, true, -4, 4);
    plot2012(1, 5, "phiJ2", 2, true, -4, 4);
    plot2012(1, 6, "phiJ2", 2, true, -4, 4);
    plot2012(1, 5, "mWR", -1, true, 0.0, 4000);

    plot2012(2, 5, "mWR");
    plot2012(3, 5, "mWR");

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
