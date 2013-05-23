//root -b -q -l upperLimit.C+
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
#include "TAxis.h"

#include "tdrstyle.C"
#include "../acceptance/acceptance_db.C"

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <bits/stl_vector.h>
#include "stdio.h"

struct XYZ
{
    int x, y;
    float z;

    bool operator ==(const XYZ & v) const
    {
        return(x == v.x) && (y == v.y);
    }
} ;

static const double lumi2010 = 0;
static const double lumi2011 = 4971.0;
static const double lumi2012 = 3585.0;

AcceptanceDB* db, *db2011;

class compareByNu
{
public:

    bool operator() (const XYZ& a, const XYZ& b)
    {
        return a.y > b.y;
    }
} ;

class compareByNu2
{
public:

    bool operator() (const XYZ& a, const XYZ& b)
    {
        return a.y < b.y;
    }
} ;

std::vector<XYZ> getSortedVector(int iwr, std::vector<XYZ> a)
{

    // First create a vector that contains only the points we want
    std::vector<XYZ> vec;
    for(unsigned int i = 0; i < a.size(); i++)
        if(a.at(i).x == iwr) vec.push_back(a.at(i));
    std::sort(vec.begin(), vec.end(), compareByNu());
    return vec;
}

std::vector<XYZ> getSortedVector(int iwr, std::vector<XYZ> a, std::vector<XYZ> b)
{

    std::vector<XYZ> vec;
    for(unsigned int i = 0; i < a.size(); i++)
        if(a.at(i).x == iwr) vec.push_back(a.at(i));
    for(unsigned int i = 0; i < b.size(); i++)
        if(b.at(i).x == iwr) vec.push_back(b.at(i));
    return getSortedVector(iwr, vec);
}

double getWgtAcceptance(int imwr, int imnu, int mode)
{
    switch(mode)
    {
        case 0:
        case 1:
        case 2:
            return db->getBestEstimate(imwr, imnu, 2012);
        case 3:
            //horrible hack here to fix lack of any > 3TeV 2011 acceptance 
            if(imwr <= 3000) return (lumi2011 * db2011->getBestEstimate(imwr, imnu, 2011) + lumi2012 * db->getBestEstimate(imwr, imnu, 2012)) / (lumi2011 + lumi2012);
            else return db->getBestEstimate(imwr, imnu, 2012);
    }
    return 0;
}

double getEfficiency(int imwr, std::vector<XYZ> v, int mode)
{
    int minDistLo = 10000;
    int minDistHi = 10000;
    unsigned int nearestLo = v.size() + 1;
    unsigned int nearestHi = v.size() + 1;
    for(unsigned int i = 0; i < v.size(); i++)
    {
        int dist = fabs(imwr - v.at(i).x);
        if(dist < minDistLo && v.at(i).x <= imwr)
        {
            nearestLo = i;
            minDistLo = dist;
        }
        if(dist < minDistHi && v.at(i).x >= imwr)
        {
            nearestHi = i;
            minDistHi = dist;
        }
    }
    // Now for efficiency, take the Nu mass closest to the midpoint
    int minNuDistLo = 10000;
    int minNuDistHi = 10000;
    for(unsigned int i = 0; i < v.size(); i++)
    {
        if(fabs(imwr - v.at(i).x) == minDistLo && v.at(i).x <= imwr)
        {
            int dist = fabs(imwr / 2 - v.at(i).y);
            if(dist < minNuDistLo)
            {
                nearestLo = i;
                minNuDistLo = dist;
            }
        }
        if(fabs(imwr - v.at(i).x) == minDistHi && v.at(i).x >= imwr)
        {
            int dist = fabs(imwr / 2 - v.at(i).y);
            if(dist < minNuDistHi)
            {
                nearestHi = i;
                minNuDistHi = dist;
            }
        }
    }

    //     std::cout << "Looking for point: " << imwr << " and found nearest points "
    //               << v.at(nearestLo).x << "," << v.at(nearestLo).y
    //               << " and " << v.at(nearestHi).x << "," << v.at(nearestHi).y
    //               << ": " << nearestLo << "," << nearestHi << std::endl ; 

    double eff = -1.;
    if(nearestHi > v.size()) std::cout << "NONONO " << imwr << std::endl;
    if(nearestLo != nearestHi)
    { // Not a reference point
        double wgtLoAcc = getWgtAcceptance(v.at(nearestLo).x, v.at(nearestLo).y, mode);
        double wgtHiAcc = getWgtAcceptance(v.at(nearestHi).x, v.at(nearestHi).y, mode);

        double loEff = v.at(nearestLo).z * wgtLoAcc;
        double hiEff = v.at(nearestHi).z * wgtHiAcc;
        eff = interpol1d(imwr, v.at(nearestLo).x, loEff, v.at(nearestHi).x, hiEff);
    }
    else
    {
        if(nearestLo > v.size()) std::cout << "Warning: Could not find point" << std::endl;
        else
        {
            // std::cout << "At a reference point" << std::endl ; 
            double wgtAcc = getWgtAcceptance(v.at(nearestLo).x, v.at(nearestLo).y, mode);
            // std::cout << "wgtAcc is " << wgtAcc << std::endl ; 
            eff = v.at(nearestLo).z * wgtAcc;
            // std::cout << "eff: " << eff << std::endl ; 
        }
    }
    return eff;
}

void plotLimits(int mode = 0, int minval = -1, int maxval = -1, float xmini = 1000.0, float xmaxi = 3100.0, float ymaxi = 2000.0, std::string limitFile = "", std::string csFile = "cs.txt", std::string csFile2 = "xsecs.txt")
{
    if(minval < 0) minval = 1000;
    if(maxval < 0) maxval = 3000;

    std::string file("");
    if(limitFile.length() == 0)
    {
        switch(mode)
        {
            case 0:
                file = "mu12_jun26_smoothed.csv";
                break;
            case 1:
                file = "el12_jun26_smoothed.csv";
                break;
            case 2:
                file = "em12_jun26_smoothed.csv";
                break;
            case 3:
                file = "mu2y_smoothed.csv";
                break;
        }
        limitFile = file;
    }

    const int MWR_STEP = 100, MNU_STEP = 10;
    const int TRANSITION_MWR = 2500;
    const float xmin = xmini, xmax = xmaxi, ymax = ymaxi;

    db = new AcceptanceDB("2012");
    db2011 = new AcceptanceDB("Muon");

    // Inizialization
    //FILE *Ptr;
    char *str1 = new char[128];
    FILE *inLimit = fopen(limitFile.c_str(), "r");
    FILE *inCS = fopen(csFile.c_str(), "r");
    FILE *inCS2 = fopen(csFile2.c_str(), "r");

    //std::string outputFileName = "obsexpSummary.txt";
    //Ptr = fopen(outputFileName.c_str(), "a");

    std::vector<XYZ> obsLimit;
    std::vector<XYZ> expLimit;
    std::vector<XYZ> exp68loLimit;
    std::vector<XYZ> exp68hiLimit;
    std::vector<XYZ> exp95loLimit;
    std::vector<XYZ> exp95hiLimit;

    std::map<std::pair<int, int>, float > ps;
    std::map<std::pair<int, int>, float > csforsf;

    char buff[4096];
    char *c;
    while(!feof(inCS) && (c = fgets(buff, 4096, inCS)) != NULL)
    {
        int energy, tmp_mWR, tmp_mNu;
        float tmp_cs, tmp_sf;
        for(char* k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
        if(sscanf(buff, "%d %d %d %f %f\n", &energy, &tmp_mWR, &tmp_mNu, &tmp_cs, &tmp_sf) == 5 && energy == 8)
        {
            csforsf[std::make_pair(tmp_mWR, tmp_mNu)] = tmp_cs * tmp_sf;
        }
    }
    fclose(inCS);

    while(!feof(inCS2) && (c = fgets(buff, 4096, inCS2)) != NULL)
    {
        int energy, tmp_mWR, tmp_mNu;
        float tmp_cs;
        for(char* k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
        if(sscanf(buff, "%d %d %d %f\n", &energy, &tmp_mWR, &tmp_mNu, &tmp_cs) == 4 && energy == 8)
        {
            if(!(tmp_mWR % MWR_STEP))
            {
                ps[std::make_pair(tmp_mWR, tmp_mNu)] = tmp_cs * 1000000000.0;  //think about the units!
                //std::cout << tmp_mWR << "\t" << tmp_mNu << std::endl;
            }
        }
    }
    fclose(inCS2);

    float nFlavorScale = 1.0;
    switch(mode)
    {
        case 0:
            nFlavorScale = 1.0;
            break;
        case 1:
            nFlavorScale = 1.0;
            break;
        case 2:
            nFlavorScale = 8.0 / 3.0;
            break;
        case 3:
            nFlavorScale = 1.0;
            break;
    }
    // Build up cross-section as a function of W_R and N mass
    TGraph2D* grCS = new TGraph2D(ps.size());
    grCS->SetNameTitle("xs", "xs;WR;nuR;#sigma");
    int ibin = 0;
    for(std::map<std::pair<int, int>, float >::const_iterator it = ps.begin(); it != ps.end(); ++it)
    {

        std::pair<int, int> mid = std::make_pair(it->first.first, it->first.first / 2);
        double midNLOcs = csforsf[mid], midLOcs = ps[mid];
        switch(mode)
        {
            case 0:
            case 1:
            case 2:
                grCS->SetPoint(ibin++, it->first.first, it->first.second, it->second * nFlavorScale * midNLOcs / midLOcs);
                break;
            case 3:
                grCS->SetPoint(ibin++, it->first.first, it->first.second, (it->second * nFlavorScale / midLOcs));
                break;
        }
    }
    
    double jscale = 1.0;
    switch(mode)
        {
            case 0:
            case 1:
            case 2:
                jscale = 1.0;
                break;
            case 3:
                jscale = 1.0;  //the combined limits are given in ratio to expectation
                break;
        }

    while(!feof(inLimit) && (c = fgets(buff, 4096, inLimit)) != NULL)
    {
        float tmp_mwr, tmp_mnu, f3, f4, f5, f6, f7, f8;
        for(char* k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
        if(sscanf(buff, "%f %f %f %f %f %f %f %f\n", &tmp_mwr, &tmp_mnu, &f3, &f4, &f5, &f6, &f7, &f8) == 8)
        {
            //printf("%f %f %f %f %f %f %f %f\n", tmp_mwr, tmp_mnu, f3, f4, f5, f6, f7, f8);
            const XYZ p = {tmp_mwr, tmp_mnu, f3 / jscale};  
            const XYZ q = {tmp_mwr, tmp_mnu, f4 / jscale};
            const XYZ h68 = {tmp_mwr, tmp_mnu, f7 / jscale};
            const XYZ l68 = {tmp_mwr, tmp_mnu, f6 / jscale};
            const XYZ h95 = {tmp_mwr, tmp_mnu, f8 / jscale};
            const XYZ l95 = {tmp_mwr, tmp_mnu, f5 / jscale};
            if(find(obsLimit.begin(), obsLimit.end(), p) == obsLimit.end()) obsLimit.push_back(p);
            if(find(expLimit.begin(), expLimit.end(), q) == expLimit.end()) expLimit.push_back(q);
            if(find(exp68hiLimit.begin(), exp68hiLimit.end(), h68) == exp68hiLimit.end()) exp68hiLimit.push_back(h68);
            if(find(exp68loLimit.begin(), exp68loLimit.end(), l68) == exp68loLimit.end()) exp68loLimit.push_back(l68);
            if(find(exp95hiLimit.begin(), exp95hiLimit.end(), h95) == exp95hiLimit.end()) exp95hiLimit.push_back(h95);
            if(find(exp95loLimit.begin(), exp95loLimit.end(), l95) == exp95loLimit.end()) exp95loLimit.push_back(l95);
        }
    }
    fclose(inLimit);

    std::map<std::pair<int, int>, double> rawObsPts;
    std::map<std::pair<int, int>, double> rawExpPts;

    for(int imwr = minval; imwr <= maxval; imwr += MWR_STEP)
    {

        if(imwr % 100 == 0) std::cout << "Computing for W_R mass : " << imwr << std::endl;

        double eff_wr = getEfficiency(imwr, obsLimit, mode);
        for(int imnu = 10; imnu < imwr; imnu += MNU_STEP)
        {
            double acc = getWgtAcceptance(imwr, imnu, mode);
            if(acc < 0) continue;
            rawObsPts[std::make_pair(imwr, imnu)] = eff_wr / acc;
            //std::cout << imwr << "\t" << imnu << "\t" << eff_wr / acc << std::endl;
        }

        eff_wr = getEfficiency(imwr, expLimit, mode);
        for(int imnu = 10; imnu < imwr; imnu += MNU_STEP)
        {
            double acc = getWgtAcceptance(imwr, imnu, mode);
            if(acc < 0) continue;
            rawExpPts[std::make_pair(imwr, imnu)] = eff_wr / acc;
            //std::cout << imwr << "\t" << imnu << "\t" << eff_wr / acc << std::endl;
        }
    }

    //
    // 1D plots
    //
    for(int iwr = minval; iwr <= maxval; iwr += MWR_STEP)
    {

        TGraph* h_exp = new TGraph();
        TGraph* h_obs = new TGraph();
        TGraph* h_theoryL = new TGraph();
        TGraphAsymmErrors* h_theory = new TGraphAsymmErrors();
        TGraphAsymmErrors* h_exp1sig = new TGraphAsymmErrors();
        TGraphAsymmErrors* h_exp2sig = new TGraphAsymmErrors();

        double eff_68lo = getEfficiency(iwr, exp68loLimit, mode);
        double eff_68hi = getEfficiency(iwr, exp68hiLimit, mode);
        double eff_95lo = getEfficiency(iwr, exp95loLimit, mode);
        double eff_95hi = getEfficiency(iwr, exp95hiLimit, mode);

        unsigned int c0 = 0, c1 = 0, c2 = 0;
        std::map<std::pair<int, int>, double >::const_iterator imp;
        for(int inu = 50; inu < iwr - 50; inu += MNU_STEP)
        {
            std::pair<int,  int> i(iwr, inu);
            if((imp = rawExpPts.find(i)) != rawExpPts.end())
            {
                double acc = getWgtAcceptance(iwr, inu, mode);
                if(acc < 0) continue;
                double d68lo = imp->second - eff_68lo / acc;
                double d68hi = eff_68hi / acc - imp->second;
                double d95lo = imp->second - eff_95lo / acc;
                double d95hi = eff_95hi / acc - imp->second;

                h_exp->SetPoint(c1, (double)inu, imp->second);
                h_exp1sig->SetPoint(c1, (double)inu, imp->second);
                h_exp2sig->SetPoint(c1, (double)inu, imp->second);

                h_exp1sig->SetPointError(c1, 0., 0., d68lo, d68hi);
                h_exp2sig->SetPointError(c1++, 0., 0., d95lo, d95hi);
            }
            if((imp = rawObsPts.find(std::make_pair(iwr, inu))) != rawObsPts.end())
                h_obs->SetPoint(c2++, (double)inu, imp->second);

            double sigmaBR = grCS->Interpolate(iwr, inu);
            // if ( sigmaBR < 0.000001 ) continue ;
            if(sigmaBR > 1e-5 && sigmaBR < 1000.0 && inu > 10)
            {
                h_theoryL->SetPoint(c0, inu, sigmaBR);
                h_theory->SetPoint(c0, inu, sigmaBR);
                h_theory->SetPointError(c0++, 0., 0., sigmaBR * 0.1, sigmaBR * 0.1);
            }
        }

        // Style for lines
        h_exp->SetLineWidth(3);
        h_exp->SetLineStyle(2);
        h_obs->SetLineWidth(3);
        h_obs->SetLineColor(1);
        h_obs->SetLineStyle(1);
        h_exp1sig->SetFillColor(kBlue - 7);
        h_exp2sig->SetFillColor(kYellow);
        h_theoryL->SetLineColor(kRed);
        h_theoryL->SetLineWidth(3);
        h_theory->SetFillColor(kGreen);

        sprintf(str1, "csUL%i", iwr);

        TCanvas* csUL = new TCanvas(str1, "c1", 600, 600);
        setTDRStyle();
        fixOverlay();
        csUL->SetLeftMargin(0.19);
        csUL->SetRightMargin(0.06);
        csUL->SetTopMargin(0.06);
        csUL->SetLogy(1);

        TMultiGraph* mg = new TMultiGraph();
        mg->SetName("");
        switch(mode)
        {
            case 0:
                mg->SetTitle(";M_{N_{#mu}} [GeV];#sigma(pp #rightarrow W_{R})#times Br(W_{R} #rightarrow #mu#mujj) [pb]");
                break;
            case 1:
                mg->SetTitle(";M_{N_{e}} [GeV];#sigma(pp #rightarrow W_{R})#times Br(W_{R} #rightarrow eejj) [pb]");
                break;
            case 2:
                mg->SetTitle(";M_{N_{e,#mu,#tau}} [GeV];#sigma(pp #rightarrow W_{R})#times Br(W_{R} #rightarrow (ee+#mu#mu+#tau#tau)jj) [pb]");
                break;
            case 3:
                mg->SetTitle(";M_{N_{#mu}} [GeV];Arbitrary scale");
                break;
        }

        mg->Add(h_exp2sig, "E3");
        mg->Add(h_exp1sig, "E3");
        //mg->Add(h_theory, "E3");
        mg->Add(h_theoryL, "L");
        mg->Add(h_exp, "L");
        mg->Add(h_obs, "L");

        TH1D* dummy = new TH1D("dummy", "dummy", 100, 0, double(iwr));
        switch(mode)
        {
            case 0:
                dummy->GetXaxis()->SetTitle("M_{N_{#mu}} [GeV]");
                dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow #mu#mujj) [pb]");
                break;
            case 1:
                dummy->GetXaxis()->SetTitle("M_{N_{e}} [GeV]");
                dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow eejj) [pb]");
                break;
            case 2:
                dummy->GetXaxis()->SetTitle("M_{N_{e,#mu,#tau}} [GeV]");
                dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow (ee+#mu#mu+#tau#tau)jj) [pb]");
                break;
            case 3:
                dummy->GetXaxis()->SetTitle("M_{N_{#mu}} [GeV]");
                dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow #mu#mujj) [pb]");
                break;
        }
        dummy->GetYaxis()->SetTitleOffset(1.45);
        // dummy->SetMinimum(minyval) ; 
        // dummy->SetMaximum(maxyval) ; 

        mg->Draw("AH");
        mg->GetXaxis()->SetRangeUser(0., double(iwr));
        // mg->GetYaxis()->SetRangeUser(0,ymax) ; 
        mg->GetXaxis()->SetNdivisions(7, 5, 0, true);
        mg->GetXaxis()->SetLabelSize(0.05);
        dummy->Draw();
        mg->Draw("A");
        dummy->Delete();

        TLatex* mark = new TLatex(0.20, 0.95, "CMS Preliminary");
        // TLatex* mark2 = new TLatex(0.6,0.65, "#int Ldt = 227 pb^{-1} at 7 TeV") ;
        TLatex* mark2;
        switch(mode)
        {
            case 0:
                mark2 = new TLatex(0.70, 0.95, "19.6 fb^{-1} at 8 TeV");
                break;
            case 1:
                mark2 = new TLatex(0.70, 0.95, "19.6 fb^{-1} at 8 TeV");
                break;
        }
            
        mark->SetNDC(kTRUE);
        mark->SetTextSize(0.035);
        mark->SetTextFont(42);
        mark2->SetNDC(kTRUE);
        mark2->SetTextSize(0.035);
        mark2->SetTextFont(42);

        //TLatex* mark  = new TLatex(0.0*mg->GetXaxis()->GetXmax(), 1.03*mg->GetYaxis()->GetXmax(), "CMS Preliminary") ;
        //mark->SetTextSize(0.04);

        //TLatex* mark2 = new TLatex(0.62*mg->GetXaxis()->GetXmax(),1.03*mg->GetYaxis()->GetXmax(), "227 pb^{-1} at 7 TeV") ;
        //mark2->SetTextSize(0.04);

        char leglabel[128];
        sprintf(leglabel, "               #lower[0.1]{M_{W_{R}} = %d GeV}", iwr);

        TLegend *leg = new TLegend(0.50, 0.68, 0.94, 0.94, leglabel, "brNDC");
        leg->SetTextSize(0.030);
        leg->SetTextFont(42);
        leg->SetBorderSize(1);
        leg->SetFillColor(10);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->AddEntry(h_obs, "95% C.L. Observed Limit", "L");
        leg->AddEntry(h_exp, "95% C.L. Expected Limit", "L");
        leg->AddEntry(h_exp1sig, "#pm1#sigma Expected Limit", "F");
        leg->AddEntry(h_exp2sig, "#pm2#sigma Expected Limit", "F");
        leg->AddEntry(h_theoryL, "NNLO Signal cross section", "L");
        //leg->AddEntry(h_theory, "#pm1#sigma PDF+scale unc.", "F");
        leg->Draw();

        mark->Draw();
        mark2->Draw();

        char cn[128], png[128], eps[128];
        switch(mode)
        {
            case 0:
                sprintf(cn, "Mu_Limit_1D_mWR_%d.pdf", iwr);
                sprintf(png, "Mu_Limit_1D_mWR_%d.png", iwr);
                sprintf(eps, "Mu_Limit_1D_mWR_%d.eps", iwr);
                break;
            case 1:
                sprintf(cn, "Elec_Limit_1D_mWR_%d.pdf", iwr);
                sprintf(png, "Elec_Limit_1D_mWR_%d.png", iwr);
                sprintf(eps, "Elec_Limit_1D_mWR_%d.eps", iwr);
                break;
            case 2:
                sprintf(cn, "El-Mu_Limit_1D_mWR_%d.pdf", iwr);
                sprintf(png, "El-Mu_Limit_1D_mWR_%d.png", iwr);
                sprintf(eps, "El-Mu_Limit_1D_mWR_%d.eps", iwr);
                break;
            case 3:
                sprintf(cn, "Mu_2011-2012_Limit_1D_mWR_%d.pdf", iwr);
                sprintf(png, "Mu_2011-2012_Limit_1D_mWR_%d.png", iwr);
                sprintf(eps, "Mu_2011-2012_Limit_1D_mWR_%d.eps", iwr);
                break;
        }
        csUL->Print(cn);
        csUL->Print(png);
        csUL->Print(eps);
    }

    //
    // Now, create a set of points and make the 2D exclusion plot
    //

    std::vector<XYZ> expLowPts;
    std::vector<XYZ> expHighPts;
    std::vector<XYZ> obsLowPts;
    std::vector<XYZ> obsHighPts;
    std::vector<XYZ> obs_right, exp_right;
    for(int iwr = minval; iwr <= maxval; iwr += MWR_STEP)
    {
        int nuExpLo = 0;
        int nuExpHi = iwr;
        int nuObsLo = 0;
        int nuObsHi = iwr;
        int expExclPts = 0;
        int obsExclPts = 0;
        for(int inu = 10; inu < iwr; inu += MNU_STEP)
        {
            std::pair<int, int> masspoint(iwr, inu);
            double sigmaBR = grCS->Interpolate(iwr, inu);
            if(sigmaBR <= 0) continue;

            double expDiff = sigmaBR - rawExpPts[masspoint];
            if(expDiff < 0)
            {
                if(inu <= (iwr / 2) && (inu > nuExpLo)) nuExpLo = inu;
                if(inu > (iwr / 2) && (inu < nuExpHi)) nuExpHi = inu;
            }
            else
            {
                expExclPts++;
            }
            if(iwr > TRANSITION_MWR && expDiff > 0 && inu < iwr - 50)
            {
                std::pair<int, int> nmasspoint(iwr + MWR_STEP, inu);
                double nextSigmaBR = grCS->Interpolate(iwr + MWR_STEP, inu);
                double nextExpDiff = nextSigmaBR - rawExpPts[nmasspoint];
                if(nextExpDiff < 0)
                {
                    double x1 = iwr;
                    double y1 = expDiff;
                    double x2 = iwr + MWR_STEP;
                    double y2 = nextExpDiff;
                    double interpol = -y1 * (x2 - x1) / (y2 - y1) + x1; // notice possible divide-by-zero!
                    XYZ p = {interpol, inu, 0};
                    if(interpol > x1 && interpol < x2)
                        exp_right.push_back(p);
                }
            }

            double obsDiff = sigmaBR - rawObsPts[masspoint];
            if(obsDiff < 0)
            {
                if(inu <= (iwr / 2) && (inu > nuObsLo)) nuObsLo = inu;
                if(inu > (iwr / 2) && (inu < nuObsHi)) nuObsHi = inu;
            }
            else
            {
                obsExclPts++;
            }
            if(iwr > TRANSITION_MWR && obsDiff > 0 && inu < iwr - 50)
            {
                std::pair<int, int> nmasspoint(iwr + MWR_STEP, inu);
                double nextSigmaBR = grCS->Interpolate(iwr + MWR_STEP, inu);
                double nextObsDiff = nextSigmaBR - rawObsPts[nmasspoint];
                if(nextObsDiff < 0)
                {
                    double x1 = iwr;
                    double y1 = obsDiff;
                    double x2 = iwr + MWR_STEP;
                    double y2 = nextObsDiff;
                    double interpol = -y1 * (x2 - x1) / (y2 - y1) + x1; // notice possible divide-by-zero!
                    XYZ p = {interpol, inu, 0};
                    if(interpol > x1 && interpol < x2)
                        obs_right.push_back(p);
                }
            }
        }
        if(iwr <= TRANSITION_MWR)
        {
            if(nuObsLo > 0 && obsExclPts > 0)
            {
                XYZ p = {iwr, nuObsLo, 0.};
                obsLowPts.push_back(p);
            }
            if(nuObsHi < iwr && obsExclPts > 0)
            {
                XYZ p = {iwr, nuObsHi, 0.};
                obsHighPts.push_back(p);
            }
            if(nuExpLo > 0 && expExclPts > 0)
            {
                XYZ p = {iwr, nuExpLo, 0.};
                expLowPts.push_back(p);
            }
            if(nuExpHi < iwr && expExclPts > 0)
            {
                XYZ p = {iwr, nuExpHi, 0.};
                expHighPts.push_back(p);
            }
        }
    }

    std::sort(obs_right.begin(), obs_right.end(), compareByNu2());
    std::sort(exp_right.begin(), exp_right.end(), compareByNu2());

    TGraph* observedPlot = new TGraph();
    TGraph* expectedPlot = new TGraph();

    unsigned int istart = 0;
    if(minval > 700)
    {
        observedPlot->SetPoint(istart, 700, 50);
        expectedPlot->SetPoint(istart, 700, 50);
        istart++;
    }
    if(minval > 800)
    {
        observedPlot->SetPoint(istart, 800, 54);
        expectedPlot->SetPoint(istart, 800, 52);
        istart++;
    }
    if(minval > 900)
    {
        observedPlot->SetPoint(istart, 900, 56);
        expectedPlot->SetPoint(istart, 900, 59);
        istart++;
    }

    /*    int ctr2d = istart;
    for(unsigned int i = 0; i < obsLowPts.size(); i++)
    {
        if(obsLowPts.at(i).x % 100 == 0)
        {
            if(i == 0) observedPlot->SetPoint(ctr2d++, obsLowPts.at(i).x, (obsLowPts.at(i).y + 0.5 * obsLowPts.at(i + 1).y) / 1.5);
            if(i == obsLowPts.size() - 1) observedPlot->SetPoint(ctr2d++, obsLowPts.at(i).x, (obsLowPts.at(i).y + 0.5 * obsLowPts.at(i - 1).y) / 1.5);
            else observedPlot->SetPoint(ctr2d++, obsLowPts.at(i).x, (obsLowPts.at(i).y + 0.5 * obsLowPts.at(i + 1).y + 0.5 * obsLowPts.at(i - 1).y) / 2.0);
        }
    }
    for(unsigned int i = 0; i < obs_right.size(); i++) observedPlot->SetPoint(ctr2d++, obs_right.at(i).x, obs_right.at(i).y);
    for(int i = obsHighPts.size() - 1; i < obsHighPts.size(); i--)
    {
        if(obsHighPts.at(i).x % 100 == 0)
        {
            //std::cout << "HERE: " << i << "\t" << obsHighPts.at(i).x << std::endl;
            if(i == 0) observedPlot->SetPoint(ctr2d++, obsHighPts.at(i).x, (obsHighPts.at(i).y + 0.5 * obsHighPts.at(i + 1).y) / 1.5);
            if(i == obsHighPts.size() - 1) observedPlot->SetPoint(ctr2d++, obsHighPts.at(i).x, (obsHighPts.at(i).y + 0.5 * obsHighPts.at(i - 1).y) / 1.5);
            else observedPlot->SetPoint(ctr2d++, obsHighPts.at(i).x, (obsHighPts.at(i).y + 0.5 * obsHighPts.at(i + 1).y + 0.5 * obsHighPts.at(i - 1).y) / 2.0);
        }
    }
    ctr2d = istart;
    for(unsigned int i = 0; i < expLowPts.size(); i++) expectedPlot->SetPoint(ctr2d++, expLowPts.at(i).x, expLowPts.at(i).y);
    for(unsigned int i = 0; i < exp_right.size(); i++) expectedPlot->SetPoint(ctr2d++, exp_right.at(i).x, exp_right.at(i).y);
    for(unsigned int i = expHighPts.size(); i > 0; i--) expectedPlot->SetPoint(ctr2d++, expHighPts.at(i - 1).x, expHighPts.at(i - 1).y);*/

    int ctr2d = istart;
    for(unsigned int i = 0; i < obsLowPts.size(); i++) observedPlot->SetPoint(ctr2d++, obsLowPts.at(i).x, obsLowPts.at(i).y);
    for(unsigned int i = 0; i < obs_right.size(); i++) observedPlot->SetPoint(ctr2d++, obs_right.at(i).x, obs_right.at(i).y);
    for(unsigned int i = obsHighPts.size(); i > 0; i--) observedPlot->SetPoint(ctr2d++, obsHighPts.at(i - 1).x, obsHighPts.at(i - 1).y);
    ctr2d = istart;
    for(unsigned int i = 0; i < expLowPts.size(); i++) expectedPlot->SetPoint(ctr2d++, expLowPts.at(i).x, expLowPts.at(i).y);
    for(unsigned int i = 0; i < exp_right.size(); i++) expectedPlot->SetPoint(ctr2d++, exp_right.at(i).x, exp_right.at(i).y);
    for(unsigned int i = expHighPts.size(); i > 0; i--) expectedPlot->SetPoint(ctr2d++, expHighPts.at(i - 1).x, expHighPts.at(i - 1).y);

    if(minval > 900)
    {
        observedPlot->SetPoint(observedPlot->GetN(), 900, 819);
        expectedPlot->SetPoint(expectedPlot->GetN(), 900, 810);
    }
    if(minval > 800)
    {
        observedPlot->SetPoint(observedPlot->GetN(), 800, 722);
        expectedPlot->SetPoint(expectedPlot->GetN(), 800, 726);
    }
    if(minval > 700)
    {
        observedPlot->SetPoint(observedPlot->GetN(), 700, 645);
        expectedPlot->SetPoint(expectedPlot->GetN(), 700, 641);
    }

    double xx, yy;
    std::cout << "OBSERVED" << std::endl;
    for(int i = 0; i < observedPlot->GetN(); i++)
    {
        observedPlot->GetPoint(i, xx, yy);
        std::cout << xx << "\t" << yy << std::endl;
    }
    std::cout << "EXPECTED" << std::endl;
    for(int i = 0; i < expectedPlot->GetN(); i++)
    {
        expectedPlot->GetPoint(i, xx, yy);
        if(xx < 1000) std::cout << xx << "\t" << yy << std::endl;
    }

    TCanvas *twoDlimits = new TCanvas("twoDlimits", "twoDlimits", 800, 800);
    twoDlimits->SetRightMargin(0.03);
    twoDlimits->SetLeftMargin(0.18);
    setTDRStyle();
    fixOverlay();

    TH2D* grid = new TH2D("grid", "grid", 100, xmin, xmax, 100, 0, ymax); //1000) ;
    grid->GetXaxis()->SetTitle("M_{W_{R}} [GeV]");
    grid->GetYaxis()->SetTitleOffset(1.45);
    grid->GetXaxis()->SetNdivisions(6, 5, 0);

    switch(mode)
    {
        case 0:
            grid->GetYaxis()->SetTitle("M_{N_{#mu}} [GeV]");
            expectedPlot->SetLineColor(kBlue);
            break;
        case 1:
            grid->GetYaxis()->SetTitle("M_{N_{e}} [GeV]");
            expectedPlot->SetLineColor(kRed);
            observedPlot->SetLineColor(kRed);
            break;
        case 2:
            grid->GetYaxis()->SetTitle("M_{N_{e,#mu,#tau}} [GeV]");
            expectedPlot->SetLineColor(kMagenta + 2);
            break;
        case 3:
            grid->GetYaxis()->SetTitle("M_{N_{#mu}} [GeV]");
            expectedPlot->SetLineColor(kBlue);
            break;
    }
    grid->Draw();

    
    observedPlot->SetLineColor(kBlack);
    observedPlot->SetFillStyle(3005);
    observedPlot->SetLineWidth(5);
    observedPlot->Draw("F SAME");
    observedPlot->Draw("L SAME");

    expectedPlot->SetLineStyle(2);
    expectedPlot->SetLineWidth(5);
    expectedPlot->Draw("L SAME");

    TLegend *leg = new TLegend(.81, .86, .96, .94);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->AddEntry(observedPlot, "Observed", "L");
    leg->AddEntry(expectedPlot, "Expected", "L");

    // ============== plot Tevatron limit region:

    //float x[] = {750, 890, 890, 750};
    //float y[] = {0, 0, 890, 750};
    //TGraph* Tevatron = new TGraph(4, x, y);
    //Tevatron->SetFillColor(kGray + 1);
    //Tevatron->SetLineColor(kGray);

    // ============== plot region M_nuR > M_WR:
    float x2[] = {xmin, ymax, xmin};
    float y2[] = {ymax, ymax, xmin};
    TGraph* wrnu = new TGraph(3, x2, y2);
    wrnu->SetLineWidth(3);
    wrnu->SetLineColor(kYellow - 4);
    wrnu->SetFillColor(kYellow);

    // ============= plot labels:
    TLatex text;
    text.SetTextSize(0.04);
    text.SetTextFont(42);

    TLatex text2;
    text2.SetTextSize(0.035);
    text2.SetTextAngle(90);
    text2.SetTextFont(42);

    TLatex* mark = new TLatex(1300, 1.015 * ymax, "CMS Preliminary");
    mark->SetNDC();
    mark->SetTextSize(0.04*1.1);
    mark->SetTextFont(42);


    //TLatex* mark2 = new TLatex(1300, 0.93 * ymax, "4.7 fb^{-1} at 7 TeV");
    //mark2->SetTextSize(0.04);
    //mark2->SetTextFont(42);


    //Tevatron->Draw("F SAME");
    wrnu->Draw("F SAME");
    switch(mode)
    {
        case 0:
            text.DrawLatex(xmin + 100, .92 * ymax, "M_{N_{#mu}} > M_{W_{R}}");
            mark->DrawLatex(0.220,0.9575, "CMS Preliminary    #sqrt{s} = 8 TeV    19.6 fb^{-1}");
            //mark->DrawLatex(0.71, 0.96, "12.1 fb^{-1} at 8 TeV");
            break;
        case 1:
            text.DrawLatex(xmin + 100, .92 * ymax, "M_{N_{e}} > M_{W_{R}}");
            mark->DrawLatex(0.220,0.9575, "CMS Preliminary    #sqrt{s} = 8 TeV    19.6 fb^{-1}");
            //mark->DrawLatex(0.71, 0.96, "12.3 fb^{-1} at 8 TeV");
            break;
        case 2:
            text.DrawLatex(xmin + 100, .92 * ymax, "M_{N_{e,#mu,#tau}} > M_{W_{R}}");
            mark->DrawLatex(0.71, 0.96, "3.6 fb^{-1} at 8 TeV");
            break;
        case 3:
            text.DrawLatex(xmin + 100, .92 * ymax, "M_{N_{#mu}} > M_{W_{R}}");
            mark->DrawLatex(0.51, 0.96, "5.0 fb^{-1} (7 TeV) + 3.6 fb^{-1} (8 TeV)");
            break;

    }
    //text2.DrawLatex(880, 20, "Excluded by Tevatron");
    //mark->DrawLatex(0.18, 0.96, "CMS Preliminary");

    fixOverlay();

    leg->Draw("same");

    switch(mode)
    {
        case 0:
            twoDlimits->Print("Mu_Lim2D.pdf");
            twoDlimits->Print("Mu_Lim2D.png");
            twoDlimits->Print("Mu_Lim2D.eps");
            break;
        case 1:
            twoDlimits->Print("Elec_Lim2D.pdf");
            twoDlimits->Print("Elec_Lim2D.png");
            twoDlimits->Print("Elec_Lim2D.eps");
            break;
        case 2:
            twoDlimits->Print("El-Mu_Lim2D.pdf");
            twoDlimits->Print("El-Mu_Lim2D.png");
            twoDlimits->Print("El-Mu_Lim2D.eps");
            break;
        case 3:
            twoDlimits->Print("Mu_2011-2012_Lim2D.pdf");
            twoDlimits->Print("Mu_2011-2012_Lim2D.png");
            twoDlimits->Print("Mu_2011-2012_Lim2D.eps");
            break;
    }

    //delete[] db;

}
