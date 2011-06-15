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
// #include "../acceptance/db_jun6.C"

#include <iostream>
#include <fstream>

struct XYZ
{
    int x, y;
    float z;

    bool operator ==(const XYZ & v) const
    {
        return(x == v.x) && (y == v.y);
    }
};

static const double lumi2010 = 36.12;
static const double lumi2011 = 204.16;

static const AcceptanceDB db;

class compareByNu
{
public:

    bool operator() (const XYZ& a, const XYZ& b)
    {
        return a.y > b.y;
    }
};

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

double getWgtAcceptance(int imwr, int imnu)
{

    double acc2010 = db.getBestEstimate(imwr, imnu, 2010);
    double acc2011 = db.getBestEstimate(imwr, imnu, 2011);
    double wgtAcc = (lumi2010 * acc2010 + lumi2011 * acc2011) / (lumi2010 + lumi2011);

    return wgtAcc;
}

double getEfficiency(int imwr, std::vector<XYZ> v)
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
    if(nearestLo != nearestHi)
    { // Not a reference point
        double wgtLoAcc = getWgtAcceptance(v.at(nearestLo).x, v.at(nearestLo).y);
        double wgtHiAcc = getWgtAcceptance(v.at(nearestHi).x, v.at(nearestHi).y);

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
            double wgtAcc = getWgtAcceptance(v.at(nearestLo).x, v.at(nearestLo).y);
            // std::cout << "wgtAcc is " << wgtAcc << std::endl ; 
            eff = v.at(nearestLo).z * wgtAcc;
            // std::cout << "eff: " << eff << std::endl ; 
        }
    }
    return eff;
}

void plotLimits(int maxval = -1, bool doMuon = true, std::string limitFile = "limitSummary.txt", std::string csFile = "cs.txt")
{

    if(maxval < 0) maxval = 1700;


    // Inizialization
    FILE *Ptr;
    char *str1 = (char*)alloca(10000);
    ifstream inLimit(limitFile.c_str());
    ifstream inCS(csFile.c_str());

    std::string outputFileName = "obsexpSummary.txt";
    Ptr = fopen(outputFileName.c_str(), "a");

    // static const int numPts = 133 ; 
    static const int numPts = 400;

    Double_t mWR[numPts];
    Double_t mNu[numPts];
    Double_t mWRcs[numPts];
    Double_t mNucs[numPts];
    Double_t theCS[numPts];
    Double_t obs[numPts];
    Double_t exp[numPts];
    Double_t exp68lo[numPts];
    Double_t exp68hi[numPts];
    Double_t exp95lo[numPts];
    Double_t exp95hi[numPts];

    double csRelError = 0.10;

    int ctr = 0;
    while(!inCS.eof())
    {
        double cs, sf;
        inCS >> mWRcs[ctr];
        inCS >> mNucs[ctr];
        inCS >> cs;
        inCS >> sf;
        if(!inCS.good()) break;
        theCS[ctr] = cs * sf;
        ctr++;
        if(ctr == numPts) break;
    }

    // Build up cross-section as a function of W_R and N mass
    std::vector<XYZ> ps;
    for(int i = 0; i < ctr; i++)
    {
        const XYZ p = {mWRcs[i], mNucs[i], theCS[i]};
        if(find(ps.begin(), ps.end(), p) == ps.end()) ps.push_back(p);
    }
    TGraph2D* grCS = new TGraph2D(ps.size());
    grCS->SetNameTitle("xs", "xs;WR;nuR;#sigma");
    for(unsigned int i = 0; i < ps.size(); i++) grCS->SetPoint(i, ps[i].x, ps[i].y, ps[i].z);

    ctr = 0;
    while(!inLimit.eof())
    {
        inLimit >> mWR[ctr];
        inLimit >> mNu[ctr];
        inLimit >> obs[ctr];
        inLimit >> exp[ctr];
        inLimit >> exp68lo[ctr];
        inLimit >> exp68hi[ctr];
        inLimit >> exp95lo[ctr];
        inLimit >> exp95hi[ctr];

        if(!inLimit.good()) break;
        ctr++;
        if(ctr >= numPts) break;
    }

    // 
    // Build up observed/expected limit as a function of W_R and N mass
    // 
    std::vector<XYZ> obsLimit;
    std::vector<XYZ> expLimit;
    std::vector<XYZ> exp68loLimit;
    std::vector<XYZ> exp68hiLimit;
    std::vector<XYZ> exp95loLimit;
    std::vector<XYZ> exp95hiLimit;
    for(int i = 0; i < ctr; i++)
    {
        const XYZ p = {mWR[i], mNu[i], obs[i]};
        const XYZ q = {mWR[i], mNu[i], exp[i]};
        const XYZ h68 = {mWR[i], mNu[i], exp68hi[i]};
        const XYZ l68 = {mWR[i], mNu[i], exp68lo[i]};
        const XYZ h95 = {mWR[i], mNu[i], exp95hi[i]};
        const XYZ l95 = {mWR[i], mNu[i], exp95lo[i]};
        if(find(obsLimit.begin(), obsLimit.end(), p) == obsLimit.end()) obsLimit.push_back(p);
        if(find(expLimit.begin(), expLimit.end(), q) == expLimit.end()) expLimit.push_back(q);
        if(find(exp68hiLimit.begin(), exp68hiLimit.end(), h68) == exp68hiLimit.end()) exp68hiLimit.push_back(h68);
        if(find(exp68loLimit.begin(), exp68loLimit.end(), l68) == exp68loLimit.end()) exp68loLimit.push_back(l68);
        if(find(exp95hiLimit.begin(), exp95hiLimit.end(), h95) == exp95hiLimit.end()) exp95hiLimit.push_back(h95);
        if(find(exp95loLimit.begin(), exp95loLimit.end(), l95) == exp95loLimit.end()) exp95loLimit.push_back(l95);

        //     const XYZ cs = {mWR[i], mNu[i], theCS[i]} ;
        //     if ( find(sigmaBR.begin(),sigmaBR.end(),cs) == sigmaBR.end() ) sigmaBR.push_back( cs ) ;
    }

    for(unsigned int i = 0; i < obsLimit.size(); i++)
    {
        if(obsLimit.at(i).x == 700)
            std::cout << " mwr: " << obsLimit.at(i).x
                << " mnu: " << obsLimit.at(i).y
                << " val: " << obsLimit.at(i).z
                << std::endl;
    }

    std::vector<XYZ> extraObsPts;
    std::vector<XYZ> extraExpPts;
    for(int imwr = 700; imwr <= 1700; imwr += 100)
    {

        if(imwr % 100 == 0) std::cout << "Computing for W_R mass : " << imwr << std::endl;

        double eff_wr = getEfficiency(imwr, obsLimit);
        for(int imnu = 10; imnu < imwr; imnu += 1)
        {
            double acc = getWgtAcceptance(imwr, imnu);
            if(acc < 0) continue;
            std::cout << "For " << imwr << ", " << imnu << " get acceptance " << acc << std::endl ;
            XYZ p = {imwr, imnu, eff_wr / acc};
            // if ( find(obsLimit.begin(),obsLimit.end(),p) == obsLimit.end() )
            extraObsPts.push_back(p);
            //           else
            //               std::cout << "Check: " << p.z << " vs. "
            //                         << (*find(obsLimit.begin(),obsLimit.end(),p)).z
            //                         << std::endl ;
        }

        eff_wr = getEfficiency(imwr, expLimit);
        for(int imnu = 10; imnu < imwr; imnu += 1)
        {
            double acc = getWgtAcceptance(imwr, imnu);
            if(acc < 0) continue;
            XYZ p = {imwr, imnu, eff_wr / acc};
            // if ( find(expLimit.begin(),expLimit.end(),p) == expLimit.end() )
            extraExpPts.push_back(p);
            //           else
            //               std::cout << "Check (Exp): " << p.z << " vs. "
            //                         << (*find(expLimit.begin(),expLimit.end(),p)).z
            //                         << std::endl ;
        }
    }




    //
    // 1D plots
    //
    std::vector<XYZ> exp68loLimitFull;
    std::vector<XYZ> exp68hiLimitFull;
    std::vector<XYZ> exp95loLimitFull;
    std::vector<XYZ> exp95hiLimitFull;
    for(int iwr = 700; iwr <= 1700; iwr += 100)
    {

        TGraph* h_exp = new TGraph();
        TGraph* h_obs = new TGraph();
        TGraph* h_theoryL = new TGraph();
        TGraphAsymmErrors* h_theory = new TGraphAsymmErrors();
        TGraphAsymmErrors* h_exp1sig = new TGraphAsymmErrors();
        TGraphAsymmErrors* h_exp2sig = new TGraphAsymmErrors();

        double eff_68lo = getEfficiency(iwr, exp68loLimit);
        double eff_68hi = getEfficiency(iwr, exp68hiLimit);
        double eff_95lo = getEfficiency(iwr, exp95loLimit);
        double eff_95hi = getEfficiency(iwr, exp95hiLimit);
        for(int inu = 10; inu < iwr; inu += 1)
        {
            double acc = getWgtAcceptance(iwr, inu);
            if(acc < 0) continue;
            XYZ pl68 = {iwr, inu, eff_68lo / acc};
            XYZ ph68 = {iwr, inu, eff_68hi / acc};
            XYZ pl95 = {iwr, inu, eff_95lo / acc};
            XYZ ph95 = {iwr, inu, eff_95hi / acc};
            // if ( find(exp68loLimit.begin(),exp68loLimit.end(),pl68) == exp68loLimit.end() )
            exp68loLimitFull.push_back(pl68);
            // if ( find(exp68hiLimit.begin(),exp68hiLimit.end(),ph68) == exp68hiLimit.end() )
            exp68hiLimitFull.push_back(ph68);
            // if ( find(exp95loLimit.begin(),exp95loLimit.end(),pl95) == exp95loLimit.end() )
            exp95loLimitFull.push_back(pl95);
            // if ( find(exp95hiLimit.begin(),exp95hiLimit.end(),ph95) == exp95hiLimit.end() )
            exp95hiLimitFull.push_back(ph95);

            double sigmaBR = grCS->Interpolate(iwr, inu);
            // if ( sigmaBR < 0.000001 ) continue ; 
            h_theoryL->SetPoint(inu - 10, inu, sigmaBR);
            h_theory->SetPoint(inu - 10, inu, sigmaBR);
            h_theory->SetPointError(inu - 10, 0., 0., sigmaBR * 0.1, sigmaBR * 0.1);
        }


        // Find points, fill values
        // std::vector<XYZ> expectedPts = getSortedVector(iwr,expLimit,extraExpPts) ; 
        // std::vector<XYZ> observedPts = getSortedVector(iwr,obsLimit,extraObsPts) ; 
        std::vector<XYZ> expectedPts = getSortedVector(iwr, extraExpPts);
        std::vector<XYZ> observedPts = getSortedVector(iwr, extraObsPts);

        std::vector<XYZ> expected68loPts = getSortedVector(iwr, exp68loLimitFull);
        std::vector<XYZ> expected68hiPts = getSortedVector(iwr, exp68hiLimitFull);
        std::vector<XYZ> expected95loPts = getSortedVector(iwr, exp95loLimitFull);
        std::vector<XYZ> expected95hiPts = getSortedVector(iwr, exp95hiLimitFull);
        
        std::cout << "SIZES: \t" << expectedPts.size() << "\t" << observedPts.size() << "\t" << expected68loPts.size()
             << "\t" << expected68hiPts.size()  << "\t" << expected95loPts.size()  << "\t" << expected95hiPts.size() << std::endl;

        std::cout << h_theory->GetN() << std::endl;
        for(int i = 0; i < h_theory->GetN(); i++)
        {
            double x, y;
            h_theory->GetPoint(i, x, y);
            if(iwr == 1000) std::cout << x << ": " << y << std::endl ;
            if(x < 10 || x > iwr || y < 1e-5 || y > 1e2)
            {
                h_theoryL->RemovePoint(i);
                h_theory->RemovePoint(i);
                i--;
            }
        }
        std::cout << h_theory->GetN() << std::endl;

        for(unsigned int i = 0; i < expectedPts.size(); i++)
        {
            h_exp->SetPoint(i, expectedPts.at(i).y, expectedPts.at(i).z);
            h_exp1sig->SetPoint(i, expectedPts.at(i).y, expectedPts.at(i).z);
            h_exp2sig->SetPoint(i, expectedPts.at(i).y, expectedPts.at(i).z);
        }
        for(unsigned int i = 0; i < observedPts.size(); i++)
            h_obs->SetPoint(i, observedPts.at(i).y, observedPts.at(i).z);

        for(unsigned int i = 0; i < expected68loPts.size(); i++)
        {
            double d68lo = expectedPts.at(i).z - expected68loPts.at(i).z;
            double d68hi = expected68hiPts.at(i).z - expectedPts.at(i).z;
            double d95lo = expectedPts.at(i).z - expected95loPts.at(i).z;
            double d95hi = expected95hiPts.at(i).z - expectedPts.at(i).z;
            h_exp1sig->SetPointError(i, 0., 0., d68lo, d68hi);
            h_exp2sig->SetPointError(i, 0., 0., d95lo, d95hi);
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
        std::string cLabel(str1);

        TCanvas* csExp = new TCanvas("csExp");
        setTDRStyle();
        fixOverlay();

        TCanvas* csUL = new TCanvas(cLabel.c_str());
        csUL->SetLeftMargin(0.19);
        csUL->SetRightMargin(0.06);
        csUL->SetTopMargin(0.06);
        csUL->SetLogy(1);

        TMultiGraph* mg = new TMultiGraph();
        mg->SetName("");
        mg->SetTitle(";M_{N_{#mu}} [GeV];#sigma(pp #rightarrow W_{R})#times Br(W_{R} #rightarrow #mu#mujj) [pb]");
        if(!doMuon) mg->SetTitle(";M_{N_{e}} [GeV];#sigma(pp #rightarrow W_{R})#times Br(W_{R} #rightarrow eejj) [pb]");

        mg->Add(h_exp2sig, "E3");
        mg->Add(h_exp1sig, "E3");
        mg->Add(h_theory, "E3");
        mg->Add(h_theoryL, "L");
        mg->Add(h_exp, "L");
        mg->Add(h_obs, "L");

        TH1D* dummy = new TH1D("dummy", "dummy", 100, 0, double(iwr));
        dummy->GetXaxis()->SetTitle("M_{N_{#mu}} [GeV]");
        dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow #mu#mujj) [pb]");
        if(!doMuon)
        {
            dummy->GetXaxis()->SetTitle("M_{N_{e}} [GeV]");
            dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow W_{R}) #times Br(W_{R} #rightarrow eejj) [pb]");
        }
        dummy->GetYaxis()->SetTitleOffset(1.45);
        // dummy->SetMinimum(minyval) ; 
        // dummy->SetMaximum(maxyval) ; 

        mg->Draw("AH");
        mg->GetXaxis()->SetRangeUser(0., double(iwr));
        // mg->GetYaxis()->SetRangeUser(0,ymax) ; 
        dummy->Draw();
        mg->Draw("A");
        dummy->Delete();

        TLegend *leg = new TLegend(0.55, 0.72, 0.94, 0.94, NULL, "brNDC");
        leg->SetTextSize(0.025);
        leg->AddEntry(h_obs, "95% C.L. Observed Limit", "L");
        leg->AddEntry(h_exp, "95% C.L. Expected Limit", "L");
        leg->AddEntry(h_exp1sig, "#pm1#sigma Expected Limit", "F");
        leg->AddEntry(h_exp2sig, "#pm2#sigma Expected Limit", "F");
        leg->AddEntry(h_theoryL, "NNLO Signal cross section", "L");
        leg->AddEntry(h_theory, "#pm1#sigma PDF+scale unc.", "F");
        leg->Draw();
        
        char cn[128];
        sprintf(cn, "Limit_1D_mWR_%d.pdf", iwr);
        csUL->Print(cn);

    }
    

    // Now, create a set of points and make the 2D exclusion plot
    std::vector<XYZ> expLowPts;
    std::vector<XYZ> expHighPts;
    std::vector<XYZ> obsLowPts;
    std::vector<XYZ> obsHighPts;
    std::vector<XYZ> exp_justExcluded, exp_justNotExcluded;
    std::vector<XYZ> obs_justExcluded, obs_justNotExcluded;
    for(int iwr = 700; iwr <= 1700; iwr += 100)
    {
        std::vector<XYZ> expectedPts = getSortedVector(iwr, extraExpPts);
        std::vector<XYZ> observedPts = getSortedVector(iwr, extraObsPts);

        int nuExpLo = 0;
        int nuExpHi = iwr;
        int nuObsLo = 0;
        int nuObsHi = iwr;
        int expExclPts = 0;
        int obsExclPts = 0;
        for(int inu = 10; inu < iwr; inu += 1)
        {
            double sigmaBR = grCS->Interpolate(iwr, inu);
            if(sigmaBR <= 0) continue;
            for(unsigned int i = 0; i < expectedPts.size(); i++)
            {
                if(expectedPts.at(i).y == inu)
                {
                    double expDiff = sigmaBR - expectedPts.at(i).z;
                    //if(iwr == 1700) std::cout << "Exp: " << iwr << ", " << inu << ": "
                    //       << sigmaBR << " - " << expectedPts.at(i).z << " = "
                    //       << expDiff << std::endl;
                    if(expDiff < 0)
                    {
                        if(inu <= (iwr / 2) && (inu > nuExpLo)) nuExpLo = inu;
                        if(inu > (iwr / 2) && (inu < nuExpHi)) nuExpHi = inu;
                        if(iwr == 1700) std::cout << "Exp vals: " << nuExpLo << ", " << nuExpHi << std::endl;
                    }
                    else
                    {
                        expExclPts++;
                    }
                    if (iwr==1500 && expDiff>0) {
                    	XYZ p = {iwr, inu, expDiff};
                    	exp_justExcluded.push_back(p);
                    } 
                    if (iwr==1600) {
                    	XYZ p = {iwr, inu, expDiff};
                    	exp_justNotExcluded.push_back(p);
                    } 

                }
            }
            for(unsigned int i = 0; i < observedPts.size(); i++)
            {
                if(observedPts.at(i).y == inu)
                {
                    double obsDiff = sigmaBR - observedPts.at(i).z;
                    if(iwr == 1700) std::cout << "Obs: " << iwr << ", " << inu << ": "
                            << sigmaBR << " - " << observedPts.at(i).z << " = "
                            << obsDiff << std::endl;
                    if(obsDiff < 0)
                    {
                        if(inu <= (iwr / 2) && (inu > nuObsLo)) nuObsLo = inu;
                        if(inu > (iwr / 2) && (inu < nuObsHi)) nuObsHi = inu;
                        if(iwr == 1700) std::cout << "Obs vals: " << nuObsLo << ", " << nuObsHi << std::endl;
                    }
                    else
                    {
                        obsExclPts++;
                    }
                    if (iwr==1500 && obsDiff>0) {
                    	XYZ p = {iwr, inu, obsDiff};
                    	obs_justExcluded.push_back(p);
                    } 
                    if (iwr==1600) {
                    	XYZ p = {iwr, inu, obsDiff};
                    	obs_justNotExcluded.push_back(p);
                    } 
                    
                }
            }
        }
        if(iwr < 1500)
        {
        	if(nuObsLo > 0 && obsExclPts > 0)
	        {
	            XYZ p = {iwr, nuObsLo, 0.};
	            obsLowPts.push_back(p);
	        }
	        else if(iwr == 700)
	        {
	        	XYZ p = {iwr, 50.0, 0.};
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
	        else if(iwr == 700)
	        {
	        	XYZ p = {iwr, 50.0, 0.};
	            expLowPts.push_back(p);
	        }
	        if(nuExpHi < iwr && expExclPts > 0)
	        {
	            XYZ p = {iwr, nuExpHi, 0.};
	            expHighPts.push_back(p);
	        }
        }

    }
    
    std::vector<XYZ> obs_right, exp_right;
    for (unsigned int i=0; i<obs_justExcluded.size(); i+=5)    
    	for (unsigned int j=0; j<obs_justNotExcluded.size(); j++)
    		if (obs_justExcluded.at(i).y == obs_justNotExcluded.at(j).y) {
    			double x1=obs_justExcluded.at(i).x;
    			double y1=obs_justExcluded.at(i).z;
    			double x2=obs_justNotExcluded.at(j).x;
    			double y2=obs_justNotExcluded.at(j).z;
    			double interpol=-y1*(x2-x1)/(y2-y1)+x1; // notice possible divide-by-zero!
    			XYZ p = { interpol, obs_justExcluded.at(i).y, 0 };
    			if (interpol>x1 && interpol<x2)
     			  obs_right.push_back(p);	
    			}

    for (unsigned int i=0; i<exp_justExcluded.size(); i+=5)    
    	for (unsigned int j=0; j<exp_justNotExcluded.size(); j++)
    		if (exp_justExcluded.at(i).y == exp_justNotExcluded.at(j).y) {
    			double x1=exp_justExcluded.at(i).x;
    			double y1=exp_justExcluded.at(i).z;
    			double x2=exp_justNotExcluded.at(j).x;
    			double y2=exp_justNotExcluded.at(j).z;
    			double interpol=-y1*(x2-x1)/(y2-y1)+x1; // notice possible divide-by-zero!
    			XYZ p = { interpol, exp_justExcluded.at(i).y, 0 };
    			exp_right.push_back(p);
    			}
    			

    

    for(unsigned int i = 0; i < expLowPts.size(); i++)
        std::cout << "exp lo: " << expLowPts.at(i).x << " " << expLowPts.at(i).y << std::endl;
    for(unsigned int i = 0; i < expHighPts.size(); i++)
        std::cout << "exp hi: " << expHighPts.at(i).x << " " << expHighPts.at(i).y << std::endl;

    for(unsigned int i = 0; i < obsLowPts.size(); i++)
        std::cout << "obs lo: " << obsLowPts.at(i).x << " " << obsLowPts.at(i).y << std::endl;
    for(unsigned int i = 0; i < obsHighPts.size(); i++)
        std::cout << "obs hi: " << obsHighPts.at(i).x << " " << obsHighPts.at(i).y << std::endl;

    TGraph* observedPlot = new TGraph();
    TGraph* expectedPlot = new TGraph();

    int ctr2d = 0;
    for(unsigned int i = 0; i < obsLowPts.size(); i++) observedPlot->SetPoint(ctr2d++, obsLowPts.at(i).x, obsLowPts.at(i).y);
   for(unsigned int i = 0; i < obs_right.size(); i++) observedPlot->SetPoint(ctr2d++, obs_right.at(i).x, obs_right.at(i).y);   
    for(unsigned int i = obsHighPts.size(); i > 0; i--) observedPlot->SetPoint(ctr2d++, obsHighPts.at(i - 1).x, obsHighPts.at(i - 1).y);
    ctr2d = 0;
    for(unsigned int i = 0; i < expLowPts.size(); i++) expectedPlot->SetPoint(ctr2d++, expLowPts.at(i).x, expLowPts.at(i).y);
    for(unsigned int i = 0; i < exp_right.size(); i++) expectedPlot->SetPoint(ctr2d++, exp_right.at(i).x, exp_right.at(i).y);   
    for(unsigned int i = expHighPts.size(); i > 0; i--) expectedPlot->SetPoint(ctr2d++, expHighPts.at(i - 1).x, expHighPts.at(i - 1).y);

    TCanvas *twoDlimits = new TCanvas("twoDlimits", "twoDlimits", 800, 800);
    twoDlimits->SetRightMargin(0.03);
    setTDRStyle();
    fixOverlay();

    TH2D* grid = new TH2D("grid", "grid", 100, 650, 1700, 100, 0, 1200); //1000) ; 
    grid->GetXaxis()->SetTitle("M_{W_{R}} [GeV]");
    grid->GetYaxis()->SetTitle("M_{N_{#mu}} [GeV]");
    if(!doMuon) grid->GetYaxis()->SetTitle("M_{N_{e}} [GeV]");
    grid->Draw();

    observedPlot->SetLineColor(kBlue);
    observedPlot->SetFillColor(kBlue - 4);
    observedPlot->SetFillStyle(3005);
    observedPlot->SetLineWidth(5);
    observedPlot->Draw("F SAME");
    observedPlot->Draw("L SAME");

    expectedPlot->SetLineColor(kBlack);
    expectedPlot->SetLineStyle(2);
    expectedPlot->SetLineWidth(5);
    expectedPlot->Draw("L SAME");

    TLegend *leg = new TLegend(.8, .8, .95, .91);
    leg->SetBorderSize(1);
    leg->SetFillColor(10);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->AddEntry(observedPlot, "Observed", "L");
    leg->AddEntry(expectedPlot, "Expected", "L");
	
	// ============== plot Tevatron limit region:

    float x[] = {650, 740, 740, 650};
    float y[] = {0, 0, 740, 650};
    TGraph* Tevatron = new TGraph(4, x, y);
    Tevatron->SetFillColor(kGray+1) ;
    Tevatron->SetLineColor(kGray) ;

    // ============== plot region M_nuR > M_WR:
    float ymax = 1200;
    float x2[] = {650, ymax, 650} ;
    float y2[] = {ymax, ymax, 650} ;
    TGraph* wrnu = new TGraph(3, x2, y2);
    wrnu->SetLineWidth(3);
    wrnu->SetLineColor(kYellow-4);
    wrnu->SetFillColor(kYellow);

    // ============= plot labels:

    TLatex text;
    text.SetTextSize(0.04);

    TLatex text2;
    text2.SetTextSize(0.035);
    text2.SetTextAngle(90);

    leg->Draw();
    Tevatron->Draw("F SAME") ;
    wrnu->Draw("F SAME") ;
    text.DrawLatex(720, .94*ymax, "M_{N_{#mu}} > M_{W_{R}}");
    text2.DrawLatex(730, 20, "Excluded by Tevatron");
    
    twoDlimits->Print("lim2D.pdf");

}


