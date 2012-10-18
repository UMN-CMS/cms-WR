#ifndef HEAVY_NU_TOP_HISTSET_INCLUDED
#define HEAVY_NU_TOP_HISTSET_INCLUDED 1

#include "HeavyNu/AnalysisModules/src/HeavyNuHistSet.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class HeavyNuTopHist : public HeavyNuHistSet
{
public:
    HeavyNuTopHist(TFileDirectory *td, const std::string& post, int histsToBook = 1);
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&);
    // fill all histos of the set with the two lepton candidates
    void fill(HeavyNuEvent& hne);

private:
    TH1 *dptMu1gen, *dptMu2gen ;
    TH1 *dRMu1gen, *dRMu2gen ;
    TH1 *qualMu1, *qualMu2 ;
    
    TH1 *muonpt, *muoneta, *muonphi, *elecpt, *eleceta, *elecphi;

    TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
    TH1 *e1trackIso, *e1hcalIso, *e1ecalIso, *e1caloIso, *e1dB;

    TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
    TH1 *e1trackRelIso, *e1hcalRelIso, *e1ecalRelIso, *e1caloRelIso;

    TH1 *mMuEBarrel, *mMuEBB;
    TH1 *mWRBarrel, *mWRBB ;
    TH2 *mWRvsminLPt, *mWRvsNPV;
};

#endif