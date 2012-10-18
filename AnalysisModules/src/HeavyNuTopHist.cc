#include "HeavyNu/AnalysisModules/src/HeavyNuTopHist.h"

static const int muonQualityFlags = 4;
static const std::string muonQuality[4] = {"All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"};

HeavyNuTopHist::HeavyNuTopHist(TFileDirectory *td, const std::string& post, int histsToBook)
{
    if(histsToBook & 1)
    {
        book(td, post);
    }
}

//book histogram set w/ common suffix inside the provided TFileDirectory

void HeavyNuTopHist::book(TFileDirectory * td, const std::string& post)
{
    HeavyNuHistSet::book(td, post);

    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    t = "Quality (#mu_{1}) " + post;
    qualMu1 = td->make<TH1D > ("qualMu1", t.c_str(), muonQualityFlags, 0, muonQualityFlags);
    for (int i = 0; i < muonQualityFlags; i++)
    {
        qualMu1->GetXaxis()->SetBinLabel(i + 1, muonQuality[i].c_str()) ;
    }

    // isolation

    t =   "trackIso(#mu_{1}) " + post;
    mu1trackIso = td->make<TH1F > ("mu1trackIso",   t.c_str(), 40, 0., 200.);
    t =    "hcalIso(#mu_{1}) " + post;
    mu1hcalIso = td->make<TH1F > ("mu1hcalIso",    t.c_str(), 40, 0., 200.);
    t =    "ecalIso(#mu_{1}) " + post;
    mu1ecalIso = td->make<TH1F > ("mu1ecalIso",    t.c_str(), 40, 0., 200.);
    t =    "caloIso(#mu_{1}) " + post;
    mu1caloIso = td->make<TH1F > ("mu1caloIso",    t.c_str(), 40, 0., 200.);
    t =        "Dxy(#mu_{1}) " + post;
    mu1dB = td->make<TH1F > ("mu1dB",         t.c_str(), 50, 0., 1.);

    t =   "trackIso(#e_{1}) " + post;
    e1trackIso = td->make<TH1F > ("mu2trackIso",   t.c_str(), 40, 0., 200.);
    t =    "hcalIso(#e_{1}) " + post;
    e1hcalIso = td->make<TH1F > ("mu2hcalIso",    t.c_str(), 40, 0., 200.);
    t =    "ecalIso(#e_{1}) " + post;
    e1ecalIso = td->make<TH1F > ("mu2ecalIso",    t.c_str(), 40, 0., 200.);
    t =    "caloIso(#e_{1}) " + post;
    e1caloIso = td->make<TH1F > ("mu2caloIso",    t.c_str(), 40, 0., 200.);
    t =        "Dxy(#e_{1}) " + post;
    e1dB = td->make<TH1F > ("mu2dB",         t.c_str(), 50, 0., 1.);

    t = "trackRelIso(#mu_{1}) " + post;
    mu1trackRelIso = td->make<TH1F > ("mu1trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#mu_{1}) " + post;
    mu1hcalRelIso = td->make<TH1F > ("mu1hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#mu_{1}) " + post;
    mu1ecalRelIso = td->make<TH1F > ("mu1ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#mu_{1}) " + post;
    mu1caloRelIso = td->make<TH1F > ("mu1caloRelIso", t.c_str(), 50, 0., 5.);

    t = "trackRelIso(#e_{1}) " + post;
    e1trackRelIso = td->make<TH1F > ("mu2trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#e_{1}) " + post;
    e1hcalRelIso = td->make<TH1F > ("mu2hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#e_{1}) " + post;
    e1ecalRelIso = td->make<TH1F > ("mu2ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#e_{1}) " + post;
    e1caloRelIso = td->make<TH1F > ("mu2caloRelIso", t.c_str(), 50, 0., 5.);

    // ----------  Composite histograms  ----------

    t = "M(W_{R}) e in Barrel"         + post;
    mWRBarrel = td->make<TH1F > ("mWReBarrel", t.c_str(), 70, 0, 2800);
    t = "M(W_{R}) e and #mu in Barrel" + post;
    mWRBB = td->make<TH1F > ("mWRemuBarrel", t.c_str(), 70, 0, 2800);

    t = "M(#mu #mu) e in Barrel"       + post;
    mMuEBarrel = td->make<TH1F > ("mMuMueBarrel", t.c_str(), 50, 0, 2000);
    t = "M(#mu #mu) e and #mu in Barrel" + post;
    mMuEBB = td->make<TH1F > ("mMuMuemuBarrel", t.c_str(), 50, 0, 2000);

    double edges[] = {0, 20, 30, 40, 50, 60, 1500};
    t = "M(W_{R}) vs. min lepton p_{t} " + post;
    mWRvsminLPt = td->make<TH2F > ("mWRvsminLPt",  t.c_str(), 70, 0, 2800, 6, edges);
    t = "M(W_{R}) vs. N primary vertices " + post;
    mWRvsNPV = td->make<TH2F > ("mWRvsNPV",  t.c_str(), 70, 0, 2800, 30, 0, 30);
    
    // -------- particular lepton flavor histograms ---------
    t = "p_{t} #mu "         + post;
    muonpt = td->make<TH1F > ("ptMu", t.c_str(), 100, 0., 1000.);
    t = "#eta #mu "         + post;
    muoneta = td->make<TH1F > ("etaMu", t.c_str(), 50, -2.5, 2.5);
    t = "#phi #mu "         + post;
    muonphi = td->make<TH1F > ("phiMu", t.c_str(), 30, -3.14159, 3.14159);
    
    t = "p_{t} e "         + post;
    elecpt = td->make<TH1F > ("ptE", t.c_str(), 100, 0., 1000.);
    t = "#eta e "         + post;
    eleceta = td->make<TH1F > ("etaE", t.c_str(), 50, -2.5, 2.5);
    t = "#phi e "         + post;
    elecphi = td->make<TH1F > ("phiE", t.c_str(), 30, -3.14159, 3.14159);
}

// fill all histos of the set with the two lepton candidates

void HeavyNuTopHist::fill(HeavyNuEvent& hne)
{
    HeavyNuHistSet::fill(hne);

    double wgt = hne.eventWgt ;

    // Muons
    double mu1pt = hne.mu1.pt();
    double e1pt  = hnu::getElectronEt(hne.e1, false);

    mu1trackIso->Fill(hne.mu1.trackIso(), wgt);
    mu1hcalIso ->Fill(hne.mu1.hcalIso(), wgt);
    mu1ecalIso ->Fill(hne.mu1.ecalIso(), wgt);
    mu1caloIso ->Fill(hne.mu1.caloIso(), wgt);
    mu1dB      ->Fill(hne.mu1.dB(), wgt);
    e1trackIso->Fill(hne.e1.trackIso(), wgt);
    e1hcalIso ->Fill(hne.e1.hcalIso(), wgt);
    e1ecalIso ->Fill(hne.e1.ecalIso(), wgt);
    e1caloIso ->Fill(hne.e1.caloIso(), wgt);
    e1dB      ->Fill(hne.e1.dB(), wgt);

    mu1trackRelIso->Fill(hne.mu1.trackIso() / mu1pt, wgt);
    mu1hcalRelIso ->Fill(hne.mu1.hcalIso() / mu1pt, wgt);
    mu1ecalRelIso ->Fill(hne.mu1.ecalIso() / mu1pt, wgt);
    mu1caloRelIso ->Fill(hne.mu1.caloIso() / mu1pt, wgt);
    e1trackRelIso->Fill(hne.e1.trackIso() / e1pt, wgt);
    e1hcalRelIso ->Fill(hne.e1.hcalIso() / e1pt, wgt);
    e1ecalRelIso ->Fill(hne.e1.ecalIso() / e1pt, wgt);
    e1caloRelIso ->Fill(hne.e1.caloIso() / e1pt, wgt);

    for (int i = 0; i < muonQualityFlags; i++)
    {
        if (hne.mu1.muonID(muonQuality[i])) qualMu1->Fill( i, wgt ) ;
    }

    if(hne.e1.isEB())
    {
        if(hne.mu1.eta() < 1.44) mMuEBB->Fill(hne.mLL, wgt );
    }

    mWRvsminLPt->Fill(hne.mWR, std::min(hne.mu1.pt(), hne.e1.pt()), wgt);
    mWRvsNPV->Fill(hne.mWR, hne.n_primary_vertex, wgt);

    if(hne.e1.isEB())
    {
        mMuEBarrel->Fill(hne.mLL, wgt );
        if(hne.mu1.eta() < 1.44) mMuEBB->Fill(hne.mLL, wgt );

        mWRBarrel->Fill(hne.mWR, wgt);
        if(hne.mu1.eta() < 1.44) mWRBB->Fill(hne.mWR, wgt);
    }
    
    muonpt->Fill(hne.mu1.pt(), wgt);
    muoneta->Fill(hne.mu1.eta(), wgt);
    muonphi->Fill(hne.mu1.phi(), wgt);
    
    elecpt->Fill(hnu::getElectronEt(hne.e1, false));
    eleceta->Fill(hnu::getElectronSCEta(hne.e1));
    elecphi->Fill(hnu::getElectronSCPhi(hne.e1));
}