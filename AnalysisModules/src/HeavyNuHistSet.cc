// Define default functions here
// Any histogram

#include "HeavyNu/AnalysisModules/src/HeavyNuHistSet.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

HeavyNuHistSet::HeavyNuHistSet(TFileDirectory *td, const std::string& post, int histsToBook)
{
    if(histsToBook & 1)
    {
        book(td, post);
    }
}

void HeavyNuHistSet::book(TFileDirectory *td, const std::string& post)
{
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    t = "event weights " + post;
    evtWeight = td->make<TH1F > ("evtWeight", t.c_str(), 1000, 0.0, 10.0);

    // ----------  Lepton histograms  ----------

    t = "p_{T}(#L_{1}) " + post;
    ptL1 = td->make<TH1F > ("ptL1", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(#L_{2}) " + post;
    ptL2 = td->make<TH1F > ("ptL2", t.c_str(), 100, 0., 1000.);
    t = "#eta(#L_{1}) " + post;
    etaL1 = td->make<TH1F > ("etaL1", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(#L_{1}, p_{T} > 30 GeV) " + post;
    etaL1pt30 = td->make<TH1F > ("etaL1pt30", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(#L_{1}, p_{T} > 40 GeV) " + post;
    etaL1pt40 = td->make<TH1F > ("etaL1pt40", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(#L_{2}) " + post;
    etaL2 = td->make<TH1F > ("etaL2", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(#L_{2}, p_{T} > 30 GeV) " + post;
    etaL2pt30 = td->make<TH1F > ("etaL2pt30", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(#L_{2}, p_{T} > 40 GeV) " + post;
    etaL2pt40 = td->make<TH1F > ("etaL2pt40", t.c_str(), 50, -2.5, 2.5);
    t = "#phi(#L_{1}) " + post;
    phiL1 = td->make<TH1F > ("phiL1", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#L_{1}, p_{T} > 30 GeV) " + post;
    phiL1pt30 = td->make<TH1F > ("phiL1pt30", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#L_{1}, p_{T} > 40 GeV) " + post;
    phiL1pt40 = td->make<TH1F > ("phiL1pt40", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#L_{2}) " + post;
    phiL2 = td->make<TH1F > ("phiL2", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#L_{2}, p_{T} > 30 GeV) " + post;
    phiL2pt30 = td->make<TH1F > ("phiL2pt30", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#L_{2}, p_{T} > 40 GeV) " + post;
    phiL2pt40 = td->make<TH1F > ("phiL2pt40", t.c_str(), 30, -3.14159, 3.14159);

    t = "p_{T}(#L_{1}) vs. p_{T}(#L_{2}) (SS) " + post + ";p_{T}(#L_{1})(GeV);p_{T}(#L_{2}(GeV))";

    ptL1VsPtL2ss = td->make<TH2F > ("ptL1VsPtL2ss", t.c_str(), 50, 0., 2000., 50, 0., 2000);

    t = "p_{T}(#L_{1}) vs. p_{T}(#L_{2}) (OS) " + post + ";p_{T}(#L_{1})(GeV);p_{T}(#L_{2}(GeV))";

    ptL1VsPtL2os = td->make<TH2F > ("ptL1VsPtL2os", t.c_str(), 50, 0., 2000., 50, 0., 2000);

    // delta angles

    t = "#Delta#eta(#L_{1},#L_{2}) " + post;
    dEtaL = td->make<TH1F > ("dEtaL", t.c_str(), 50, 0, 5);
    t = "#Delta#phi(#L_{1},#L_{2}) " + post;
    dPhiL = td->make<TH1F > ("dPhiL", t.c_str(), 30, 0, 3.14159);
    t = "#L #Delta#eta vs. #Delta#phi " + post;
    t += ";#Delta#eta; #Delta#phi";
    dEtaPhiL = td->make<TH2F > ("dEtaPhiL", t.c_str(), 50, 0, 5, 30, 0, 3.14159);

    // ----------  Jet histograms ----------

    t = "p_{T}(j_{1}) " + post;
    ptJet1 = td->make<TH1F > ("ptJet1", t.c_str(), 50, 0., 500.);
    t = "p_{T}(j_{2}) " + post;
    ptJet2 = td->make<TH1F > ("ptJet2", t.c_str(), 50, 0., 500.);
    t = "#eta(j_{1}) " + post;
    etaJet1 = td->make<TH1F > ("etaJet1", t.c_str(), 100, -5, 5);
    t = "#eta(j_{2}) " + post;
    etaJet2 = td->make<TH1F > ("etaJet2", t.c_str(), 100, -5, 5);
    t = "#phi(j_{1}) " + post;
    phiJet1 = td->make<TH1F > ("phiJet1", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(j_{2}) " + post;
    phiJet2 = td->make<TH1F > ("phiJet2", t.c_str(), 30, -3.14159, 3.14159);

    t = "#Delta#eta(j_{1},j_{2}) " + post;
    dEtaJet = td->make<TH1F > ("dEtaJet", t.c_str(), 100, 0, 5);
    t = "#Delta#phi(j_{1},j_{2}) " + post;
    dPhiJet = td->make<TH1F > ("dPhiJet", t.c_str(), 30, 0, 3.14159);

    t = "btag(j_{1}) " + post;
    btagJet1 = td->make<TH1F > ("btagJet1", t.c_str(), 40, 0, 5);
    t = "btag(j_{2}) " + post;
    btagJet2 = td->make<TH1F > ("btagJet2", t.c_str(), 40, 0, 5);

    t = "# B-tagged Jets in Event " + post;
    numBjets = td->make<TH1F > ("numBjets", t.c_str(), 3, -0.5, 2.5);
    t = "#  Jets in Event " + post;
    njets = td->make<TH1F > ("njets", t.c_str(), 200, -0.5, 199.5);

    t = "Jet #Delta#eta vs. #Delta#phi " + post + ";#Delta#eta; #Delta#phi";
    dEtaPhiJet = td->make<TH2F > ("dEtaPhiJet", t.c_str(), 50, 0, 5, 30, 0, 3.14159);
    t = "Jet ID(j_{2}) vs. ID(j_{1}) " + post + "; ID(j_{1}); ID(j_{2})";
    jetID2d = td->make<TH2I > ("jetID2d", t.c_str(), 3, 0, 3, 3, 0, 3);
    jetID2d->GetXaxis()->SetBinLabel(1, "Neither");
    jetID2d->GetXaxis()->SetBinLabel(2, "Loose");
    jetID2d->GetXaxis()->SetBinLabel(3, "Tight");
    jetID2d->GetYaxis()->SetBinLabel(1, "Neither");
    jetID2d->GetYaxis()->SetBinLabel(2, "Loose");
    jetID2d->GetYaxis()->SetBinLabel(3, "Tight");

    t = "min #DeltaR(jet, Wr genjet) " + post;
    mindRjet_genjet = td->make<TH1F > ("mindRjet_genjet", t.c_str(), 70, 0.0, 7.0);
    t = "max #DeltaR(jet, Wr genjet) " + post;
    maxdRjet_genjet = td->make<TH1F > ("maxdRjet_genjet", t.c_str(), 70, 0.0, 7.0);
    t = "number of Nu matched Jets " + post;
    nuLMatchedJets = td->make<TH1F > ("nuLMatchedJets", t.c_str(), 3, -0.5, 2.5);
    // ----------  MET histograms     ----------

    t = "MET distribution " + post;
    met = td->make<TH1F > ("met", t.c_str(), 100, 0, 2000);

    t = "MC Type " + post;
    mc_type = td->make<TH1F > ("mc_type", "MC Type Code", 100, -0.5, 99.5);

    t = "n(Pileup) " + post;
    n_pileup = td->make<TH1F > ("n_pileup", "Number of (MC) pileup", 50, -0.5, 49.5);
    t = "n(Vertex) " + post;
    n_vertex = td->make<TH1F > ("n_vertex", "Number of reconstructed vertices", 50, -0.5, 49.5);
    t = "n(Vertex), no weights " + post;
    n_vertex_noWgt = td->make<TH1F > ("n_vertex_noWgt", "Number of reconstructed vertices", 50, -0.5, 49.5);


    // ----------  L/Jet histograms  ----------

    t = "MiniLm #Delta R(#L,jet) " + post;
    dRminLJet = td->make<TH1F > ("dRminLJet", t.c_str(), 50, 0, 5.);
    t = "MiniLm #Delta R(#L_{1},jet) " + post;
    dRminL1jet = td->make<TH1F > ("dRminL1jet", t.c_str(), 50, 0, 5.);
    t = "MiniLm #Delta R(#L_{2},jet) " + post;
    dRminL2jet = td->make<TH1F > ("dRminL2jet", t.c_str(), 50, 0, 5.);

    t = "p_{T,rel}(#L_{1},jet)" + post;
    hptrelL1 = td->make<TH1F > ("ptrelL1", t.c_str(), 50, 0, 1000.);
    t = "p_{T,rel}(#L_{2},jet)" + post;
    hptrelL2 = td->make<TH1F > ("ptrelL2", t.c_str(), 50, 0, 1000.);

    t = "p_{T,rel}(#L_{1},jet) vs #Delta R(#L_{1},jet)" + post;
    t += "; #Delta R(#L_{1},jet); p_{T,rel}(#L_{1},jet)";
    ptrelVsdRminL1jet = td->make<TH2F > ("ptrelVsdRminL1jet", t.c_str(), 50, 0, 5., 50, 0, 1000);
    t = "p_{T,rel}(#L_{2},jet) vs #Delta R(#L_{2},jet)" + post;
    t += "; #Delta R(#L_{2},jet); p_{T,rel}(#L_{2},jet)";
    ptrelVsdRminL2jet = td->make<TH2F > ("ptrelVsdRminL2jet", t.c_str(), 50, 0, 5., 50, 0, 1000);

    // ----------  Composite histograms  ----------

    t = "M(W_{R}) " + post;
    mWR = td->make<TH1F > ("mWR", t.c_str(), 100, 0, 4000);
    t = "M(N_{R}) with #L_{1} " + post;
    mNuR1 = td->make<TH1F > ("mNuR1", t.c_str(), 100, 0, 4000);
    t = "M(N_{R}) with #L_{2} " + post;
    mNuR2 = td->make<TH1F > ("mNuR2", t.c_str(), 100, 0, 4000);
    t = "M(N_{R}) #L_{1} vs. #L_{2} " + post;
    mNuR2D = td->make<TH2F > ("mNuR2D", t.c_str(), 100, 0, 4000, 100, 0, 4000);

    t = "L_{1} p_{T} / M(W_{R}) " + post;
    L1ptFracWRmass = td->make<TH1F > ("L1ptFracWRmass", t.c_str(), 100, 0, 1.);

    t = "M(jj) " + post;
    mJJ = td->make<TH1F > ("mJJ", t.c_str(), 50, 0, 2000);
    t = "M(LL) " + post;
    mLL = td->make<TH1F > ("mLL", t.c_str(), 100, 0, 2000);
    t = "M(LL)(OS) " + post;
    mLLOS = td->make<TH1F > ("mLLOS", t.c_str(), 100, 0, 2000);
    t = "M(LL)(SS) " + post;
    mLLSS = td->make<TH1F > ("mLLSS", t.c_str(), 100, 0, 2000);
    t = "M(LL) " + post;
    mLLZoom = td->make<TH1F > ("mLLZoom", t.c_str(), 400, 0, 400);
    t = "DiLon Charge " + post;
    diLCharge = td->make<TH1F > ("diLCharge", t.c_str(), 2, -1, 1);

    t = "M(W_{R}) 1 BJet " + post;
    mWR_1b  = td->make<TH1F > ("mWR_1b", t.c_str(), 100, 0, 4000);
    t = "M(W_{R}) 2 BJet " + post;
    mWR_2b  = td->make<TH1F > ("mWR_2b", t.c_str(), 100, 0, 4000);
    t = "M(LL) 1 BJet " + post;
    mLL_1b = td->make<TH1F > ("mLL_1b", t.c_str(), 100, 0, 2000);
    t = "M(LL) 2 BJet " + post;
    mLL_2b = td->make<TH1F > ("mLL_2b", t.c_str(), 100, 0, 2000);
    t = "M(LL) 1 BJet " + post;
    mLLZoom_1b = td->make<TH1F > ("mLLZoom_1b", t.c_str(), 400, 0, 400);
    t = "M(LL) 2 BJet " + post;
    mLLZoom_2b = td->make<TH1F > ("mLLZoom_2b", t.c_str(), 400, 0, 400);

    diLCharge->GetXaxis()->SetBinLabel(1, "OS");
    diLCharge->GetXaxis()->SetBinLabel(2, "SS");

    t = "M(#L #L) vs. M_{W_{R}} " + post;
    mLLvsmWR = td->make<TH2F > ("mLLvsmWR", t.c_str(), 100, 0, 1000, 100, 0, 4000);
    t = "M_{W_{R}} vs. M_{N_{1}} " + post;
    mWRvsmNuR1 = td->make<TH2F > ("mWRvsmNuR1", t.c_str(), 100, 0, 4000, 100, 0, 4000);
    t = "M_{W_{R}} vs. M_{N_{2}} " + post;
    mWRvsmNuR2 = td->make<TH2F > ("mWRvsmNuR2", t.c_str(), 100, 0, 4000, 100, 0, 4000);
    t = "M_{W_{R}} vs. NPV " + post;
    mWRvsNPV = td->make<TH2F > ("mWRvsNPV", t.c_str(), 100, 0, 4000, 70, 0, 70);
    t = "M(#L #L) vs. NPV " + post;
    mLLZoomvsNPV = td->make<TH2F > ("mLLZoomvsNPV", t.c_str(), 400, 0, 400, 70, 0, 70);

    t = "mLLZoomvsptL1vsptL2 " + post;
    double xbins[] = {0.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 2000.0};
    double ybins[] = {0.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 2800.0};
    double zbins[] = {0.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 2800.0};
    mLLZoomvsptL1vsptL2 = td->make<TH3F > ("mLLZoomvsptL1vsptL2", t.c_str(), 14, xbins, 9, ybins, 7, zbins);

    // vertex histograms
    t = "LL vtx DeltaZ " + post;
    vtx_LL = td->make<TH1F > ("vtxMM", t.c_str(), 200, 0., 2.);
    t = "jj vtx DeltaZ " + post;
    vtx_jj = td->make<TH1F > ("vtxJJ", t.c_str(), 200, 0., 2.);
    t = "m1j vtx min DeltaZ " + post;
    vtx_min_L1j = td->make<TH1F > ("vtxM1Jmin", t.c_str(), 200, 0., 2.);
    t = "m2j vtx min DeltaZ " + post;
    vtx_min_L2j = td->make<TH1F > ("vtxM2Jmin", t.c_str(), 200, 0., 2.);
    t = "L-j vtx min DeltaZ " + post;
    vtx_min_Lj = td->make<TH1F > ("vtxLJmin", t.c_str(), 200, 0., 2.);
    t = "mmjj vtx max DeltaZ " + post;
    vtx_max_dist = td->make<TH1F > ("vtxDistmax", t.c_str(), 400, 0., 5.);

}// end of book()

void HeavyNuHistSet::fill(HeavyNuEvent& hne)
{
    double wgt = hne.eventWgt;

    if(hne.l1 == NULL && hne.l2 == NULL && hne.nLeptons >= 2)
    {
        std::cout << "INVALID LEPTON in HeavyNuHistSet fill" << std::endl;
    }

    evtWeight->Fill(wgt) ;
    mc_type->Fill(hne.mc_class, wgt);

    n_pileup->Fill(hne.n_pue, wgt) ;
    n_vertex->Fill(hne.n_primary_vertex, wgt) ;
    n_vertex_noWgt->Fill(hne.n_primary_vertex) ;

    // Leptons
    if(hne.nLeptons >= 1)
    {
        ptL1->Fill(hne.l1pt, wgt);
        etaL1->Fill(hne.l1eta, wgt);
        phiL1->Fill(hne.l1phi, wgt);
        if ( hne.l1pt > 30. )
        {
            etaL1pt30->Fill(hne.l1eta, wgt);
            phiL1pt30->Fill(hne.l1phi, wgt);
            if ( hne.l1pt > 40. )
            {
                etaL1pt40->Fill(hne.l1eta, wgt);
                phiL1pt40->Fill(hne.l1phi, wgt);
            }
        }
    }

    if(hne.nLeptons >= 2)
    {
        ptL2->Fill(hne.l2pt, wgt);
        etaL2->Fill(hne.l2eta, wgt);
        phiL2->Fill(hne.l2phi, wgt);
        if ( hne.l2pt > 30. )
        {
            etaL2pt30->Fill(hne.l2eta, wgt);
            phiL2pt30->Fill(hne.l2phi, wgt);
            if ( hne.l2pt > 40. )
            {
                etaL2pt40->Fill(hne.l2eta, wgt);
                phiL2pt40->Fill(hne.l2phi, wgt);
            }
        }

        dPhiL->Fill(fabs(deltaPhi(hne.l1phi, hne.l2phi)), wgt);
        dEtaL->Fill(fabs(hne.l1eta - hne.l2eta), wgt);
        dEtaPhiL->Fill(fabs(hne.l1eta - hne.l2eta),
                       fabs(deltaPhi(hne.l1phi, hne.l2phi)), wgt);

        vtx_LL->Fill(fabs(hne.l1->vertex().Z() - hne.l2->vertex().Z()), wgt);
    }

    int jet1id = 0;
    int jet2id = 0;

    // Jets
    njets->Fill(hne.nJets, wgt);
    if(hne.nJets > 0)
    {
        if(!hne.pfJets) jet1id = hnu::jetID(hne.j1);

        double j1bdisc = hne.j1.bDiscriminator(hne.btagName);

        ptJet1->Fill(hne.j1.pt(), wgt);
        etaJet1->Fill(hne.j1.eta(), wgt);
        phiJet1->Fill(hne.j1.phi(), wgt);
        btagJet1->Fill(j1bdisc, wgt);

        if(hne.nJets >= 2 && hne.nLeptons >= 2)
        {
            double j2bdisc = hne.j2.bDiscriminator(hne.btagName);

            ptJet2->Fill(hne.j2.pt(), wgt);
            etaJet2->Fill(hne.j2.eta(), wgt);
            phiJet2->Fill(hne.j2.phi(), wgt);
            btagJet2->Fill(j2bdisc, wgt);

            numBjets->Fill((int)hne.isBJet1 + (int)hne.isBJet2, wgt);

            if(!hne.pfJets) jet2id = hnu::jetID(hne.j2);

            dPhiJet->Fill(fabs(deltaPhi(hne.j1.phi(), hne.j2.phi())), wgt);
            dEtaJet->Fill(fabs(hne.j1.eta() - hne.j2.eta()), wgt);
            dEtaPhiJet->Fill(fabs(hne.j1.eta() - hne.j2.eta()),
                             fabs(deltaPhi(hne.j1.phi(), hne.j2.phi())), wgt);

            mWR->Fill(hne.mWR, wgt);
	    // Special work-around for top, where the 1st lepton is always the muon
	    // This ensures that the neutrino candidate with the highest pT lepton is always NuR1
	    if ( hne.mode == HeavyNuEvent::TOP ) { 
	      if ( hne.l1pt > hne.l2pt ) { 
		mNuR1->Fill(hne.mNuR1, wgt);
		mNuR2->Fill(hne.mNuR2, wgt);
		mNuR2D->Fill(hne.mNuR1, hne.mNuR2, wgt);
	      } else { 
		mNuR1->Fill(hne.mNuR2, wgt);
		mNuR2->Fill(hne.mNuR1, wgt);
		mNuR2D->Fill(hne.mNuR2, hne.mNuR1, wgt);
	      }	      
	    } else { 
	      mNuR1->Fill(hne.mNuR1, wgt);
	      mNuR2->Fill(hne.mNuR2, wgt);
	      mNuR2D->Fill(hne.mNuR1, hne.mNuR2, wgt);
	    }
            mJJ->Fill(hne.mJJ, wgt);

            L1ptFracWRmass->Fill(hne.l1pt / hne.mWR, wgt);

            float deltaVzJ1J2 = fabs(hne.tjV1 - hne.tjV2);
            float deltaVzJ1M1 = fabs(hne.tjV1 - hne.l1->vertex().Z());
            float deltaVzJ2M2 = fabs(hne.tjV2 - hne.l2->vertex().Z());
            float deltaVzJ1M2 = fabs(hne.tjV1 - hne.l2->vertex().Z());
            float deltaVzJ2M1 = fabs(hne.tjV2 - hne.l1->vertex().Z());
            float deltaVzM1M2 = fabs(hne.l1->vertex().Z() - hne.l2->vertex().Z());

            vtx_jj->Fill(deltaVzJ1J2, wgt);
            float minDeltaVzL1J = std::min(deltaVzJ1M1, deltaVzJ2M1);
            float minDeltaVzL2J = std::min(deltaVzJ2M2, deltaVzJ2M2);
            vtx_min_L1j->Fill(minDeltaVzL1J, wgt);
            vtx_min_L2j->Fill(minDeltaVzL2J, wgt);
            vtx_min_Lj->Fill(std::min(minDeltaVzL1J, minDeltaVzL2J), wgt);

            float maxDeltaVzLJ1 = std::max(deltaVzJ1M1, deltaVzJ1M2);
            float maxDeltaVzLJ2 = std::max(deltaVzJ2M1, deltaVzJ2M2);
            float maxDeltaVzMMJJ = std::max(deltaVzM1M2, deltaVzJ1J2);
            float maxDeltaVzLJ = std::max(maxDeltaVzLJ1, maxDeltaVzLJ2);
            vtx_max_dist->Fill(std::max(maxDeltaVzMMJJ, maxDeltaVzLJ), wgt);
        }

        dRminL1jet->Fill(hne.dRminL1jet, wgt);
        dRminL2jet->Fill(hne.dRminL2jet, wgt);

        dRminLJet->Fill(std::min(hne.dRminL1jet, hne.dRminL2jet), wgt);

        hptrelL1->Fill(hne.ptrelL1, wgt);
        hptrelL2->Fill(hne.ptrelL2, wgt);

        ptrelVsdRminL1jet->Fill(hne.dRminL1jet, hne.ptrelL1, wgt);
        ptrelVsdRminL2jet->Fill(hne.dRminL2jet, hne.ptrelL2, wgt);

        if(hne.isMC)
        {
            double dR1 = deltaR(hne.j1.p4(), hne.gj1.p4());
            double dR2 = deltaR(hne.j2.p4(), hne.gj2.p4());
            mindRjet_genjet->Fill(std::min(dR1, dR2), wgt);
            maxdRjet_genjet->Fill(std::max(dR1, dR2), wgt);
            nuLMatchedJets->Fill(hne.numNuLJetsMatched, wgt);
        }
    }

    if(!hne.pfJets) jetID2d->Fill(jet1id, jet2id, wgt);

    // met
    met->Fill(hne.met1.pt(), wgt);

    mLL->Fill(hne.mLL, wgt);

    if(hne.nLeptons >= 2)
    {
        if(hne.l1->charge() == hne.l2->charge())
        {
            mLLSS->Fill(hne.mLL, wgt);
            ptL1VsPtL2ss->Fill(hne.l1pt, hne.l2pt, wgt);
        }
        else
        {
            mLLOS->Fill(hne.mLL, wgt);
            ptL1VsPtL2os->Fill(hne.l1pt, hne.l2pt, wgt);
        }

        diLCharge->Fill(0.5 * hne.l1->charge() * hne.l2->charge(), wgt);
        mLLZoomvsptL1vsptL2->Fill(hne.mLL, hne.l1pt, hne.l2pt, wgt);
    }

    if(hne.isBJet1 || hne.isBJet2)
    {
        mWR_1b->Fill(hne.mWR, wgt);
        mLL_1b->Fill(hne.mLL, wgt);
        mLLZoom_1b->Fill(hne.mLL, wgt);
    }
    if(hne.isBJet1 && hne.isBJet2)
    {
        mWR_2b->Fill(hne.mWR, wgt);
        mLL_2b->Fill(hne.mLL, wgt);
        mLLZoom_2b->Fill(hne.mLL, wgt);
    }

    mLLZoom->Fill(hne.mLL, wgt);
    mLLvsmWR->Fill(hne.mLL, hne.mWR, wgt);
    // Special work-around for top, where the 1st lepton is always the muon
    // This ensures that the neutrino candidate with the highest pT lepton is always NuR1
    if ( hne.mode == HeavyNuEvent::TOP ) { 
      if ( hne.l1pt > hne.l2pt ) { 
	mWRvsmNuR1->Fill(hne.mWR, hne.mNuR1, wgt);
	mWRvsmNuR2->Fill(hne.mWR, hne.mNuR2, wgt);
      } else { 
	mWRvsmNuR1->Fill(hne.mWR, hne.mNuR2, wgt);
	mWRvsmNuR2->Fill(hne.mWR, hne.mNuR1, wgt);
      }
    } else { 
      mWRvsmNuR1->Fill(hne.mWR, hne.mNuR1, wgt);
      mWRvsmNuR2->Fill(hne.mWR, hne.mNuR2, wgt);
    }
    mWRvsNPV->Fill(hne.mWR, hne.n_primary_vertex, wgt);
    mLLZoomvsNPV->Fill(hne.mLL, hne.n_primary_vertex, wgt);


}// end of fill()
