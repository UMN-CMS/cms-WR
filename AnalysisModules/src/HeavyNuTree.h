#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"

#ifndef HeavyNuTree_h_included
#define HeavyNuTree_h_included

class HeavyNuTree
{
private:
    TDirectory& m_file;
    TTree* m_tree;
    TBranch* branch_;

public:

    struct HNuSlopeFitInfo
    {
        Float_t mlljj, mll, weight, l1pt, l2pt;
        Short_t flavor, cutlevel, n_pileup, n_primaryVertex;
    } event_;

    HeavyNuTree(TDirectory& f, bool writable) : m_file(f)
    {
        if(writable)
        {
            m_file.cd();
            m_tree = new TTree("HeavyNuTuple", "HeavyNuTuple");
            branch_ = m_tree->Branch("slopefit", &event_, "mlljj/F:mll:weight:l1pt:l2pt:flavor/S:cutlevel:npu:npv");
        }
        else
        {
            m_tree = (TTree*) m_file.Get("HeavyNuTuple");
            m_tree->SetBranchAddress("slopefit", &event_);
        }
    }

    void clear()
    {
        event_.mlljj = event_.mll = event_.weight = 0.0;
        event_.flavor = event_.cutlevel = event_.n_pileup = event_.n_primaryVertex = 0;
    }

    void fill()
    {
        m_tree->Fill();
    }
};

#endif
