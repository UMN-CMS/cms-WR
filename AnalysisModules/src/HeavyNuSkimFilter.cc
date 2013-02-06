// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNu
//
/**\class HeavyNu HeavyNu.cc HeavyNu/AnalyzerModules/src/HeavyNu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNuSkimFilter.cc,v 1.7 2012/12/21 22:21:17 pastika Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
//#include <c++/4.1.2/bits/stl_vector.h>

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
// Order valid for 38X only.  Can be moved after Frameworkfwd.h in 39X
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
// Needed for 39X
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

class HeavyNuSkimFilter : public edm::EDFilter
{
public:
    explicit HeavyNuSkimFilter(const edm::ParameterSet&);
    ~HeavyNuSkimFilter();


private:

    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    edm::InputTag triggerResults_;
    std::vector<std::string> filterNames_;
} ;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HeavyNuSkimFilter::HeavyNuSkimFilter(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    triggerResults_ = iConfig.getParameter< edm::InputTag > ("triggerResults");
    filterNames_    = iConfig.getParameter< std::vector<std::string> > ("filterOnPaths");
}

HeavyNuSkimFilter::~HeavyNuSkimFilter(){

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------

bool HeavyNuSkimFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<edm::TriggerResults> trigResults;
    if (!iEvent.getByLabel(triggerResults_, trigResults))
        throw cms::Exception("TriggerResults collection does not exist") ; 

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*trigResults);

    bool trigger_fired = false;

    for (unsigned int i=0; i<filterNames_.size(); i++) { 
        int itrig = triggerNames.triggerIndex(filterNames_.at(i));
        if (trigResults->accept(itrig)) { trigger_fired = true ; break ; }
    }
    
    return trigger_fired ; 
}

// ------------ method called once each job just before starting event loop  ------------

void HeavyNuSkimFilter::beginJob(){ }

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNuSkimFilter::endJob(){ }

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuSkimFilter);
