// -*- C++ -*-
//
// Package:    LumiList
// Class:      LumiList
// 
/**\class LumiList LumiList.cc HeavyNu/LumiList/src/LumiList.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nathaniel Pastika
//         Created:  Fri Nov 18 08:08:01 CST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class LumiList : public edm::EDAnalyzer {
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      struct lumiStruct
      {
          unsigned int run;
          unsigned int llumi;
          unsigned int ulumi;
      } lumientry;
      
      std::vector<lumiStruct > prelumijson;
      bool first;
      
   public:
      explicit LumiList(const edm::ParameterSet&);
      ~LumiList();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      friend bool compLumiStruct(LumiList::lumiStruct a, LumiList::lumiStruct b);
};

bool compLumiStruct(LumiList::lumiStruct a, LumiList::lumiStruct b)
{
    return (a.run == b.run)?(a.ulumi < b.ulumi):(a.run < b.run);
}

LumiList::LumiList(const edm::ParameterSet& iConfig)
{
    first = true;
    lumientry.run = 0;
    lumientry.llumi = 0;
    lumientry.ulumi = 0;
}


LumiList::~LumiList()
{
}


void LumiList::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}


void LumiList::beginJob()
{
}

void LumiList::endJob() 
{
    std::sort(prelumijson.begin(), prelumijson.end(), compLumiStruct);
    std::cout << "~~LUMIDUMP~~\n";
    for(std::vector<LumiList::lumiStruct>::const_iterator i = prelumijson.begin(); i != prelumijson.end(); i++)
    {
        std::cout << "\"" << i->run << ":" << i->llumi << "-" << i->ulumi << "\"" << std::endl;
    }
    std::cout << "~~NO MORE LUMIDUMP~~\n";
}

void LumiList::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

void LumiList::endRun(edm::Run const&, edm::EventSetup const&)
{
    prelumijson.push_back(lumientry);
    first = true; 
}

void LumiList::beginLuminosityBlock(edm::LuminosityBlock const& lb, edm::EventSetup const&)
{
    unsigned int locRunNum = (unsigned int)lb.run();
    unsigned int locLumiSec = (unsigned int)lb.luminosityBlock();
    if(first)
    {
        lumientry.run = locRunNum;
        lumientry.llumi = locLumiSec;
        first = false;
    }
    else if(locLumiSec - lumientry.ulumi > 1)
    {
        prelumijson.push_back(lumientry);
        lumientry.run = locRunNum;
        lumientry.llumi = locLumiSec;
    }
    lumientry.ulumi = locLumiSec;
}

void LumiList::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void LumiList::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LumiList);
