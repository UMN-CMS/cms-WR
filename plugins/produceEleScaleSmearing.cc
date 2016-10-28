#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <memory>
#include <iomanip>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"

//#include "TLorentzVector.h"
#include <vector>
using namespace std;

class produceEleScaleSmearing : public edm::EDProducer
{


  typedef float v_t;
  typedef edm::ValueMap<v_t> scale_factors_Map;

public:
  produceEleScaleSmearing(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  edm::EDGetToken srcToken_;
  std::string outputCollName1;     ///<label name of collection made by this producer
  std::string outputCollName2;     ///<label name of collection made by this producer
  std::string outputCollName3;     ///<label name of collection made by this producer
  std::string outputCollName4;     ///<label name of collection made by this producer
  std::string outputCollName5;     ///<label name of collection made by this producer
  std::string outputCollName6;     ///<label name of collection made by this producer

};

produceEleScaleSmearing::produceEleScaleSmearing(const edm::ParameterSet& cfg)
  : srcToken_ ( consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("electrons_src"))),
    outputCollName1(cfg.getParameter<std::string>("OutputCollectionName1")),
    outputCollName2(cfg.getParameter<std::string>("OutputCollectionName2")),
    outputCollName3(cfg.getParameter<std::string>("OutputCollectionName3")),
    outputCollName4(cfg.getParameter<std::string>("OutputCollectionName4")),
    outputCollName5(cfg.getParameter<std::string>("OutputCollectionName5")),
    outputCollName6(cfg.getParameter<std::string>("OutputCollectionName6"))
{
  produces<scale_factors_Map>(outputCollName1);
  produces<scale_factors_Map>(outputCollName2);
  produces<scale_factors_Map>(outputCollName3);
  produces<scale_factors_Map>(outputCollName4);
  produces<scale_factors_Map>(outputCollName5);
  produces<scale_factors_Map>(outputCollName6);
}

void produceEleScaleSmearing::produce(edm::Event& event, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<pat::Electron> > electrons;
  event.getByToken(srcToken_, electrons);

  EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV1_2016_ele");

  std::vector<v_t>  scaleError;
  std::vector<v_t>  smearingSigma;
  std::vector<v_t>  smearingSigma_phi_up;
  std::vector<v_t>  smearingSigma_phi_down;
  std::vector<v_t>  smearingSigma_rho_up;
  std::vector<v_t>  smearingSigma_rho_down;

  std::auto_ptr<scale_factors_Map> scaleErrorMap(new scale_factors_Map());
  std::auto_ptr<scale_factors_Map> smearingSigmaMap(new scale_factors_Map());
  std::auto_ptr<scale_factors_Map> smearingSigma_phi_upMap(new scale_factors_Map());
  std::auto_ptr<scale_factors_Map> smearingSigma_phi_downMap(new scale_factors_Map());
  std::auto_ptr<scale_factors_Map> smearingSigma_rho_upMap(new scale_factors_Map());
  std::auto_ptr<scale_factors_Map> smearingSigma_rho_downMap(new scale_factors_Map());

  for (const auto& ele: *electrons) {
    if(event.isRealData()) {
      float error_scale = eScaler.ScaleCorrectionUncertainty(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et());

      scaleError.push_back(error_scale);
      smearingSigma.push_back(0.);
      smearingSigma_phi_up.push_back(0.);
      smearingSigma_phi_down.push_back(0.);
      smearingSigma_rho_up.push_back(0.);
      smearingSigma_rho_down.push_back(0.);

    } else {
      float sigma = eScaler.getSmearingSigma(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et(), 0, 0);
      float sigma_phi_up = eScaler.getSmearingSigma(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et(), 0, 1);
      float sigma_phi_down = eScaler.getSmearingSigma(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et(), 0, -1);
      float sigma_rho_up = eScaler.getSmearingSigma(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et(), 1, 0);
      float sigma_rho_down = eScaler.getSmearingSigma(event.id().run(), ele.isEB(), ele.r9(), ele.superCluster()->eta(), ele.et(), -1, 0);

      scaleError.push_back(0.);
      smearingSigma.push_back(sigma);
      smearingSigma_phi_up.push_back(sigma_phi_up);
      smearingSigma_phi_down.push_back(sigma_phi_down);
      smearingSigma_rho_up.push_back(sigma_rho_up);
      smearingSigma_rho_down.push_back(sigma_rho_down);

    }
  }

  scale_factors_Map::Filler scaleError_filler(*scaleErrorMap);
  scale_factors_Map::Filler smearingSigma_filler(*smearingSigmaMap);
  scale_factors_Map::Filler smearingSigma_phi_up_filler(*smearingSigma_phi_upMap);
  scale_factors_Map::Filler smearingSigma_phi_down_filler(*smearingSigma_phi_downMap);
  scale_factors_Map::Filler smearingSigma_rho_up_filler(*smearingSigma_rho_upMap);
  scale_factors_Map::Filler smearingSigma_rho_down_filler(*smearingSigma_rho_downMap);

  scaleError_filler.insert(electrons, scaleError.begin(), scaleError.end());
  smearingSigma_filler.insert(electrons, smearingSigma.begin(), smearingSigma.end());
  smearingSigma_phi_up_filler.insert(electrons, smearingSigma_phi_up.begin(), smearingSigma_phi_up.end());
  smearingSigma_phi_down_filler.insert(electrons, smearingSigma_phi_down.begin(), smearingSigma_phi_down.end());
  smearingSigma_rho_up_filler.insert(electrons, smearingSigma_rho_up.begin(), smearingSigma_rho_up.end());
  smearingSigma_rho_down_filler.insert(electrons, smearingSigma_rho_down.begin(), smearingSigma_rho_down.end());

  scaleError_filler.fill();
  smearingSigma_filler.fill();
  smearingSigma_phi_up_filler.fill();
  smearingSigma_phi_down_filler.fill();
  smearingSigma_rho_up_filler.fill();
  smearingSigma_rho_down_filler.fill();

  event.put(scaleErrorMap, outputCollName1);
  event.put(smearingSigmaMap, outputCollName2);
  event.put(smearingSigma_phi_upMap, outputCollName3);
  event.put(smearingSigma_phi_downMap, outputCollName4);
  event.put(smearingSigma_rho_upMap, outputCollName5);
  event.put(smearingSigma_rho_downMap, outputCollName6);

  scaleError.clear();
  smearingSigma.clear();
  smearingSigma_phi_up.clear();
  smearingSigma_phi_down.clear();
  smearingSigma_rho_up.clear();
  smearingSigma_rho_down.clear();

}

DEFINE_FWK_MODULE(produceEleScaleSmearing);
