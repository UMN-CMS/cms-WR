#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <string>
#include "TLorentzVector.h"
#include <vector>
using namespace std;

int FindBin(double Grid[],double Value,const int size);
const int Bins=14;
double Grid[Bins]={-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,-0.2,0.2,0.3,0.9,1.2,1.6,2.1,2.4};

class produceLepIdIsoScaleFactors : public edm::EDProducer{


typedef float v_t;
typedef edm::ValueMap<v_t> scale_factors_Map;
  
public:
	produceLepIdIsoScaleFactors(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);
	TLorentzVector Mu_Original;
private:
	const edm::InputTag src;	///<input particle objects
        std::string outputCollName1;     ///<label name of collection made by this producer
        std::string outputCollName2;     ///<label name of collection made by this producer
        std::string outputCollName3;     ///<label name of collection made by this producer
        std::string outputCollName4;     ///<label name of collection made by this producer

        std::vector<double> Scale_Factor_ID_Central;
        std::vector<double> Scale_Factor_ID_Error; 
        std::vector<double> Scale_Factor_ISO_Central;
        std::vector<double> Scale_Factor_ISO_Error;
//        double Scale_Factor_Central;
//        double Scale_Factor_Error;

};

produceLepIdIsoScaleFactors::produceLepIdIsoScaleFactors(const edm::ParameterSet& cfg)
	: src(cfg.getParameter<edm::InputTag>("src")),
          outputCollName1(cfg.getParameter<std::string>("OutputCollectionName1")),
          outputCollName2(cfg.getParameter<std::string>("OutputCollectionName2")),
          outputCollName3(cfg.getParameter<std::string>("OutputCollectionName3")),
          outputCollName4(cfg.getParameter<std::string>("OutputCollectionName4")),
          Scale_Factor_ID_Central(cfg.getParameter<std::vector<double>>("Scale_Factor_ID_Central")),
          Scale_Factor_ID_Error(cfg.getParameter<std::vector<double>>("Scale_Factor_ID_Error")),
          Scale_Factor_ISO_Central(cfg.getParameter<std::vector<double>>("Scale_Factor_ISO_Central")),
          Scale_Factor_ISO_Error(cfg.getParameter<std::vector<double>>("Scale_Factor_ISO_Error"))
{
        produces<scale_factors_Map>(outputCollName1);
        produces<scale_factors_Map>(outputCollName2);
        produces<scale_factors_Map>(outputCollName3);
        produces<scale_factors_Map>(outputCollName4);
}

void produceLepIdIsoScaleFactors::produce(edm::Event& event, const edm::EventSetup&)
{
	edm::Handle<pat::MuonCollection> muons;
	event.getByLabel(src, muons);
	std::auto_ptr<pat::MuonCollection> mus(new pat::MuonCollection);
        std::vector<v_t>  scale_factor_ID_central;
        std::vector<v_t>  scale_factor_ID_error;
        std::vector<v_t>  scale_factor_ISO_central;
        std::vector<v_t>  scale_factor_ISO_error;

       std::auto_ptr<scale_factors_Map> scale_factor_ID_centralMap(new scale_factors_Map());
       std::auto_ptr<scale_factors_Map> scale_factor_ID_errorMap(new scale_factors_Map());
       std::auto_ptr<scale_factors_Map> scale_factor_ISO_centralMap(new scale_factors_Map());
       std::auto_ptr<scale_factors_Map> scale_factor_ISO_errorMap(new scale_factors_Map());

	for(auto mu : *muons) {
            Mu_Original.SetPtEtaPhiE(mu.pt(),mu.eta(),mu.phi(),mu.energy());
            if(!event.isRealData()){
               if(FindBin(Grid,mu.eta(),Bins) != -1){
                  scale_factor_ID_central.push_back(Scale_Factor_ID_Central[FindBin(Grid,mu.eta(),Bins)]);
                  scale_factor_ID_error.push_back(Scale_Factor_ID_Error[FindBin(Grid,mu.eta(),Bins)]);
                  scale_factor_ISO_central.push_back(Scale_Factor_ISO_Central[FindBin(Grid,mu.eta(),Bins)]);
                  scale_factor_ISO_error.push_back(Scale_Factor_ISO_Error[FindBin(Grid,mu.eta(),Bins)]);

                 }
               else{
                    scale_factor_ID_central.push_back(-1);
                    scale_factor_ID_error.push_back(-1);
                    scale_factor_ISO_central.push_back(-1);
                    scale_factor_ISO_error.push_back(-1);
                     } 
              }		
	}

       scale_factors_Map::Filler scale_factor_ID_central_filler(*scale_factor_ID_centralMap);
       scale_factors_Map::Filler scale_factor_ID_error_filler(*scale_factor_ID_errorMap);
       scale_factors_Map::Filler scale_factor_ISO_central_filler(*scale_factor_ISO_centralMap);
       scale_factors_Map::Filler scale_factor_ISO_error_filler(*scale_factor_ISO_errorMap);


       scale_factor_ID_central_filler.insert(muons,scale_factor_ID_central.begin(),scale_factor_ID_central.end());
       scale_factor_ID_error_filler.insert(muons,scale_factor_ID_error.begin(),scale_factor_ID_error.end());
       scale_factor_ISO_central_filler.insert(muons,scale_factor_ISO_central.begin(),scale_factor_ISO_central.end());
       scale_factor_ISO_error_filler.insert(muons,scale_factor_ISO_error.begin(),scale_factor_ISO_error.end());


       event.put(scale_factor_ID_centralMap, outputCollName1);
       event.put(scale_factor_ID_errorMap, outputCollName2);
       event.put(scale_factor_ISO_centralMap, outputCollName3);
       event.put(scale_factor_ISO_errorMap, outputCollName4);

       scale_factor_ID_central.clear();
       scale_factor_ID_error.clear();
       scale_factor_ISO_central.clear();
       scale_factor_ISO_error.clear();
     
}

int FindBin(double Grid[],double Value,const int size){
for(int iii=0;iii<size-1;iii++){
       if(Value>=Grid[iii] && Value<Grid[iii+1]){
          if(iii < 7) return (6-iii);
          else return iii; 
         }
     }
    return -1;
   }


DEFINE_FWK_MODULE(produceLepIdIsoScaleFactors);
