class miniTreeEvent{
    unsigned run;
    unsigned lumi;
    unsigned long long event;

    std::vector<TLorentzVector> electrons_p4;
    std::vector<TLorentzVector> muons_p4;
    std::vector<TLorentzVector> jets_p4;
    
    std::vector<float> jet_uncertainty;
    std::vector<float> electron_scale;
    std::vector<float> electron_smearing;

    float PU;
    float weight;
        
    miniTreeEvent(void) { clear(); };

  };


    void clear() {
      run = lumi = event = 0;
    
      electrons_p4.clear();
      muons_p4.clear();
      jets_p4.clear();

      jet_uncertainty.clear();
      electron_scale.clear();
      electron_smearing.clear();

      PU = -999.;
      weight = 0.0;
      
    }
