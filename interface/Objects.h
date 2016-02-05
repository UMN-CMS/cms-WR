class myElectron{
 public:
  TLorentzVector p4;
  Float_t scale;
  Float_t smearing;
  Int_t charge;
  
  myElectron(); ///< default contructor (empty)

};

myElectron::myElectron(){
}

class myMuon{
 public:
  TLorentzVector p4;
  Float_t IDSF;
  Float_t IsoSF;
  Float_t IDSF_error;
  Float_t IsoSF_error;
  Int_t charge;
  
  myMuon(); ///< default contructor (empty)

};

myMuon::myMuon(){
}

class myJet{
 public:
  TLorentzVector p4;
  Float_t jec_uncertainty;

  myJet(); ///< default contructor (empty)

};

myJet::myJet(){
}
