#ifndef myobjects_h
#define myobjects_h

class myElectron
{
public:
	TLorentzVector p4;
	Float_t scale;
	Float_t smearing;
	Int_t charge;
	Float_t r9;
	Float_t weight;
	Float_t IDSF;
	Float_t RecoSF;
	Float_t HltSF;
	Float_t IDSF_error;
	Float_t RecoSF_error;
	Float_t HltSF_error;
	Float_t scale_error;
	Float_t smearing_error;

	myElectron() {}; ///< default contructor (empty)

};

/* myElectron::myElectron(){ */
/* } */

class myMuon
{
public:
	TLorentzVector p4;
	Float_t IDSF;
	Float_t IsoSF;
	Float_t IDSF_error;
	Float_t IsoSF_error;
	Int_t charge;
	Float_t weight;

	myMuon() {}; ///< default contructor (empty)

};

class myJet
{
public:
	TLorentzVector p4;
	Float_t jec_uncertainty;
	Float_t weight;
	Float_t JER;
	Float_t JER_sf;
	Float_t JER_sf_up;
	Float_t JER_sf_down;


	myJet() {}; ///< default contructor (empty)

};



typedef std::vector<myElectron> myElectronCollection;
typedef std::vector<myMuon> myMuonCollection;
typedef std::vector<myJet> myJetCollection;


#endif
