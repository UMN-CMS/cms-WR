// How to execute:
// "linux> root fall10zjets.root
//  root [0] 
//  Attaching file fall10zjets.root as _file0...
//  root [1] .x residualJESMES.C(_file0)"
//
void residualJESMES(TFile* zjetsf0)
{
  // To get one or the other, global replace JES <-> MES
  // from here on down and re-execute
  //
  TH1* hnm6 = zjetsf0->Get("hNu/cut6_Mu1HighPt/mMuMuZoom");
  TH1* hhi6 = zjetsf0->Get("hNuJEShi/cut6_Mu1HighPt/mMuMuZoom");
  TH1* hlo6 = zjetsf0->Get("hNuJESlo/cut6_Mu1HighPt/mMuMuZoom");

  TH1* hnm7 = zjetsf0->Get("hNu/cut7_diLmass/mMuMu");
  TH1* hhi7 = zjetsf0->Get("hNuJEShi/cut7_diLmass/mMuMu");
  TH1* hlo7 = zjetsf0->Get("hNuJESlo/cut7_diLmass/mMuMu");

  TH1* hnm8 = zjetsf0->Get("hNu/cut8_mWRmass/mMuMu");
  TH1* hhi8 = zjetsf0->Get("hNuJEShi/cut8_mWRmass/mMuMu");
  TH1* hlo8 = zjetsf0->Get("hNuJESlo/cut8_mWRmass/mMuMu");

  // Integrals of nominal/hi/lo at cut levels 6,7,8
  //
  double cut6nm = hnm6->Integral(), cut6hi = hhi6->Integral(), cut6lo = hlo6->Integral();
  double cut7nm = hnm7->Integral(), cut7hi = hhi7->Integral(), cut7lo = hlo7->Integral();
  double cut8nm = hnm8->Integral(), cut8hi = hhi8->Integral(), cut8lo = hlo8->Integral();

  // Ratios of hi/lo to nominal
  //
  double rHi2nom6 = cut6hi/cut6nm, rLo2nom6 = cut6lo/cut6nm; //  at cut level 6
  double rHi2nom7 = cut7hi/cut7nm, rLo2nom7 = cut7lo/cut7nm; //  at cut level 7
  double rHi2nom8 = cut8hi/cut8nm, rLo2nom8 = cut8lo/cut8nm; //  at cut level 8

  // Because we chi2-fit the Z-peak to data in the nominal case, we
  // can theoretically do the same for the high and low-fluctuated
  // systematic cases. Therefore the important quantity for Z+Jets is
  // not the ratio of #events remaining (lo/nom,hi/nom), since a
  // renormalization at the Z-peak would largely eliminate that;
  // rather it is the change in the shape of the distribution caused by
  // the energy scale fluctuation. For that we look at how the cut flow
  // changes for high and low fluctuations ("cutn/cut6") and take that
  // as our systematic.
  //
  // cut6 acts as our baseline, because it includes all of our jet/mu selections,
  //     as well as the Zpeak before the MLL cut (cut7) of 200GeV;
  //     it is the Zpeak we would use to normalize to data.
  // 
  printf("cut:   IntNom\tJESlo\t(lo/nom) (cutn/cut6)\tJEShi\t(hi/nom) (cutn/cut6)\n");

  printf("6  :   %6.1f\t%6.1f\t(%6.4f)\t\t%6.1f\t(%6.4f)\n",
	 cut6nm,cut6lo,rLo2nom6,cut6hi,rHi2nom6);

  printf("7  :   %6.4f\t%6.4f\t(%6.4f) (%6.4f)\t%6.4f\t(%6.4f) (%6.4f)\n",
	 cut7nm,cut7lo,rLo2nom7,rLo2nom7/rLo2nom6,cut7hi,rHi2nom7,rHi2nom7/rHi2nom6);

  printf("8  :   %6.4f\t%6.4f\t(%6.4f) (%6.4f)\t%6.4f\t(%6.4f) (%6.4f)\n",
	 cut8nm,cut8lo,rLo2nom8,rLo2nom8/rLo2nom6,cut8hi,rHi2nom8,rHi2nom8/rHi2nom6);
}
