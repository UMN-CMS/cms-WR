void plotAllObs(const char* f) {

  gROOT->SetStyle("Plain");

  int nm=0;
  int mw[10],nmw[10];
  double mn[10][50],obs[10][50],exp[10][50],exp_p1s[10][50],exp_p2s[10][50];

  int rmw,rmn;
  float rexp,rexp_p1,rexp_p2,rexp_m1,rexp_m2,robs;
  char buffer[100];

  FILE *ff=fopen(f,"r");
  while (!feof(ff)) {
    buffer[0]=0;
    fgets(buffer,90,ff);
    int igot=sscanf(buffer,"%d,%d,%f,%f,%f,%f,%f,%f",&rmw,&rmn,
		    &rexp,&rexp_p2,&rexp_p1,&rexp_m1,&rexp_m2,&robs);
    
    if (igot==8) {
      int ix=-1;
      for (int ij=0; ij<nm; ij++) {
	if (mw[ij]==rmw) { 
	  ix=ij; 
	  break;
	}
      }
      if (ix==-1) {
	ix=nm;
	mw[ix]=rmw;
	nmw[ix]=0;
	nm++;
      }
      int n=nmw[ix];
      mn[ix][n]=rmn;
      obs[ix][n]=robs;
      exp[ix][n]=rexp;
      exp_p1s[ix][n]=rexp_p1;
      exp_p2s[ix][n]=rexp_p2;
      nmw[ix]++;
    }
  }

  fclose(ff);

  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetLogy();

  TGraph* obsg[10];
  for (int i=0; i<nm; i++)
    obsg[i]=new TGraph(nmw[i],mn[i],obs[i]);

  TH2* dummy=new TH2F("dummy","",40,100,1500,30,0.1,2.0);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR( W_{R} #rightarrow #mu + #mu + 2j ) [pb]");
  dummy->GetXaxis()->SetTitle("M(#nu_{heavy}) [GeV]");
  dummy->Draw();

  TLegend* tl=new TLegend(0.6,0.9,0.9,0.6);

  for (int i=0; i<nm; i++) {
    obsg[i]->SetLineWidth(2);
    obsg[i]->SetLineColor(i+1);
    obsg[i]->Draw("L");
    char* text=new char[30];
    sprintf(text,"M_{WR} = %d GeV",mw[i]);
    tl->AddEntry(obsg[i], text,"L");
  }

  tl->Draw("SAME");

  char ofn[120];
  sprintf(ofn,"lim_obs.png",mw);
  c1->Print(ofn);

}
