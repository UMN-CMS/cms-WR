
double getSignalXSec(int mw, int mnu) {
  FILE *ff=fopen("signalXSec.csv","r");
  char buffer[100];
  double rv=-1;
  const double kfactor=1.13;
  float ax;
  int rmw,rmnu;

  while (!feof(ff)) {
    buffer[0]=0;
    fgets(buffer,90,ff);
    int igot=sscanf(buffer,"%d,%d,%f",&rmw,&rmnu,&ax);
    if (igot==3 && mw==rmw && rmnu==mnu) {
      fclose(ff);
      return ax*1e9*kfactor; // mb -> pb
    }
  } 
  fclose(ff);
  return rv;
  

}

void plotLimits(const char* f,int mw) {

  gROOT->SetStyle("Plain");

  int n=0;
  double mn[50],obs[50],exp[50],exp_p1s[50],exp_p2s[50],exp_m1s[50],exp_m2s[50];
  double theory[50];

  int rmw,rmn;
  float rexp,rexp_p1,rexp_p2,rexp_m1,rexp_m2,robs;
  char buffer[100];

  FILE *ff=fopen(f,"r");

  while (!feof(ff)) {
    buffer[0]=0;
    fgets(buffer,90,ff);
    int igot=sscanf(buffer,"%d,%d,%f,%f,%f,%f,%f,%f",&rmw,&rmn,
		    &rexp,&rexp_p2,&rexp_p1,&rexp_m1,&rexp_m2,&robs);
    
    if (igot==8 && rmw==mw) {
      // this part gets the bits in order
      int pos;
      for (pos=0; pos<n; pos++) {
	if (rmn<mn[pos]) break;
      }
      // move items as needed
      for (int j=n-1; j>=pos; j--) {
	mn[j+1]=mn[j];
	obs[j+1]=obs[j];
	exp[j+1]=exp[j];
	exp_p2s[j+1]=exp_p2s[j];
	exp_p1s[j+1]=exp_p1s[j];
	exp_m1s[j+1]=exp_m1s[j];
	exp_m2s[j+1]=exp_m2s[j];
	theory[j+1]=theory[j];
      }
      
      mn[pos]=rmn;
      obs[pos]=robs;
      exp[pos]=rexp;
      exp_p1s[pos]=rexp_p1;
      exp_p2s[pos]=rexp_p2;
      exp_m1s[pos]=rexp_m1;
      exp_m2s[pos]=rexp_m2;
      theory[pos]=getSignalXSec(mw,rmn);
      n++;
    }
  }

  fclose(ff);

  TCanvas* c1=new TCanvas("c1","c1",500,500);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetLogy();

  double mn_s[120],exp_1s[120], exp_2s[120];

  mn[n]=mn[n-1];
  mn[n+1]=mn[0];

  for (int ii=0; ii<n; ii++) {
    mn_s[ii]=mn[ii];
    exp_1s[ii]=exp_p1s[ii];
    exp_2s[ii]=exp_p2s[ii];
  }

  for (int ii=0; ii<n; ii++) {
    mn_s[2*n-ii-1]=mn[ii];
    exp_1s[2*n-ii-1]=exp_m1s[ii];
    exp_2s[2*n-ii-1]=exp_m2s[ii];
  }

  TGraph* graphObs=new TGraph(n,mn,obs);
  TGraph* graphExp=new TGraph(n,mn,exp);
  TGraph* graphExp2s=new TGraph(2*n,mn_s,exp_2s);
  TGraph* graphExp1s=new TGraph(2*n,mn_s,exp_1s);
  TGraph* graphTh=new TGraph(n,mn,theory);

  TH2* dummy=new TH2F("dummy","",40,100,1500,30,0.01,2.0);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR( W_{R} #rightarrow #mu + #mu + 2j ) [pb]");
  dummy->GetXaxis()->SetTitle("M(#nu_{heavy}) [GeV]");
  dummy->Draw();

  graphExp2s->SetFillColor(kYellow);
  graphExp2s->Draw("SAMEF");

  graphExp1s->SetFillColor(kGreen);
  graphExp1s->Draw("SAMEF");
  
  graphObs->SetLineWidth(2);
  graphObs->Draw("SAMEL");

  graphExp->SetLineWidth(2);
  graphExp->SetLineStyle(2);
  graphExp->SetMarkerStyle(20);
  graphExp->Draw("SAMELP");

  graphTh->SetLineWidth(2);
  graphTh->SetLineColor(kRed);
  graphTh->Draw("SAMEL");

  char text[20];
  sprintf(text,"M_{WR} = %.1f TeV",mw/1000.0f);
  TLegend* tl=new TLegend(0.6,0.98,0.98,0.75,text);
  tl->AddEntry(graphExp,"95% CL Expected Limit","L");
  tl->AddEntry(graphExp1s,"#pm1#sigma Expected Limit","F");
  tl->AddEntry(graphExp2s,"#pm2#sigma Expected Limit","F");
  tl->AddEntry(graphObs,"95% CL Observed Limit","L");
  tl->AddEntry(graphTh,"Pythia Signal Xsec","L");
  tl->Draw();

  char ofn[120];
  sprintf(ofn,"lim_mw%d.png",mw);
  c1->Print(ofn);

}
