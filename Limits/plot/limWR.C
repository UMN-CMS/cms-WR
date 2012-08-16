double accRatio(int mwr) {
  return 1.0;
  // (from the acceptance db) for generated versus MN=MW/2
  int roundmwr=(mwr/100)*100;
  switch (mwr) {
  case (1000) : return 0.605/0.630;
  case (1100) : return 0.637/0.657;
  case (1200) : return 0.678/0.691;
  case (1300) : return 0.696/0.708;
  case (1400) : return 0.713/0.729;
  case (1500) : return 0.729/0.742;
  case (1600) : return 0.740/0.759;
  case (1700) : return 0.749/0.767;
  case (1800) : return 0.766/0.774;
  case (1900) : return 0.772/0.774;
  case (2000) : return 0.778/0.789;
  case (2100) : return 0.791/0.794;
  case (2200) : return 0.794/0.793;
  case (2300) : return 0.787/0.800;
  case (2400) : return 0.799/0.808;
  case (2500) : return 0.803/0.805;
  default: 
    return accRatio(roundmwr)+(mwr-(roundmwr+100))*(accRatio(roundmwr+100)-accRatio(roundmwr))/100.0;
    //return 1.0;
  }
}

void smoothSG(float pts[], int n, int mode=4) {
  // − 3y1 + 12y2 + 17y3 + 12y4 − 3y5) / 35

  float ptsOut[1000];
 
  //  const int mode=19;
  mode+=10;

  if (mode==1) {
    for (int i=0; i<n; i++) {
      double sum=0;
      
      if (i>1) sum+=-3*pts[i-2];
      else if (i>0) sum+=-3*pts[i-1];
      else sum+=-3*pts[i];
      
      if (i>0) sum+=12*pts[i-1];
      else sum+=12*pts[i];
      
      sum+=17*pts[i];
      
      if (i<n-1) sum+=12*pts[i+1];
      else sum+=12*pts[i];
      
      if (i<n-2) sum+=-3*pts[i+2];
      else if (i<n-1) sum+=-3*pts[i+1];
      else sum+=-3*pts[i];
      
      sum/=35;
      ptsOut[i]=sum;
    }
  }
  if (mode==2) {
    for (int i=0; i<n; i++) {
      double sum=0;
      
      if (i>1) sum+=pts[i-2];
      else if (i>0) sum+=pts[i-1];
      else sum+=pts[i];
      
      if (i>0) sum+=pts[i-1];
      else sum+=pts[i];
      
      sum+=pts[i];
      
      if (i<n-1) sum+=pts[i+1];
      else sum+=pts[i];

      
      if (i<n-2) sum+=pts[i+2];
      else if (i<n-1) sum+=pts[i+1];
      else sum+=pts[i];
      
      sum/=5;
      ptsOut[i]=sum;  
    }
  }
  if ((mode/10)==1) {
    int window=mode%10;


    float x[100],y[100];
    for (int i=0; i<=2*window; i++) x[i]=i-window;

    TF1 f1("f1","pol1",-window+0.5,window-0.5);

    for (int i=0; i<n; i++) {

      for (int j=-window; j<=window; j++) {
	if ((i+j)<0) y[j+window]=pts[0];
	else if ((i+j)>=n) y[j+window]=pts[n-1];
	else y[j+window]=pts[i+j];
      }

      TGraph tg(window*2+1,x,y);

      tg.Fit(&f1,"QR","",-window+0.5,+window-0.5);
      
      ptsOut[i]=f1.Eval(0);
    }
  }


  for (int i=0; i<n; i++) pts[i]=ptsOut[i];
}

#include "tdrstyle.C"

void limWR(const char* fname,int which, int smooth=0,const char* asUsed=0,double lumi=5.0) {
  float mw[100], mn[100], mwt[100];
  float obs[100], exp[100], expm2s[100],expm1s[100],expp1s[100],expp2s[100],xsec[100];
  float obs_ns[100],exp_ns[100];
  int n=0,nx=0;
  char buffer[1024];
  double pbr=1e3; // pb ratio
  double corr=1.0;

  setTDRStyle();

  int ecm;
  if (which==0) { // 2011 Muons
    ecm=7;
  } else if (which==1) { // 2012 Muons
    ecm=8;
  } else if (which==2) { // 2012 Electrons
    ecm=8;
  } else if (which==3) { // 2012 Muons+Electrons
    ecm=8;
  } else if (which==4) { // USe the 8 TeV xsec as master
    ecm=8;
  }


  bool correctToMWMN2=true;

  FILE* f1=fopen(fname,"rt");
  int dummyi;

  
  while (!feof(f1)) {
    buffer[0]=0;
    fgets(buffer,1000,f1);
    if (sscanf(buffer,"%f %f %f %f %f %f %f %f ",mw+n,mn+n,
	       obs+n,exp+n,expm2s+n,expm1s+n,expp1s+n,expp2s+n)==8) {
      if (correctToMWMN2) {
	corr=1.0/accRatio(mw[n]);
	mn[n]=mw[n]/2.0;
      }

      if (which==3) {
	corr=3.0/2.0; // correct to all leptons
      }


      obs[n]*=pbr*corr;
      exp[n]*=pbr*corr;
      expm2s[n]*=pbr*corr;
      expm1s[n]*=pbr*corr;
      expp1s[n]*=pbr*corr;
      expp2s[n]*=pbr*corr;
      n++;
    }
  }
  fclose(f1);


  FILE* fx=fopen("Limits/plot/cs.txt","rt");
  float rawx, kf,xmw,xmn;
  int ecmx;
  while (!feof(fx)) {
    buffer[0]=0;
    fgets(buffer,1000,fx);
    if (sscanf(buffer,"%d %f %f %f %f",&ecmx, &xmw,&xmn,&rawx,&kf)==5) {
      if (ecmx!=ecm) continue;
      for (int search=0; search<n; search++) {
	if (fabs(xmw-mw[search])<1.0 && fabs(xmn-mn[search])<1.0) {
	  xsec[nx]=rawx*kf*pbr;
	  mwt[nx]=xmw;

	  if (which==3) {
	    xsec[nx]*=8.0/3.0; // from Pythia
	  }

	  //      xsec[nx]/=0.75;
	  nx++;
	  break;
	}
      }
    }
  }
  fclose(fx);


  if (which==4) {

    for (int i=0; i<n; i++) {

      double ax=1;
      for (int search=0; search<nx; search++) {
	if (fabs(mw[i]-mwt[search])<1.0) {
	  ax=xsec[search];
	  xsec[search]=1.0;
	  break;
	}
      }
  
      obs[i]/=ax;
      exp[i]/=ax;
      expm2s[i]/=ax;
      expm1s[i]/=ax;
      expp1s[i]/=ax;
      expp2s[i]/=ax;
      
    }

  }
  
  //  gROOT->SetStyle("Plain");
  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetTopMargin(0.04);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetTicks(1,1);
  c1->SetLogy();

  //  TH1* dummy=new TH1F("dummy","",30,1000,3000);
  TH1* dummy;
  if (which!=0) dummy=new TH1F("dummy","",30,1000,3000);
  else dummy=new TH1F("dummy","",30,1000,2600);

  //  dummy->SetMinimum(0.5);
  dummy->SetMinimum(0.5);
  if (which==4) {
   dummy->SetMaximum(10);
   dummy->SetMinimum(0.005);
  } else if (which==3) {
    dummy->SetMaximum(1000);
  } else {
    dummy->SetMaximum(500);
  }
  dummy->SetStats(0);
  dummy->GetXaxis()->SetNdivisions(507);
  dummy->GetXaxis()->SetTitle("M_{W_{R}} [GeV]");
  
  if (which==0 || which==1) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R})#times BR(W_{R}#rightarrow #mu#mu jj) [fb]");
  } else if (which==2) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R})#times BR(W_{R}#rightarrow ee jj) [fb]");
  } else if (which==3) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R})#times BR(W_{R}#rightarrow (ee+#mu#mu+#tau#tau) jj) [fb]");
    dummy->GetYaxis()->SetTitleSize(0.055)
  } else if (which==4) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R})#times BR(W_{R}#rightarrow #mu#mu jj)/#sigma_{g_{R}=g_{L}}");
  }
  dummy->GetYaxis()->SetTitleOffset(1.1);
  dummy->Draw("HIST");

  TLegend* tl;
  if (which==0 || which==1) {
    tl=new TLegend(0.50,0.93,0.95,0.60);//,"CL_{S} Method          M_{N_{#mu}}=M_{W_{R}}/2");
    tl->SetTextFont(42);
    tl->SetHeader("CL_{S} Method          M_{N_{#mu}}=M_{W_{R}}/2");
  } else if (which==3) {
    tl=new TLegend(0.50,0.93,0.95,0.62,"CL_{S} Method          M_{N}=M_{W_{R}}/2");
  } else if (which==2) {
    tl=new TLegend(0.50,0.93,0.95,0.60,"CL_{S} Method          M_{N_{e}}=M_{W_{R}}/2");
  } else if (which==4) {
    tl=new TLegend(0.16,0.93,0.90,0.75,"CL_{S} Method          M_{N_{#mu}}=M_{W_{R}}/2");
    tl->SetNColumns(2);
  }
  tl->SetTextFont(42);
  tl->SetFillStyle(0);


  if ((smooth/100)%10) {
    smoothSG(expm2s,n);
    smoothSG(expp2s,n);
  }

  if ((smooth/10)%10) {
    smoothSG(expm1s,n);
    smoothSG(expp1s,n);
  }

  if ((smooth/1)%10) {
    smoothSG(exp,n);
  }

  if ((smooth/1000)%10) {
    smoothSG(obs,n,(smooth/1000)%10);
  }

  float mw2[200],e2s[200],e1s[200];
  for (int i=0; i<n; i++) {
    mw2[i]=mw[i];
    mw2[n+i]=mw[n-i-1];
    e2s[i]=expm2s[i];
    e2s[n+i]=expp2s[n-i-1];
    e1s[i]=expm1s[i];
    e1s[n+i]=expp1s[n-i-1];
    //    xsec[i]*=0.75; // 0.75 is average reduced xsection for a single channel * k factor (determined at MW=1.3 TeV)
  }

  TGraph* tg_e2s=new TGraph(n*2,mw2,e2s);
  tg_e2s->SetLineWidth(0);
  tg_e2s->SetFillColor(kYellow);
  tg_e2s->Draw("F SAME");

  TGraph* tg_e1s=new TGraph(n*2,mw2,e1s);
  tg_e1s->SetLineWidth(0);
  tg_e1s->SetFillColor(kGreen);
  tg_e1s->Draw("F SAME");

  TGraph* tg_exp=new TGraph(n,mw,exp);
  tg_exp->SetLineWidth(2);
  tg_exp->SetLineColor(kBlue);
  tg_exp->SetLineStyle(2);
  tg_exp->Draw("L SAME");


  TGraph* tg_theory=new TGraph(nx,mwt,xsec);
  tg_theory->SetLineWidth(2);
  tg_theory->SetLineColor(kRed);
  tg_theory->Draw("L SAME");

  TGraph* tg_obs=new TGraph(n,mw,obs);
  tg_obs->SetLineWidth(2);
  tg_obs->Draw("L SAME");

  tl->AddEntry(tg_obs,"Observed Limit (95%CL)","L");
  tl->AddEntry(tg_exp,"Expected Limit (95%CL)","L");
  tl->AddEntry(tg_e1s,"Expected #pm 1 #sigma","F");
  tl->AddEntry(tg_e2s,"Expected #pm 2 #sigma","F");


  if (which==3) {
    tl->AddEntry(tg_theory,"Theory Expectation","L");
    tl->AddEntry((TObject*)0,"(g_{R}=g_{L}, M_{N}_{e}=M_{N}_{#mu}=M_{N}_{#tau})","");
  } else {
    tl->AddEntry(tg_theory,"Theory Expectation (g_{R}=g_{L})","L");
  }

  TText *text;
  if (which==0) {
    
    text = new TText(0.89,0.97,"CMS");
  } else {
    text = new TText(0.72,0.97,"CMS Preliminary");
  }
  //TText *text = new TText(0.89,0.94,"CMS");
   text->SetNDC();
   text->SetTextFont(42);
   //   text->SetTextSize(0.03);
   text->SetTextSize(0.03);
   text->Draw();

   /*
   TDatime now;
   sprintf(buffer,"%d/%d/%d",now.GetDay(),now.GetMonth(),now.GetYear());

   text = new TText(0.80,0.88,buffer);
   text->SetNDC();
   text->SetTextSize(0.03);
   text->Draw();
   */
   /*
   if (which==3) {
     text = new TText(0.16,0.16,"No systematics");
     text->SetNDC();
     text->SetTextColor(kBlue);
     text->SetTextSize(0.04);
     text->Draw();
   }
   */

   

   
   if (which==4) {
     TLatex *   texa = new TLatex(0.152,0.972,"#int");
     texa->SetTextFont(42);
     texa->SetNDC();
     texa->SetTextSize(0.017);
     texa->SetLineWidth(2);
     texa->Draw();
     sprintf(buffer,"L dt = %.1f fb^{-1} at #sqrt{s}=%d TeV +",lumi,ecm);
     TLatex *   tex = new TLatex(0.162,0.970,buffer);
     tex->SetNDC();
     tex->SetTextSize(0.026);
     tex->SetLineWidth(2);
     tex->Draw();

     sprintf(buffer,"%.1f fb^{-1} at #sqrt{s}=%d TeV",5.0,7);
     tex = new TLatex(0.47,0.970,buffer);
     tex->SetNDC();
     tex->SetTextSize(0.026);
     tex->SetLineWidth(2);
     tex->Draw();
   } else {
     TLatex *   texa = new TLatex(0.152,0.972,"#int");
     texa->SetNDC();
     texa->SetTextFont(42);
     texa->SetTextSize(0.017);
     texa->SetLineWidth(2);
     texa->Draw();
     sprintf(buffer,"L dt = %.1f fb^{-1} at #sqrt{s}=%d TeV",lumi,ecm);
     TLatex *   tex = new TLatex(0.162,0.97,buffer);
     tex->SetNDC();
     tex->SetTextFont(42);
     tex->SetTextSize(0.025);
     tex->SetLineWidth(2);
     tex->Draw();
   }

   /*

   tex = new TLatex(0.75,0.5,"M(N_{R})=#frac{1}{2} M(W_{R})");
   tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->Draw();
   */

     
   tl->Draw("SAME");
  c1->RedrawAxis();
  char feps[100];
  char fpng[100];
  char* fbase="limWR";

  if (which==4) fbase="limWR_2y";
  if (which==3) fbase="limWR_emu";
  if (which==2) fbase="limWR_e12";
  if (which==1) fbase="limWR_mu12";
  if (which==0) fbase="limWR_mu11";

  sprintf(feps,"%s.eps",fbase);
  sprintf(fpng,"%s.png",fbase);
  c1->Print(feps);
  c1->Print(fpng);
  char cmd[100];
  sprintf(cmd,"eps2pdf %s",feps);
  system(cmd);

  if (asUsed!=0) {
    FILE* f=fopen(asUsed,"w");

    for (int i=0; i<n; i++) {
      fprintf(f,"%4d %4d %7f %7f %7f %7f %7f %7f 0 \n",
	      int(tg_obs->GetX()[i]),int(tg_obs->GetX()[i])/2,
	      tg_obs->GetY()[i],tg_exp->GetY()[i],
	      expm2s[i],expm1s[i],expp1s[i],expp2s[i]);

    }
    fclose(f);
  }

}
