double findEstimate(TTree* tuple, double target, double llimit=0.05) {

  TGraph* tg;

  /*
  double* x=new double[tuple->GetSelectedRows()];
  double* y=new double[tuple->GetSelectedRows()];
  */

  // get the data in, and sort it.

  // next, remove points which are very far off being continous.

  double* ey=new double[tuple->GetSelectedRows()];
  for (int i=0; i<tuple->GetSelectedRows(); i++) {
    ey[i]=0.05*tuple->GetV1()[i];
    if (i!=0 && i!=tuple->GetSelectedRows()-1) {
      double nave=(tuple->GetV1()[i-1]+tuple->GetV1()[i+1])/2.0;
      if (tuple->GetV1()[i]<0.05 && tuple->GetV1()[i]/nave<0.5)
	ey[i]=1.0;
    }
    if (ey[i]<0.01) ey[i]=0.01;
  }
 

  tg=new TGraphErrors(tuple->GetSelectedRows(),tuple->GetV2(),tuple->GetV1(),0,ey);

  TF1* f1=new TF1("exp","exp([1]*x+[0])",0.01,2.5);

  f1->SetParameter(0,1);
  f1->SetParLimits(0,-10,10.0);
  f1->SetParameter(1,-1);

  tg->Fit(f1,"Q","",llimit,1.5);
  tg->Draw("APL");

  return log(target)/f1->GetParameter(1)-f1->GetParameter(0)/f1->GetParameter(1);

  /*
  double near_low_x, near_low_y;
  double near_high_x, near_high_y;
  
  int n=tuple->GetSelectedRows();
  double* xvec=tuple->GetV1();
  double* yvec=tuple->GetV2();



  near_low_x=xvec[0]; near_high_x=xvec[0];
  near_low_y=yvec[0]; near_high_y=yvec[0];

  for (int i=1; i<n; i++) {
    if (fabs(xvec[i]-target)<fabs(near_low_x-target) && xvec[i]<target) {
      near_low_x=xvec[i];
      near_low_y=yvec[i];
    }

    if (fabs(xvec[i]-target)<fabs(near_high_x-target) && xvec[i]>target) {
      near_high_x=xvec[i];
      near_high_y=yvec[i];
    }
  }

  //  printf(" %f %f %f %f \n",near_low_x, near_low_y, near_high_x, near_high_y);

  // interpolate
  double esti=(target-near_low_x)*near_low_y+(near_high_x-target)*near_high_y;
  esti/=(near_high_x-near_low_x);
  */
  

  return esti;
}

double limit_exp, limit_exp_p1s, limit_exp_p2s, limit_exp_m1s, limit_exp_m2s, limit_obs;

void fitLimitValue(TFile* file, double mw, double mn, bool batch=false) {

  TTree* tr=file->Get("WRLimit");
  
  char cut[120];
  sprintf(cut,"TMath::Abs(info.mnu-%.0f)<2 && TMath::Abs(info.mwr-%0.f)<2",mn,mw);

  TCanvas *c1=0;
  if (!batch) c1=new TCanvas("c1","c1",800,600);

  if (c1!=0) c1->Divide(3,2);
  TGraph* tg;

  if (c1!=0) c1->cd(1);
  tr->Draw("info.cl_sb_exp:info.xsec",cut,"BOX");
  limit_exp=findEstimate(tr,0.05);
  /*
  tg=new TGraph(tr->GetSelectedRows(),tr->GetV2(),tr->GetV1());
  tg->Draw("AL");
  */

  if (c1!=0) c1->cd(2);
  tr->Draw("info.cl_sb_exp_p1s:info.xsec",cut,"BOX");
  limit_exp_p1s=findEstimate(tr,0.05);
  /*
  tg=new TGraph(tr->GetSelectedRows(),tr->GetV2(),tr->GetV1());
  tg->Draw("AL");
  */

  if (c1!=0) c1->cd(3);
  tr->Draw("info.cl_sb_exp_p2s:info.xsec",cut,"BOX");
  limit_exp_p2s=findEstimate(tr,0.05);
  /*
  tg=new TGraph(tr->GetSelectedRows(),tr->GetV2(),tr->GetV1());
  tg->Draw("AL");
  */

  if (c1!=0) c1->cd(4);
  tr->Draw("info.cl_sb_exp_m1s:info.xsec",cut,"BOX");
  limit_exp_m1s=findEstimate(tr,0.05,0.01);
  if (limit_exp_m1s<0) limit_exp_m1s=0;

  if (c1!=0) c1->cd(5);
  tr->Draw("info.cl_sb_exp_m2s:info.xsec",cut,"BOX");
  limit_exp_m2s=findEstimate(tr,0.05,0.01);
  if (limit_exp_m2s<0) limit_exp_m2s=0;
  
  if (c1!=0) c1->cd(6);
  tr->Draw("info.cl_sb_obs:info.xsec",cut,"BOX");
  limit_obs=findEstimate(tr,0.05);

  /*
  tg->Draw("AL");
  */


  printf("Limits for %.0f/%0.f : %.3f pb expected, %.3f pb observed\n", mw, mn, limit_exp, limit_obs);
  printf("  ranges +2 %.3f +1 %.3f -1 %.3f -2 %.3f \n",limit_exp_p2s,limit_exp_p1s,limit_exp_m1s,limit_exp_m2s);



}
void fitLimitValues(TFile* file, const char* ofname=0) {
  int pairs[100];
  int np=0;

  TTree* tr=file->Get("WRLimit");
  
  tr->Draw("info.mwr:info.mnu");

  int nv=tr->GetSelectedRows();
  
  FILE* plotf=0;
  if (ofname!=0) plotf=fopen(ofname,"w");

  int j;
  for (int i=0; i<nv; i++) {
    int code=int(tr->GetV1()[i])*10000+int(tr->GetV2()[i]);
    if (code<10000) continue;
    for (j=0; j<np; j++)
      if (pairs[j]==code) break;
    if (j==np) {
      pairs[np]=code;
      np++;

      printf("%3d ",np);
      double mw=tr->GetV1()[i];
      double mnu=tr->GetV2()[i];

      plotFit(file,mw,mnu,true);

      if (plotf!=0) 
	fprintf(plotf,"%.0lf,%.0lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf\n",mw,mnu,
		limit_exp,limit_exp_p2s,limit_exp_p1s,limit_exp_m1s,limit_exp_m2s,limit_obs);
    }
  }

  if (plotf!=0) fclose(plotf);
}
