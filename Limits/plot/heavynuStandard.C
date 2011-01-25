#include "TText.h"
#include "TLatex.h"
#include "TDatime.h"

static bool dotime=true;

void label_Prelim(double px, double py, double tx=-1, double ty=-1) {
  const char* time_;
  TDatime mytime;
  time_ = mytime.AsString();

  TText *plabel = new TText();
  plabel-> SetNDC();
  plabel -> SetTextFont(42);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.035);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);

  
  //Then for each plot, pick a nice spot and draw
  plabel -> DrawText(px, py, "CMS 2010 PRELIMINARY");
  if (tx>=0 && ty>=0 && dotime) {
    TText *tlabel = new TText();
    tlabel-> SetNDC();
    tlabel -> SetTextFont(42);
    tlabel -> SetTextColor(1);
    tlabel -> SetTextSize(0.02);
    tlabel -> SetTextAlign(22);
    tlabel -> SetTextAngle(0);
    
    tlabel -> DrawText(tx, ty, Form("%s",time_));
  }
}

void label_Lumi(double px, double py, int val) {
  char text[40];
  sprintf(text,"%d pb^{-1}",val);

  TLatex *plabel = new TLatex(px,py,text);
  plabel-> SetNDC();
  plabel -> SetTextFont(42);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.040);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);


  //Then for each plot, pick a nice spot and draw
  plabel -> Draw();

}
