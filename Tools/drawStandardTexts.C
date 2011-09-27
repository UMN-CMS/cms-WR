#include <iostream>
#include <string>
#include "TText.h"
#include "TLatex.h"
#include "TDatime.h"

using namespace std;
static bool dotime=false;

void drawStandardText(const std::string& text,
		      double px,    double py,
		      double tx=-1, double ty=-1,
		      double textsize=0.035)
{
  const char* time_;
  TDatime mytime;
  time_ = mytime.AsString();

  TText *plabel = new TText();
  plabel-> SetNDC();
  plabel -> SetTextFont(522);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(textsize);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);

  
  //Then for each plot, pick a nice spot and draw
  plabel -> DrawText(px, py, text.c_str());
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

void drawLumiText(double px, double py, float val) {
  char text[40];
  sprintf(text,"%.1f pb^{-1}",val);

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
