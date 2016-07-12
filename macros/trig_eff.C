void trig_eff()
{
	//=========Macro generated from canvas: c/c
	////=========  (Mon Jul 11 18:44:14 2016) by ROOT version6.02/05
	TCanvas *c = new TCanvas("c", "c",1,1,800,776);
	gStyle->SetOptStat(0);
	c->SetHighLightColor(2);
	c->Range(-992.3078,0.9387654,7469.231,1.025185);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetBorderSize(2);
	c->SetLeftMargin(0.2);
	c->SetRightMargin(0.15);
	c->SetTopMargin(0.06);
	c->SetBottomMargin(0.13);
	c->SetFrameFillStyle(0);
	c->SetFrameBorderMode(0);
	c->SetFrameFillStyle(0);
	c->SetFrameBorderMode(0);

	TH1F *h3 = new TH1F("h3","",1,700,6200);
	h3->SetMinimum(0.95);
	h3->SetMaximum(1.02);
	h3->SetStats(0);
	h3->SetLineWidth(3);
	h3->SetMarkerStyle(20);
	h3->GetXaxis()->SetTitle("W_{R} Mass Hypothesis [GeV]");
	h3->GetXaxis()->SetLabelFont(42);
	h3->GetXaxis()->SetLabelOffset(0.007);
	h3->GetXaxis()->SetLabelSize(0.05);
	h3->GetXaxis()->SetTitleSize(0.05);
	h3->GetXaxis()->SetTitleFont(42);
	h3->GetYaxis()->SetTitle("Trigger Efficiency | Offline Selection");
	h3->GetYaxis()->SetLabelFont(42);
	h3->GetYaxis()->SetLabelOffset(0.007);
	h3->GetYaxis()->SetLabelSize(0.05);
	h3->GetYaxis()->SetTitleSize(0.05);
	h3->GetYaxis()->SetTitleOffset(1.5);
	h3->GetYaxis()->SetTitleFont(42);
	h3->GetZaxis()->SetLabelFont(42);
	h3->GetZaxis()->SetLabelOffset(0.007);
	h3->GetZaxis()->SetLabelSize(0.05);
	h3->GetZaxis()->SetTitleSize(0.05);
	h3->GetZaxis()->SetTitleFont(42);
	h3->Draw("");

	Double_t Graph0_fx9[25] = {
		800,
		1000,
		1200,
		1400,
		1600,
		1800,
		2000,
		2200,
		2400,
		2600,
		2800,
		3000,
		3200,
		3600,
		3800,
		4000,
		4200,
		4400,
		4600,
		4800,
		5000,
		5200,
		5600,
		5800,
		6000};
	Double_t Graph0_fy9[25] = {
		0.9828543,
		0.9838024,
		0.986543,
		0.9876861,
		0.9886193,
		0.9892568,
		0.990447,
		0.9909273,
		0.9912287,
		0.9907072,
		0.9913204,
		0.9921334,
		0.99108,
		0.9914352,
		0.9921561,
		0.9918741,
		0.9926426,
		0.992206,
		0.9915519,
		0.9925439,
		0.9928637,
		0.9928227,
		0.9922062,
		0.99233,
		0.9916502};
	TGraph *graph = new TGraph(25,Graph0_fx9,Graph0_fy9);
	graph->SetName("Graph0");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	Int_t ci;      // for color index setting
	TColor *color; // for color definition with alpha
	ci = TColor::GetColor("#00ff00");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","Graph",100,280,6520);
	Graph_Graph9->SetMinimum(0.9818534);
	Graph_Graph9->SetMaximum(0.9938646);
	Graph_Graph9->SetDirectory(0);
	Graph_Graph9->SetStats(0);
	Graph_Graph9->SetLineWidth(3);
	Graph_Graph9->SetMarkerStyle(20);
	Graph_Graph9->GetXaxis()->SetLabelFont(42);
	Graph_Graph9->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph9->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph9->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph9->GetXaxis()->SetTitleFont(42);
	Graph_Graph9->GetYaxis()->SetLabelFont(42);
	Graph_Graph9->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph9->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph9->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph9->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph9->GetYaxis()->SetTitleFont(42);
	Graph_Graph9->GetZaxis()->SetLabelFont(42);
	Graph_Graph9->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph9->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph9->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph9->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph9);

	graph->Draw("l");

	Double_t Graph1_fx10[25] = {
		800,
		1000,
		1200,
		1400,
		1600,
		1800,
		2000,
		2200,
		2400,
		2600,
		2800,
		3000,
		3200,
		3600,
		3800,
		4000,
		4200,
		4400,
		4600,
		4800,
		5000,
		5200,
		5600,
		5800,
		6000};
	Double_t Graph1_fy10[25] = {
		0.9968987,
		0.9971195,
		0.9960687,
		0.9966613,
		0.9965108,
		0.9965275,
		0.9955773,
		0.9962382,
		0.99635,
		0.9953004,
		0.9948061,
		0.9945437,
		0.9950997,
		0.9941074,
		0.9947418,
		0.9944833,
		0.9941475,
		0.9940306,
		0.9932445,
		0.9940471,
		0.9934309,
		0.993231,
		0.9919418,
		0.9924854,
		0.9917312};
	graph = new TGraph(25,Graph1_fx10,Graph1_fy10);
	graph->SetName("Graph1");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	ci = TColor::GetColor("#ff0000");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","Graph",100,280,6520);
	Graph_Graph10->SetMinimum(0.9911924);
	Graph_Graph10->SetMaximum(0.9976584);
	Graph_Graph10->SetDirectory(0);
	Graph_Graph10->SetStats(0);
	Graph_Graph10->SetLineWidth(3);
	Graph_Graph10->SetMarkerStyle(20);
	Graph_Graph10->GetXaxis()->SetLabelFont(42);
	Graph_Graph10->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph10->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph10->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph10->GetXaxis()->SetTitleFont(42);
	Graph_Graph10->GetYaxis()->SetLabelFont(42);
	Graph_Graph10->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph10->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph10->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph10->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph10->GetYaxis()->SetTitleFont(42);
	Graph_Graph10->GetZaxis()->SetLabelFont(42);
	Graph_Graph10->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph10->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph10->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph10->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph10);

	graph->Draw("l");

	TLegend *leg = new TLegend(0.3,0.8,0.85,0.9,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.032);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	TLegendEntry *entry=leg->AddEntry("Graph0","Electron Channel","L");

	ci = TColor::GetColor("#00ff00");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph1","Muon Channel","L");

	ci = TColor::GetColor("#ff0000");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	leg->Draw();
	c->Modified();
	c->cd();
	c->SetSelected(c);
}
