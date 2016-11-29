void mu_eff()
{
	//=========Macro generated from canvas: c/c
	//=========  (Thu Jun 30 07:52:05 2016) by ROOT version6.02/05
	TCanvas *c = new TCanvas("c", "c",1,1,800,776);
	gStyle->SetOptStat(0);
	c->SetHighLightColor(2);
	c->Range(-204.1095,-0.1604938,7330.137,1.074074);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetBorderSize(2);
	c->SetLeftMargin(0.12);
	c->SetRightMargin(0.15);
	c->SetTopMargin(0.06);
	c->SetBottomMargin(0.13);
	c->SetFrameFillStyle(0);
	c->SetFrameBorderMode(0);
	c->SetFrameFillStyle(0);
	c->SetFrameBorderMode(0);

	TH1F *h2 = new TH1F("h2","",1,700,6200);
	h2->SetMinimum(0);
	h2->SetMaximum(1);
	h2->SetStats(0);
	h2->SetLineWidth(3);
	h2->SetMarkerStyle(20);
	h2->GetXaxis()->SetTitle("W_{R} Mass [GeV]");
	h2->GetXaxis()->SetLabelFont(42);
	h2->GetXaxis()->SetLabelOffset(0.007);
	h2->GetXaxis()->SetLabelSize(0.05);
	h2->GetXaxis()->SetTitleSize(0.05);
	h2->GetXaxis()->SetTitleFont(42);
	h2->GetYaxis()->SetTitle("Full Selection Eff*Accept");
	h2->GetYaxis()->SetLabelFont(42);
	h2->GetYaxis()->SetLabelOffset(0.007);
	h2->GetYaxis()->SetLabelSize(0.05);
	h2->GetYaxis()->SetTitleSize(0.05);
	h2->GetYaxis()->SetTitleOffset(1.1);
	h2->GetYaxis()->SetTitleFont(42);
	h2->GetZaxis()->SetLabelFont(42);
	h2->GetZaxis()->SetLabelOffset(0.007);
	h2->GetZaxis()->SetLabelSize(0.05);
	h2->GetZaxis()->SetTitleSize(0.05);
	h2->GetZaxis()->SetTitleFont(42);
	h2->Draw("");

	Double_t Graph0_fx5[25] = {
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
	Double_t Graph0_fy5[25] = {
		0.225654,
		0.29196,
		0.34404,
		0.39171,
		0.42122,
		0.45018,
		0.480687,
		0.49612,
		0.51334,
		0.52252,
		0.534379,
		0.54326,
		0.552249,
		0.564649,
		0.5684,
		0.57064,
		0.578565,
		0.57924,
		0.57816,
		0.57796,
		0.5802,
		0.58602,
		0.579003,
		0.58464,
		0.576451};
	TGraph *graph = new TGraph(25,Graph0_fx5,Graph0_fy5);
	graph->SetName("Graph0");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	Int_t ci;      // for color index setting
	TColor *color; // for color definition with alpha
	ci = TColor::GetColor("#00ff00");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,280,6520);
	Graph_Graph5->SetMinimum(0.1896174);
	Graph_Graph5->SetMaximum(0.6220566);
	Graph_Graph5->SetDirectory(0);
	Graph_Graph5->SetStats(0);
	Graph_Graph5->SetLineWidth(3);
	Graph_Graph5->SetMarkerStyle(20);
	Graph_Graph5->GetXaxis()->SetLabelFont(42);
	Graph_Graph5->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph5->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph5->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph5->GetXaxis()->SetTitleFont(42);
	Graph_Graph5->GetYaxis()->SetLabelFont(42);
	Graph_Graph5->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph5->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph5->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph5->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph5->GetYaxis()->SetTitleFont(42);
	Graph_Graph5->GetZaxis()->SetLabelFont(42);
	Graph_Graph5->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph5->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph5->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph5->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph5);

	graph->Draw("l");

	Double_t Graph1_fx6[25] = {
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
	Double_t Graph1_fy6[25] = {
		0.119425,
		0.14566,
		0.15378,
		0.15813,
		0.16122,
		0.15806,
		0.15029,
		0.14904,
		0.14296,
		0.1438,
		0.135777,
		0.1345,
		0.130784,
		0.126297,
		0.1241,
		0.123182,
		0.119963,
		0.1215,
		0.11816,
		0.11684,
		0.11966,
		0.11322,
		0.117102,
		0.11198,
		0.116919};
	graph = new TGraph(25,Graph1_fx6,Graph1_fy6);
	graph->SetName("Graph1");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	ci = TColor::GetColor("#ff0000");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,280,6520);
	Graph_Graph6->SetMinimum(0.107056);
	Graph_Graph6->SetMaximum(0.166144);
	Graph_Graph6->SetDirectory(0);
	Graph_Graph6->SetStats(0);
	Graph_Graph6->SetLineWidth(3);
	Graph_Graph6->SetMarkerStyle(20);
	Graph_Graph6->GetXaxis()->SetLabelFont(42);
	Graph_Graph6->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph6->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph6->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph6->GetXaxis()->SetTitleFont(42);
	Graph_Graph6->GetYaxis()->SetLabelFont(42);
	Graph_Graph6->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph6->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph6->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph6->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph6->GetYaxis()->SetTitleFont(42);
	Graph_Graph6->GetZaxis()->SetLabelFont(42);
	Graph_Graph6->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph6->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph6->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph6->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph6);

	graph->Draw("l");

	Double_t Graph2_fx7[25] = {
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
	Double_t Graph2_fy7[25] = {
		0.0131742,
		0.01276,
		0.0119,
		0.00944123,
		0.00996,
		0.00912,
		0.00805946,
		0.00762,
		0.00788,
		0.00752,
		0.00823267,
		0.00782,
		0.00778453,
		0.00768132,
		0.00754,
		0.00769621,
		0.00801892,
		0.00708,
		0.00768,
		0.0076,
		0.00738,
		0.00788,
		0.00740018,
		0.00752,
		0.00820628};
	graph = new TGraph(25,Graph2_fx7,Graph2_fy7);
	graph->SetName("Graph2");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	ci = TColor::GetColor("#0000ff");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,280,6520);
	Graph_Graph7->SetMinimum(0.00647058);
	Graph_Graph7->SetMaximum(0.01378362);
	Graph_Graph7->SetDirectory(0);
	Graph_Graph7->SetStats(0);
	Graph_Graph7->SetLineWidth(3);
	Graph_Graph7->SetMarkerStyle(20);
	Graph_Graph7->GetXaxis()->SetLabelFont(42);
	Graph_Graph7->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph7->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph7->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph7->GetXaxis()->SetTitleFont(42);
	Graph_Graph7->GetYaxis()->SetLabelFont(42);
	Graph_Graph7->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph7->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph7->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph7->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph7->GetYaxis()->SetTitleFont(42);
	Graph_Graph7->GetZaxis()->SetLabelFont(42);
	Graph_Graph7->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph7->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph7->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph7->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph7);

	graph->Draw("l");

	Double_t Graph3_fx8[25] = {
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
	Double_t Graph3_fy8[25] = {
		0.399554,
		0.49848,
		0.56248,
		0.616739,
		0.65116,
		0.67726,
		0.695807,
		0.70974,
		0.72064,
		0.72854,
		0.732811,
		0.74004,
		0.744447,
		0.750906,
		0.75294,
		0.754778,
		0.759472,
		0.75934,
		0.75572,
		0.75478,
		0.75614,
		0.75714,
		0.753354,
		0.75546,
		0.752725};
	graph = new TGraph(25,Graph3_fx8,Graph3_fy8);
	graph->SetName("Graph3");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Graph",100,280,6520);
	Graph_Graph8->SetMinimum(0.3635622);
	Graph_Graph8->SetMaximum(0.7954638);
	Graph_Graph8->SetDirectory(0);
	Graph_Graph8->SetStats(0);
	Graph_Graph8->SetLineWidth(3);
	Graph_Graph8->SetMarkerStyle(20);
	Graph_Graph8->GetXaxis()->SetLabelFont(42);
	Graph_Graph8->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph8->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph8->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph8->GetXaxis()->SetTitleFont(42);
	Graph_Graph8->GetYaxis()->SetLabelFont(42);
	Graph_Graph8->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph8->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph8->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph8->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph8->GetYaxis()->SetTitleFont(42);
	Graph_Graph8->GetZaxis()->SetLabelFont(42);
	Graph_Graph8->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph8->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph8->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph8->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph8);

	graph->Draw("l");

	TLegend *leg = new TLegend(0.3,0.8,0.85,0.9,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.032);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	TLegendEntry *entry=leg->AddEntry("Graph0","Muon BB","L");

	ci = TColor::GetColor("#00ff00");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph1","Muon EB","L");

	ci = TColor::GetColor("#ff0000");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph2","Muon EE","L");

	ci = TColor::GetColor("#0000ff");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph3","Muon Global","L");
	entry->SetLineColor(1);
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
