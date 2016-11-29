void ele_eff()
{
	//=========Macro generated from canvas: c/c
	//=========  (Thu Jun 30 07:52:05 2016) by ROOT version6.02/05
	TCanvas *c = new TCanvas("c", "c",0,0,800,800);
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

	TH1F *h1 = new TH1F("h1","",1,700,6200);
	h1->SetMinimum(0);
	h1->SetMaximum(1);
	h1->SetStats(0);
	h1->SetLineWidth(3);
	h1->SetMarkerStyle(20);
	h1->GetXaxis()->SetTitle("W_{R} Mass [GeV]");
	h1->GetXaxis()->SetLabelFont(42);
	h1->GetXaxis()->SetLabelOffset(0.007);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitleFont(42);
	h1->GetYaxis()->SetTitle("Full Selection Eff*Accept");
	h1->GetYaxis()->SetLabelFont(42);
	h1->GetYaxis()->SetLabelOffset(0.007);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.1);
	h1->GetYaxis()->SetTitleFont(42);
	h1->GetZaxis()->SetLabelFont(42);
	h1->GetZaxis()->SetLabelOffset(0.007);
	h1->GetZaxis()->SetLabelSize(0.05);
	h1->GetZaxis()->SetTitleSize(0.05);
	h1->GetZaxis()->SetTitleFont(42);
	h1->Draw("");

	Double_t Graph0_fx1[25] = {
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
	Double_t Graph0_fy1[25] = {
		0.18164,
		0.24126,
		0.287457,
		0.32276,
		0.352793,
		0.375318,
		0.396509,
		0.412382,
		0.42528,
		0.4351,
		0.44626,
		0.461,
		0.459311,
		0.471313,
		0.472735,
		0.47522,
		0.48294,
		0.478787,
		0.481565,
		0.48618,
		0.483493,
		0.48192,
		0.482911,
		0.48134,
		0.485052};
	TGraph *graph = new TGraph(25,Graph0_fx1,Graph0_fy1);
	graph->SetName("Graph0");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	Int_t ci;      // for color index setting
	TColor *color; // for color definition with alpha
	ci = TColor::GetColor("#00ff00");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,280,6520);
	Graph_Graph1->SetMinimum(0.151186);
	Graph_Graph1->SetMaximum(0.516634);
	Graph_Graph1->SetDirectory(0);
	Graph_Graph1->SetStats(0);
	Graph_Graph1->SetLineWidth(3);
	Graph_Graph1->SetMarkerStyle(20);
	Graph_Graph1->GetXaxis()->SetLabelFont(42);
	Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph1->GetXaxis()->SetTitleFont(42);
	Graph_Graph1->GetYaxis()->SetLabelFont(42);
	Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph1->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph1->GetYaxis()->SetTitleFont(42);
	Graph_Graph1->GetZaxis()->SetLabelFont(42);
	Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph1->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph1->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph1);

	graph->Draw("l");

	Double_t Graph1_fx2[25] = {
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
	Double_t Graph1_fy2[25] = {
		0.09058,
		0.11086,
		0.118399,
		0.1228,
		0.121764,
		0.119821,
		0.118706,
		0.116376,
		0.11386,
		0.11124,
		0.10816,
		0.1029,
		0.102323,
		0.0988842,
		0.100177,
		0.09782,
		0.09488,
		0.093943,
		0.0951683,
		0.09394,
		0.092428,
		0.09318,
		0.0909613,
		0.0920536,
		0.0903946};
	graph = new TGraph(25,Graph1_fx2,Graph1_fy2);
	graph->SetName("Graph1");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	ci = TColor::GetColor("#ff0000");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,280,6520);
	Graph_Graph2->SetMinimum(0.08715406);
	Graph_Graph2->SetMaximum(0.1260405);
	Graph_Graph2->SetDirectory(0);
	Graph_Graph2->SetStats(0);
	Graph_Graph2->SetLineWidth(3);
	Graph_Graph2->SetMarkerStyle(20);
	Graph_Graph2->GetXaxis()->SetLabelFont(42);
	Graph_Graph2->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph2->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph2->GetXaxis()->SetTitleFont(42);
	Graph_Graph2->GetYaxis()->SetLabelFont(42);
	Graph_Graph2->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph2->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph2->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph2->GetYaxis()->SetTitleFont(42);
	Graph_Graph2->GetZaxis()->SetLabelFont(42);
	Graph_Graph2->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph2->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph2->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph2->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph2);

	graph->Draw("l");

	Double_t Graph2_fx3[25] = {
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
	Double_t Graph2_fy3[25] = {
		0.00752,
		0.00866,
		0.00783646,
		0.00682,
		0.00643218,
		0.00597167,
		0.00639563,
		0.00587799,
		0.00556,
		0.0059,
		0.00522,
		0.00616,
		0.00512236,
		0.00523433,
		0.00577034,
		0.00554,
		0.00502,
		0.00599376,
		0.00542212,
		0.0056,
		0.00591121,
		0.00588,
		0.00601075,
		0.00570226,
		0.00542411};
	graph = new TGraph(25,Graph2_fx3,Graph2_fy3);
	graph->SetName("Graph2");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);

	ci = TColor::GetColor("#0000ff");
	graph->SetLineColor(ci);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,280,6520);
	Graph_Graph3->SetMinimum(0.004656);
	Graph_Graph3->SetMaximum(0.009024);
	Graph_Graph3->SetDirectory(0);
	Graph_Graph3->SetStats(0);
	Graph_Graph3->SetLineWidth(3);
	Graph_Graph3->SetMarkerStyle(20);
	Graph_Graph3->GetXaxis()->SetLabelFont(42);
	Graph_Graph3->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph3->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph3->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph3->GetXaxis()->SetTitleFont(42);
	Graph_Graph3->GetYaxis()->SetLabelFont(42);
	Graph_Graph3->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph3->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph3->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph3->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph3->GetYaxis()->SetTitleFont(42);
	Graph_Graph3->GetZaxis()->SetLabelFont(42);
	Graph_Graph3->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph3->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph3->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph3->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph3);

	graph->Draw("l");

	Double_t Graph3_fx4[25] = {
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
	Double_t Graph3_fy4[25] = {
		0.27974,
		0.36078,
		0.413693,
		0.45238,
		0.480988,
		0.501111,
		0.521611,
		0.534636,
		0.5447,
		0.55224,
		0.55964,
		0.57006,
		0.566757,
		0.575431,
		0.578683,
		0.57858,
		0.58284,
		0.578724,
		0.582155,
		0.58572,
		0.581832,
		0.58098,
		0.579883,
		0.579096,
		0.58087};
	graph = new TGraph(25,Graph3_fx4,Graph3_fy4);
	graph->SetName("Graph3");
	graph->SetTitle("Graph");
	graph->SetFillColor(1);
	graph->SetLineWidth(3);
	graph->SetMarkerStyle(20);

	TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,280,6520);
	Graph_Graph4->SetMinimum(0.249142);
	Graph_Graph4->SetMaximum(0.616318);
	Graph_Graph4->SetDirectory(0);
	Graph_Graph4->SetStats(0);
	Graph_Graph4->SetLineWidth(3);
	Graph_Graph4->SetMarkerStyle(20);
	Graph_Graph4->GetXaxis()->SetLabelFont(42);
	Graph_Graph4->GetXaxis()->SetLabelOffset(0.007);
	Graph_Graph4->GetXaxis()->SetLabelSize(0.05);
	Graph_Graph4->GetXaxis()->SetTitleSize(0.05);
	Graph_Graph4->GetXaxis()->SetTitleFont(42);
	Graph_Graph4->GetYaxis()->SetLabelFont(42);
	Graph_Graph4->GetYaxis()->SetLabelOffset(0.007);
	Graph_Graph4->GetYaxis()->SetLabelSize(0.05);
	Graph_Graph4->GetYaxis()->SetTitleSize(0.05);
	Graph_Graph4->GetYaxis()->SetTitleOffset(1.1);
	Graph_Graph4->GetYaxis()->SetTitleFont(42);
	Graph_Graph4->GetZaxis()->SetLabelFont(42);
	Graph_Graph4->GetZaxis()->SetLabelOffset(0.007);
	Graph_Graph4->GetZaxis()->SetLabelSize(0.05);
	Graph_Graph4->GetZaxis()->SetTitleSize(0.05);
	Graph_Graph4->GetZaxis()->SetTitleFont(42);
	graph->SetHistogram(Graph_Graph4);

	graph->Draw("l");

	TLegend *leg = new TLegend(0.3,0.8,0.85,0.9,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.032);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	TLegendEntry *entry=leg->AddEntry("Graph0","Electron BB","L");

	ci = TColor::GetColor("#00ff00");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph1","Electron EB","L");

	ci = TColor::GetColor("#ff0000");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph2","Electron EE","L");

	ci = TColor::GetColor("#0000ff");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(3);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(21);
	entry->SetMarkerSize(1);
	entry->SetTextFont(42);
	entry=leg->AddEntry("Graph3","Electron Global","L");
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
