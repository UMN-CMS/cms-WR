
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
import ROOT 

class limit1d:
	def __init__(self,xs=None):
		self.plus1sig = []
		self.min1sig = []
		self.plus2sig = []
		self.min2sig = []
		self.medians = []
		self.errors = []
		self.masses = []
		self.observed = []
		self.theory = None
		if xs is None:
			xs = 1
		self.xs = xs
	def add(self, mass, point):
		median, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = point
		self.masses.append(float(mass))
		self.medians.append(float(median))
		self.plus1sig.append(float(onesig_plus))
		self.min1sig.append(float(onesig_minus))
		self.plus2sig.append(float(twosig_plus))
		self.min2sig.append(float(twosig_minus))
	def addObserved(self, mass, point):
		self.observed.append(float(point))
	def addTheory(self, theory):
		self.theory = theory

	def plot(self, filename, x_title="", y_title="",
			x_limits=(600,6000), y_limits=(1e-3, 1e-1), 
			leg_x = .55, leg_y=.66, SetLogx=False, scale=1):
		customROOTstyle()
		#ROOT.gROOT.LoadMacro("scripts/tdrStyle.C")
		c1 = ROOT.TCanvas("c1","c1",800,800);
		c1.SetTopMargin(0.05);
		c1.SetLeftMargin(0.15);
		c1.SetRightMargin(0.05);
		c1.SetBottomMargin(0.13);
		c1.SetTicks(1,1);
		c1.SetLogy();
		if SetLogx:
			c1.SetLogx()

		mass_array = np.array(self.masses, dtype=float)
		expected_limit_array = np.array(self.medians, dtype=float)
		expected_limit_error_array = np.array(self.errors, dtype=float)
		onesig_array = np.array(self.plus1sig + self.min1sig[::-1] + [self.plus1sig[0]], dtype=float)
		twosig_array = np.array(self.plus2sig + self.min2sig[::-1] + [self.plus2sig[0]], dtype=float)
		mass_band_array = np.array(self.masses + self.masses[::-1] + [self.masses[0]], dtype=float)

		if self.theory:
			#uncomment the following three lines if unscale_by_xs in scripts/makeDatacards.py is False
			#expected_limit_array = np.array([ l*self.theory[mass] for mass, l in zip(mass_array, expected_limit_array)], dtype=float)
			#onesig_array         = np.array([ l*self.theory[mass] for mass, l in zip(mass_band_array, onesig_array)   ], dtype=float) 
			#twosig_array         = np.array([ l*self.theory[mass] for mass, l in zip(mass_band_array, twosig_array)   ], dtype=float) 
			#uncomment the following three lines if unscale_by_xs in scripts/makeDatacards.py is True
			expected_limit_array *= self.xs
			onesig_array *= self.xs
			twosig_array *= self.xs

		#rescale expected limits and expected limit sigma bands
		expected_limit_array *= scale
		onesig_array *= scale
		twosig_array *= scale
		dummy = ROOT.TH1F("dummy","",30,x_limits[0], x_limits[1])
		dummy.SetMinimum(y_limits[0])
		dummy.SetMaximum(y_limits[1])
		dummy.SetStats(0)
		dummy.GetXaxis().SetNdivisions(507)
		dummy.GetXaxis().SetTitle(x_title)
		dummy.GetYaxis().SetTitle(y_title)

		dummy.GetYaxis().SetTitleOffset(1.1)
		dummy.Draw("HIST")
		leg_w = .44
		leg_h = .21
		leg = ROOT.TLegend(leg_x,leg_y,leg_x + leg_w, leg_y + leg_h)
		#for muon channel 1D limit plot
		#latex = ROOT.TLatex(leg_x + 0.02, leg_y + 0.23, "M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
		#for ele channel 1D limit plot
		latex = ROOT.TLatex(leg_x + 0.02, leg_y + 0.23, "M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
		latex.SetNDC();
		latex.SetTextSize(0.032);
		latex.SetTextFont(42);
		latex.Draw();

		leg.SetTextFont(42);
		leg.SetTextSize(0.032);
		leg.SetFillStyle(0);
		leg.SetBorderSize(0);

		n = len(self.masses)

		g_twosig = ROOT.TGraph(n*2+1, mass_band_array, twosig_array)
		#g_twosig.SetLineWidth(2)
		#g_twosig.SetLineColor(ROOT.kBlack)
		#g_twosig.Draw("L SAME")
		g_twosig.SetFillColor(ROOT.kYellow)
		g_twosig.Draw("F SAME")

		#comment this to remove 1sigma band
		g_onesig = ROOT.TGraph(n*2+1, mass_band_array, onesig_array)
		g_onesig.SetFillColor(ROOT.kGreen)
		g_onesig.SetLineWidth(0)
		g_onesig.Draw("F SAME")

		g_exp = ROOT.TGraph(n,mass_array, expected_limit_array);
		g_exp.SetLineWidth(2);
		g_exp.SetLineColor(ROOT.kBlue); g_exp.SetLineStyle(2);
		g_exp.Draw("L SAME");

		if self.theory:
			ntheory = len(self.theory)
			theory_mass = np.array(sorted(self.theory.keys()), dtype=float)
			theory_limit = np.array( [ self.theory[mass] for mass in theory_mass], dtype=float)
			#rescale theory limit
			theory_limit *= scale
			#theory_mass_unc = onesig_array*0
			#theory_limit_unc = np.array( [ self.theory[mass]*1.010 for mass in theory_mass ], dtype=float)
			#g_theory = ROOT.TGraphErrors(ntheory, theory_mass, theory_limit, theory_mass_unc, theory_limit_unc)  #TGraphErrors assumes errors are symmetric in X and Y at each point
			g_theory = ROOT.TGraph(ntheory, theory_mass, theory_limit)
			g_theory.SetLineWidth(3);
			g_theory.SetLineColor(ROOT.kRed+2);
			g_theory.SetLineStyle(0);
			g_theory.SetFillStyle(3002);
			g_theory.SetFillColor(ROOT.kRed);
			#g_theory.Draw("F SAME")
			g_theory.Draw("L SAME")

		if self.observed:
			nobs = len(self.observed)
			observed_limit_array = np.array( self.observed, dtype=float)
			if self.theory:
				#observed_limit_array = np.array([ l*self.theory[mass] for mass, l in zip(mass_array, observed_limit_array)], dtype=float)
				observed_limit_array *= self.xs
			#rescale observed limit
			observed_limit_array *= scale
			g_obs = ROOT.TGraph(nobs, mass_array, observed_limit_array)
			g_obs.SetLineWidth(3);
			g_obs.SetLineColor(ROOT.kBlue+2);
			g_obs.SetLineStyle(0);
			g_obs.SetFillStyle(3002);
			g_obs.SetFillColor(ROOT.kBlue+2);
			g_obs.Draw("L SAME")
		
		if self.observed:
			leg.AddEntry(g_obs,"Observed limit","L")
		leg.AddEntry(g_exp, "Expected limit","L")
		leg.AddEntry(g_onesig, "Expected #pm 1 #sigma", "F")
		leg.AddEntry(g_twosig, "Expected #pm 2 #sigma", "F")
		if self.theory:
			leg.AddEntry(g_theory,"Theory (g_{R}= g_{L})","L")
		leg.Draw("same")


		#text = ROOT.TText(0.32,0.97,"CMS Preliminary #surds = 13 TeV #int lumi = 2.6 fb^{-1}")
		#text.SetNDC();
		text = ROOT.TLatex()
		text.SetTextFont(42)
		text.SetTextSize(0.043)
		text.DrawLatexNDC(0.22,0.96,"CMS Preliminary        2.6 fb^{-1} (13 TeV)")
		#text.Draw()

		c1.RedrawAxis()
		c1.SaveAs(filename + ".png")
		c1.SaveAs(filename + ".pdf")
		c1.SaveAs(filename + ".C")
