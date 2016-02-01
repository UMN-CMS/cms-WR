
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
import ROOT 

class limit1d:
	plus1sig = []
	min1sig = []
	plus2sig = []
	min2sig = []
	means = []
	medians = []
	errors = []
	masses = []
	def __init__(self,name, xs=None):
		self.name = name
		if xs is None:
			xs = 1
		self.xs = xs
	def add(self, mass, point):
		median, (mean, meanError), (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = point
		self.masses.append(mass)
		self.means.append(mean)
		self.medians.append(median)
		self.errors.append(meanError)
		self.plus1sig.append(onesig_plus)
		self.min1sig.append(onesig_minus)
		self.plus2sig.append(twosig_plus)
		self.min2sig.append(twosig_minus)
	def addTheory(self, theory):
		self.theory = theory

	def plot(self, filename, x_title="", y_title="", x_limits=(600,6000), y_limits=(1e-3, 1e-1)):
		#customROOTstyle()
		ROOT.gROOT.LoadMacro("scripts/tdrStyle.C")
		c1 = ROOT.TCanvas("c1","c1",800,800);
		c1.SetTopMargin(0.05);
		c1.SetLeftMargin(0.15);
		c1.SetRightMargin(0.05);
		c1.SetBottomMargin(0.13);
		c1.SetTicks(1,1);
		c1.SetLogy();

		mass_array = np.array(self.masses, dtype=float)
		expected_limit_array = np.array(self.medians, dtype=float)*self.xs
		expected_limit_error_array = np.array(self.errors, dtype=float)*self.xs
		onesig_array = np.array(self.plus1sig + self.min1sig[::-1] + [self.plus1sig[0]], dtype=float)*self.xs
		twosig_array = np.array(self.plus2sig + self.min2sig[::-1] + [self.plus2sig[0]], dtype=float)*self.xs
		mass_band_array = np.array(self.masses + self.masses[::-1] + [self.masses[0]], dtype=float)

		dummy = ROOT.TH1F("dummy","",30,x_limits[0], x_limits[1])
		dummy.SetMinimum(y_limits[0])
		dummy.SetMaximum(y_limits[1])
		dummy.SetStats(0)
		dummy.GetXaxis().SetNdivisions(507)
		dummy.GetXaxis().SetTitle(x_title)
		dummy.GetYaxis().SetTitle(y_title)

		dummy.GetYaxis().SetTitleOffset(1.1)
		dummy.Draw("HIST")
		leg = ROOT.TLegend(0.55,0.66,0.99,0.87)
		latex = ROOT.TLatex(0.57, 0.89, "M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
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
		g_twosig.SetFillColor(ROOT.kYellow)
		g_twosig.SetLineWidth(0)
		g_twosig.Draw("F SAME")

		g_onesig = ROOT.TGraph(n*2+1, mass_band_array, onesig_array)
		g_onesig.SetFillColor(ROOT.kGreen)
		g_onesig.SetLineWidth(0)
		g_onesig.Draw("F SAME")

		print self.min2sig, self.plus2sig
		print twosig_array

		if self.theory:
			ntheory = len(self.theory)
			theory_mass = np.array(sorted(self.theory.keys()), dtype=float)
			theory_limit = np.array( [ self.theory[mass] for mass in theory_mass], dtype=float)
			g_theory = ROOT.TGraph(ntheory, theory_mass, theory_limit)
			g_theory.SetLineWidth(3);
			g_theory.SetLineColor(ROOT.kRed+2);
			g_theory.SetLineStyle(0);
			g_theory.SetFillStyle(3002);
			g_theory.SetFillColor(ROOT.kRed);
			g_theory.Draw("L SAME")

		g_exp = ROOT.TGraphErrors(n,mass_array, expected_limit_array, expected_limit_error_array);
		g_exp.SetLineWidth(2);
		g_exp.SetLineColor(ROOT.kBlue); g_exp.SetLineStyle(2);
		g_exp.Draw("L SAME");
		
		leg.AddEntry(g_exp, "Expected limit","L")
		leg.AddEntry(g_onesig, "Expected #pm 1 #sigma", "F")
		leg.AddEntry(g_twosig, "Expected #pm 2 #sigma", "F")
		leg.AddEntry(g_theory,"Theory (g_{R}= g_{L})","L")
		leg.Draw("same")


		text = ROOT.TText(0.72,0.97,"CMS Preliminary")
		text.SetNDC();
		text.SetTextFont(42);
		text.SetTextSize(0.05);
		text.Draw()

		c1.RedrawAxis()
		c1.SaveAs(filename)

