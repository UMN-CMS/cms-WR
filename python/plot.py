
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
import ROOT 

class limit1d:
	def __init__(self,channel,xs=None):
		self.channel = channel
		self.plus1sig = []
		self.min1sig = []
		self.plus2sig = []
		self.min2sig = []
		self.medians = []
		self.errors = []
		self.masses = []
		self.observed = []
		self.theory = None
		self.obslines = None
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
	def addObsLines(self, obs):
		self.obslines = obs

	def plot(self, filename, x_title="", y_title="",
			x_limits=(600,6000), y_limits=(1e-3, 1e-1), 
			leg_x = .55, leg_y=.66, leg_h=.32, SetLogx=False):

		customROOTstyle()
		#ROOT.gROOT.LoadMacro("scripts/tdrStyle.C")
		c1 = ROOT.TCanvas("c1","c1",800,800);
		c1.SetTopMargin(0.05);
		c1.SetLeftMargin(0.17);
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
			#uncomment the next three lines if the limits were calculated with unscale_by_xs in makeDatacards.py is False
			expected_limit_array = np.array([ l*self.theory[mass] for mass, l in zip(mass_array, expected_limit_array)], dtype=float)
			onesig_array         = np.array([ l*self.theory[mass] for mass, l in zip(mass_band_array, onesig_array)   ], dtype=float) 
			twosig_array         = np.array([ l*self.theory[mass] for mass, l in zip(mass_band_array, twosig_array)   ], dtype=float) 
			#expected_limit_array *= self.xs
			#onesig_array *= self.xs
			#twosig_array *= self.xs

		dummy = ROOT.TH1F("dummy","",30,x_limits[0], x_limits[1])
		dummy.SetMinimum(y_limits[0])
		dummy.SetMaximum(y_limits[1])
		dummy.SetStats(0)
		dummy.GetXaxis().SetNdivisions(507)
		dummy.GetXaxis().SetTitle(x_title)
		dummy.GetYaxis().SetTitle(y_title)

		dummy.GetYaxis().SetTitleOffset(1.6)
		dummy.Draw("HIST")
		leg_w = .44
		leg = ROOT.TLegend(leg_x,leg_y,leg_x + leg_w, leg_y + leg_h)
		latex = ROOT.TLatex(leg_x + 0.02, leg_y + leg_h + 0.02, "M_{#scale[1.25]{N_{#scale[1.5]{%s}}}}= M_{#scale[1.25]{W_{R}}}/2" % self.channel);
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

		g_exp = ROOT.TGraph(n,mass_array, expected_limit_array);
		g_exp.SetLineWidth(2);
		g_exp.SetLineColor(ROOT.kBlue); g_exp.SetLineStyle(2);
		g_exp.Draw("L SAME");

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

		if self.observed:
			nobs = len(self.observed)
			observed_limit_array = np.array( self.observed, dtype=float)
			if self.theory:
				#observed_limit_array = np.array([ l*self.theory[mass] for mass, l in zip(mass_array, observed_limit_array)], dtype=float)
				observed_limit_array *= self.xs
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

		if self.obslines:
			#this makes several lines on the 1D expected limit plot showing
			#the observed limit from 3.8 to 6.0 TeV for different numbers (d) of observed events
			lines = []
			for i in range(len(self.obslines)):
				obs = self.obslines[i]
				print i, obs
				line = ROOT.TLine(3800, obs*self.xs, 6000, obs*self.xs)
				line.SetLineWidth(2)
				line.SetLineStyle( [1, 9, 10][i % 3] )
				line.Draw()
				lines.append(line)
				leg.AddEntry(line, "Observed = %d" % i, "L")

		leg.Draw("same")

		latex = ROOT.TLatex()
		latex.SetNDC(True)
		latex.SetTextSize(0.03)
		latex.SetTextFont(42)
		latex.DrawLatex(0.72, 0.96, "2.64 fb^{-1} (13 TeV)")

		latex.DrawLatex(leg_x + .01, leg_y + leg_h + .01, "M_{W_{R}} = 2 M_{N_{l}}")

		text = ROOT.TText(0.17,0.96,"CMS Preliminary")
		text.SetNDC();
		text.SetTextFont(42);
		text.SetTextSize(0.03);
		text.Draw()

		c1.RedrawAxis()
		c1.SaveAs(filename + ".png")
		c1.SaveAs(filename + ".pdf")
		c1.SaveAs(filename + ".C")

from collections import defaultdict
class limit2d:
	def __init__(self, ratio_file, name,channelname,  mwr_range=(27 , 800, 6200) , mnu_range = (124, 0, 6200 ), xs=None):
		self.eff = defaultdict(dict)
		self.effratio = ROOT.TH2F(name + "_eff_ratio", name + " eff ratio" , *(mwr_range + mnu_range))
		self.limits = ROOT.TH2F(name + "_limits", name + " limits", *(mwr_range + mnu_range))
		self.r = ROOT.TH2F(name + "_r", name + " r", *(mwr_range + mnu_range))
		self.exclusion = ROOT.TH2F(name + "_exc", name + " exclusion", *(mwr_range + mnu_range))
		self.crosssection = ROOT.TH2F(name + "_xs", name + " xs", *(mwr_range + mnu_range))
		self.channelname = channelname
		self.exclusionTwo = ROOT.TH2F(name + "_excTwo", name + " exclusionTwo", *(mwr_range + mnu_range))
	
		with open(ratio_file,"r") as inf:
			for line in inf:
				line = line.strip()
				if "#" == line[0]: continue
				mwr, mnu, eff = map(float,line.split())
				self.eff[int(mwr)][int(mnu)] = eff
		if xs is None:
			xs = 1
		self.xs = xs

	def addTheory(self, theory):
		self.theory = theory
	
	def addObserved(self, mwr, point):
		mwr = int(mwr)
		median = float(point)
		#median *= self.xs
		for mnu in self.eff[mwr]:

			ratio = self.eff[mwr][mnu]/self.eff[mwr][mwr/2]
			if np.isnan(ratio): continue
			if ratio == 0:
				ratio = .00001

			limit = median * self.xs / ratio
			#self.effratio.Fill(mwr, mnu, ratio)
			#self.limits.Fill(mwr, mnu, limit)
			#self.r.Fill(mwr, mnu, median * self.xs)
			#self.exclusion.Fill(mwr, mnu, limit/self.theory[(mwr,mnu)])
			self.exclusionTwo.Fill(mwr, mnu, limit/self.theory[(mwr,mnu)])
			#self.crosssection.Fill(mwr, mnu, self.theory[(mwr,mnu)])

	def add(self, mwr, point):
		mwr = int(mwr)
		median, (onesig_minus,onesig_plus), (twosig_minus,twosig_plus) = point
		for mnu in self.eff[mwr]:

			ratio = self.eff[mwr][mnu]/self.eff[mwr][mwr/2]
			if np.isnan(ratio): continue
			if ratio == 0:
				ratio = .00001

			limit = median * self.xs / ratio
			self.effratio.Fill(mwr, mnu, ratio)
			self.limits.Fill(mwr, mnu, limit)
			#self.r.Fill(mwr, mnu, median * self.xs)
			self.exclusion.Fill(mwr, mnu, limit/self.theory[(mwr,mnu)])
			self.crosssection.Fill(mwr, mnu, self.theory[(mwr,mnu)])
	def draw(self, h, filename, zrange, ztitle, logz=False, cont=None, isObserved=False):
		customROOTstyle()
		ROOT.gStyle.SetOptTitle(0)
	
		c1 = ROOT.TCanvas("c1","c1",800,800);
		c1.SetTopMargin(0.05);
		c1.SetLeftMargin(0.2);
		c1.SetRightMargin(0.2);
		c1.SetBottomMargin(0.13);
		c1.SetTicks(1,1);
		#c1.SetLogz()

		c1.SetLogz(logz)
		h.Draw("colz")
		h.SetAxisRange(zrange[0], zrange[1],"Z")
		h.SetXTitle("W_{R} Mass [GeV]")
		h.SetYTitle("N_{l} Mass [GeV]")
		h.SetZTitle(ztitle)
		h.GetYaxis().SetTitleOffset(2)
		h.GetZaxis().SetTitleOffset(1.7)

		if cont:
			x = np.array([700,4050,700,700],dtype=float)
			y = np.array([700,4050,4050,700],dtype=float)
			#x = np.array([700,3700,700,700],dtype=float)
			#y = np.array([700,3700,3700,700],dtype=float)
			area = ROOT.TPolyLine(4,x,y)
			area.SetFillColor(ROOT.kYellow)
			area.SetLineWidth(0)
			area.Draw("F")
			latex2 = ROOT.TLatex()
			latex2.SetTextSize(0.05)
			#specify the lower left corner in x, y coordinates
			latex2.DrawLatex(1400,3200, "M_{N_{l}} > M_{W_{R}} ")
			cont.Draw("L")

			x1 = 900
			y1 = 3700
			xw = 1200
			yw = 200

			leg = ROOT.TLegend(x1, y1, x1 + xw, y1 + yw,"","")
			if not isObserved: leg.AddEntry(cont, "Expected","l")
			if isObserved: leg.AddEntry(cont, "Observed","l")
			leg.SetTextFont(42);
			leg.SetTextSize(0.032);
			leg.SetFillStyle(0);
			leg.SetBorderSize(0);
			
			leg.Draw()


		latex = ROOT.TLatex()
		latex.SetNDC(True)
		latex.SetTextSize(0.03)
		latex.SetTextFont(42)
		latex.DrawLatex(0.62, 0.96, "2.64 fb^{-1} (13 TeV)")

		text = ROOT.TText(0.2,0.96,"CMS Preliminary")
		text.SetNDC();
		text.SetTextFont(42);
		text.SetTextSize(0.03);
		text.Draw()
		c1.RedrawAxis()
		c1.SaveAs(filename + ".png")
		c1.SaveAs(filename + ".pdf")
		c1.SaveAs(filename + ".C")

	def plot(self,filename):

		#exclusion and exclusionTwo have the same binning structure
		nx = self.exclusion.GetNbinsX()
		ny = self.exclusion.GetNbinsY()
		for ix in range(1, nx+1):
			for iy in range(1,ny+1):
				if not self.exclusionTwo.GetBinContent(ix,iy):
					self.exclusionTwo.SetBinContent(ix, iy, 1000);
				if not self.exclusion.GetBinContent(ix,iy):
					self.exclusion.SetBinContent(ix, iy, 1000);
				mwr = self.exclusion.GetXaxis().GetBinCenter(ix)
				mnr = self.exclusion.GetYaxis().GetBinCenter(iy)
				mwrTwo = self.exclusionTwo.GetXaxis().GetBinCenter(ix)
				mnrTwo = self.exclusionTwo.GetYaxis().GetBinCenter(iy)
				if mnr > mwr:
					self.exclusion.SetBinContent(ix, iy, 2);
				elif mnr == mwr:
					self.exclusion.SetBinContent(ix, iy, 1.5);
				if mnrTwo > mwrTwo:
					self.exclusionTwo.SetBinContent(ix, iy, 2);
				elif mnrTwo == mwrTwo:
					self.exclusionTwo.SetBinContent(ix, iy, 1.5);

		
		graphs = contourFromTH2(self.exclusion, 1)
		graphsTwo = contourFromTH2(self.exclusionTwo, 1)
		c1 = ROOT.TCanvas("c1","c1",800,800);

		self.draw(self.exclusion,    filename + "_exclusion", (0   , 3), "Limit / #sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow %sqq)" % self.channelname, logz=False, cont = graphs[0], isObserved=False)
		self.draw(self.exclusionTwo,     filename + "_exclusionTwo", (0   , 3), "Limit / #sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow %sqq)" % self.channelname, logz=False, cont = graphsTwo[0], isObserved=True)
		#self.draw(self.limits,       filename + "_limit",     (1e-3, 1), "#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow %sjj) [pb]" % self.channelname,    logz=True,  cont = graphs[0])
		#self.draw(self.effratio,     filename + "_effratio",  (0   , 2), "efficiency #times acceptance (W_{R}, N_{l}) / (W_{R}, W_{R}/2) ",                         logz=False )
		#self.draw(self.crosssection, filename + "_xs",        (1e-3, 1), "#sigma(pp#rightarrow W_{R}) #times BR(%sjj) [pb]" % self.channelname,                     logz=True  )





def contourFromTH2(h2in, threshold, minPoints=20):
	print "Getting contour at threshold ", threshold, " from ", h2in.GetName()
	#http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
	contours = np.array([threshold], dtype=float)
	
	if (h2in.GetNbinsX() * h2in.GetNbinsY() > 10000): minPoints = 50;
	if (h2in.GetNbinsX() * h2in.GetNbinsY() <= 100): minPoints = 10;
	
	h2 = frameTH2D(h2in,threshold);
	
	h2.SetContour(1, contours);
	
	# Draw contours as filled regions, and Save points
	h2.Draw("CONT Z LIST");
	ROOT.gPad.Update(); # Needed to force the plotting and retrieve the contours in TGraphs
	
	
	# Get Contours
	conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
	
	if not conts: 
		 print "*** No Contours Were Extracted!"
		 return None
	else:
		 print "found contors"
	
	ret = []
	for i in range(conts.GetSize()):
		contLevel = conts.At(i)
		#printf("Contour %d has %d Graphs\n", i, contLevel.GetSize());
		for j in range(contLevel.GetSize()):
			gr1 = contLevel.At(j)
			#printf("\t Graph %d has %d points\n", j, gr1.GetN());
			if (gr1.GetN() > minPoints):
				ret.append(gr1.Clone())
			#break;
		return ret;

def frameTH2D(h, threshold):
	# NEW LOGIC:
	#   - pretend that the center of the last bin is on the border if the frame
	#   - add one tiny frame with huge values
	frameValue = 1000;
	#if (TString(h.GetName()).Contains("bayes")) frameValue = -1000;
	
	xw = h.GetXaxis().GetBinWidth(1);
	yw = h.GetYaxis().GetBinWidth(1);
	
	nx = h.GetNbinsX();
	ny = h.GetNbinsY();
	
	x0 = h.GetXaxis().GetXmin();
	x1 = h.GetXaxis().GetXmax();
	
	y0 = h.GetYaxis().GetXmin();
	y1 = h.GetYaxis().GetXmax();
	xbins = np.zeros(999,dtype=float)
	ybins = np.zeros(999,dtype=float) 
	eps = 0.1;
	
	xbins[0] = x0 - eps*xw - xw
	xbins[1] = x0 + eps*xw - xw
	#for (int ix = 2; ix <= nx; ++ix):
	for ix in range(2,nx+1):
		xbins[ix] = x0 + (ix-1)*xw

	xbins[nx+1] = x1 - eps*xw + 0.5*xw
	xbins[nx+2] = x1 + eps*xw + xw
	
	ybins[0] = y0 - eps*yw - yw;
	ybins[1] = y0 + eps*yw - yw;
	#for (int iy = 2; iy <= ny; ++iy)
	for iy in range(2,ny+1):
		ybins[iy] = y0 + (iy-1)*yw;

	ybins[ny+1] = y1 - eps*yw + yw;
	ybins[ny+2] = y1 + eps*yw + yw;
	
	framed = ROOT.TH2D( h.GetName() + "framed", h.GetTitle() + "framed", nx + 2, xbins, ny + 2, ybins )
	
	#Copy over the contents
	#for(int ix = 1; ix <= nx ; ix++){
	for ix in range(1, nx+1):
		#for(int iy = 1; iy <= ny ; iy++){
		for iy in range(1,ny+1):
			framed.SetBinContent(1+ix, 1+iy, h.GetBinContent(ix,iy));

	#Frame with huge values
	nx = framed.GetNbinsX();
	ny = framed.GetNbinsY();
	#for(int ix = 1; ix <= nx ; ix++){
	for ix in xrange(1, nx+1):
		framed.SetBinContent(ix,  1, frameValue);
		framed.SetBinContent(ix, ny, frameValue);

	#for(int iy = 2; iy <= ny-1 ; iy++){
	for iy in xrange(2, ny):
		framed.SetBinContent( 1, iy, frameValue);
		framed.SetBinContent(nx, iy, frameValue);
	
	return framed
