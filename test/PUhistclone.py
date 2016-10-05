import ROOT

f = ROOT.TFile('puAna.root')
g = ROOT.TFile('myMCPileup.root','RECREATE')

#t = f.Get('miniTree_noSelections/t')
#h = ROOT.TH1F('pileup','pileup',50,0,50)
#for e in t:
#    h.Fill(e.nPU)

h = f.Get('ana/pileup').Clone()

h.Draw()
g.Write()
