#!/usr/bin/env python

import os,sys
sys.argv.append('-b')
import ROOT

def print_histos_from_tree(root_file,dir,tree,formats):
    """Print all the histograms from a ROOT file
    in each of the formats."""
    
    f = ROOT.TFile(root_file)
    for x in f.GetListOfKeys():
        os.system('mkdir -p plots/print_histos/'+x.GetName()+'/'+root_file.split('/')[-1][:-5])
        d = f.GetDirectory(x.GetName())
        t = d.Get(tree)
        print t
        for y in t.GetListOfBranches():
            c1 = ROOT.TCanvas("c1","Canvas",600,600)        
            a = ''+y.GetName()
            t.Draw(a)
            for i in formats:
                c1.Print('plots/print_histos/'+x.GetName()+'/'+root_file.split('/')[-1][:-5]+'/'+ a +'.'+i)
            c1 = 0



print_histos_from_tree('flat_ttree.root','ft_nocuts','t',['pdf','png'])    
