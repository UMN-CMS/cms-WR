

tex= "\\begin{{table}}[h]\n\\caption{{{caption}}}\n\\label{{tab:{label}}}\n\\centering\n\\begin{{tabular}}{{{col}}}\\hline\n{table}\\end{{tabular}}\n\\end{{table}}\n"
header = " & \\multicolumn{2}{c|}{Electrons}  & \\multicolumn{2}{c|}{Muons}  \\\\  \\hline\nMass & Low(GeV) & High(GeV) & Low(GeV) & High(GeV) \\\\  \\hline\n"
table = {}
masses = set()
with open("configs/mass_cuts.txt") as f:
	for line in f:
		if line[0] == "#": continue
		ch, mass, low, hi = line.split()
		table[(mass , ch)] = "& {low:<8} & {hi:<9}".format(low=low, hi=hi)
		masses.add(mass)

tex_table = header
for mass in sorted(masses):
	if int(mass) == 0: continue
	tex_table  += "{mass:4d} {ee} {mumu} \\\\ \\hline\n".format(mass=int(mass), ee=table[(mass, "EE")], mumu = table[(mass, "MuMu")])



print tex.format(caption="\\mW cuts that contain 90\% of signal events.", label="masscuts", table=tex_table,   col="|c|c|c|c|c|c|c|")
