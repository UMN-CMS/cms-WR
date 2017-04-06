import os
import re


print "#M_WR M_NU XS(pb) XS_ERR(pb)"
xspattern = re.compile("Before Filtrer: total cross section = (.*) \+- (.*) pb")
#epattern  = re.compile(".* (.*)   (.*)    0      -11  9900012")
#mupattern = re.compile(".* (.*)   (.*)    0      -13  9900014")
logdirs = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment/newestTwoDimLimitJobs/CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR/scripts/genAndFullOfflineAnalysisJobs/batchSubmFilesAndLogDirs/"

for ch, channeldir in  [("ee","eeChannel/"), ("mumu","mumuChannel/")]:
	for logdir in os.listdir(logdirs + channeldir):
		if "LSFJOB" in logdir:
			status = 0
			with open(logdirs + channeldir + logdir + "/STDOUT") as log:
				for line in log:
					split = line.split()
					if not split: continue
					if split[0] == "9900012":
						mnu = int(float(split[5]))
						status += 1
					elif split[0] == "9900024":
						mwr = int(float(split[6]))
						status += 2
					m = xspattern.match(line)
					if m:
						xs = m.group(1)
						xserr = m.group(2)
						status += 4

					if status == 1 + 2 + 4:
	
						print ch, mwr, mnu, xs, xserr
						break
				if status != 1 + 2 + 4:
					print logdirs + logdir






