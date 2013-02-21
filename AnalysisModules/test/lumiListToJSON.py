#!/usr/bin/env python 

import os
import sys
import commands
import fileinput
import string

from optparse import OptionParser                                                                                         
parser = OptionParser(usage="usage: %prog [options]")                                                                     
parser.add_option("--infile",  dest="infile",  help="Input lumiList file", type="string", default="lumilist.log")
(options, args) = parser.parse_args()

llfile = open(options.infile,'r')
goodLS = []

#--- Start of JSON ---#
sys.stdout.write("{")

lastRun = 0
for line in llfile:
    newRun  = line[line.find("'")+1:line.find(":")]
    firstLS = line[line.find(":")+1:line.find("-")]
    lastLS  = line[line.rfind(":")+1:line.rfind("'")]

    if newRun != lastRun:
        if lastRun != 0:
            print "], ",
        lastRun = newRun
        sys.stdout.write("\"" + newRun + "\": [")
    else:
        sys.stdout.write(", ")

    sys.stdout.write("[" + firstLS + ", " + lastLS + "]")


#--- End of JSON ---#
print "]}"
