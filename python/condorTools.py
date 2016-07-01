#!/usr/bin/python

#  Copyright (C) 2013  Alexander Gude - gude@physics.umn.edu
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#  The most recent version of this program is avaible at:
#  https://github.com/UMN-CMS/CondorUtilities


from os import environ, makedirs
from os.path import isfile, isdir, basename
from sys import exit
from math import ceil
from subprocess import call


class CondorFile:
    """ Generates a file to submit to condor """
    def __init__(self, condorFile, executable, logDir, niceUser):
        self.condorFile = condorFile
        self.executable = executable
        self.logDir = logDir
        self.niceUser = niceUser
        self.__setHeader()
        self.cont = self.header

    def __setHeader(self):
        """ Just a way to serperate the hardcoded strings from everything else.
        """
        # List of machines not to run on
        self.banned_machines = [
            "zebra01.spa.umn.edu",
            "zebra02.spa.umn.edu",
            "zebra03.spa.umn.edu",
            "zebra04.spa.umn.edu",
            "caffeine.spa.umn.edu",
            "gc1-ce.spa.umn.edu",
            "gc1-hn.spa.umn.edu",
            "gc1-se.spa.umn.edu",
            "red.spa.umn.edu",
            "hadoop-test.spa.umn.edu"
            ]

        # Set up our string
        self.header = ""  # Start blank so that we can always += append
        if self.niceUser:
            self.header += "nice_user = True\n"
        self.header += "Executable = %(executable)s\n" % {"executable": self.executable}
        self.header += "Universe = vanilla\n"
        self.header += "Output = %(logDir)s/output\n" % {"logDir": self.logDir}
        self.header += "Error = %(logDir)s/error\n" % {"logDir": self.logDir}
        self.header += "request_memory = 400\n"
        self.header += 'Requirements = (Arch=="X86_64")'
        # Insert banned machines
        for machine in self.banned_machines:
            self.header += ' && (Machine != "%s")' % machine
        self.header += '\n\n+CondorGroup="cmsfarm"\n\n'

    def addJob(self, localRT, jobDir, jobName, arguments):
        """ Add an 'Arguments' and a 'Queue' command to the condorfile. """
        self.cont += "Arguments = %(localRT)s %(jobDir)s %(jobName)s %(arguments)s\n" % {
            "localRT": localRT,
            "jobDir": jobDir,
            "jobName": jobName,
            "arguments": arguments,
            }
        self.cont += "Queue\n\n"

    def write(self):
        """ Save the condor file """
        # Check save location
        saveLocation = self.condorFile
        # Save file
        f = open(saveLocation, 'w')
        f.write(self.cont)
        f.flush()
        f.close()


class Job:
	def __init__(self, executable, jobName, prodSpace="/local/cms/user/" + environ["USER"], niceUser = True):
		# Check for critical environment variables, exit with an error if we don't
		# find them
		try:
			self.localRT = environ["LOCALRT"]
		except KeyError:
			exit("$LOCALRT not set. Remember to run 'cmsenv' in the right release area.")
		try:
			scramArch = environ["SCRAM_ARCH"]
		except KeyError:
			exit("$SCRAM_ARCH not set. Remember to run 'cmsenv' in the right release area.")

		self.executable = executable
		self.jobName = jobName
		self.prodSpace = prodSpace
		self.commands = []
		# Set up directories
		jobDir = prodSpace + '/' + jobName + '/'
		self.jobDir = jobDir
		logDir = jobDir + "log/"  # For logs from each job
		condorDir = jobDir + "condor/"
		condorFile = condorDir + jobName + ".condor"
		# Make directories if they do not already exist
		for d in [jobDir, logDir, condorDir]:
			if not isdir(d):
				makedirs(d)

		# Open files
		self.cf = CondorFile(condorFile, executable, logDir, niceUser)
		self.written = False
	

	def addJob(self, arguments, id_str):
		# Add job to condor file
		self.cf.addJob(self.localRT, self.jobDir, self.jobName + id_str , arguments)

		self.commands += ["%(exec)s %(localRT)s %(jobDir)s %(jobName)s %(arguments)s" % {
			"exec": self.executable,
			"localRT": self.localRT,
			"jobDir": self.jobDir,
			"jobName": self.jobName + id_str,
			"arguments": arguments,
			}]

	# Write Condor file
	def write(self):
		self.cf.write()
		self.written = True

	def submit(self, mode = "condor"):
		if mode == "interactive":
			for command in self.commands:
				call( command.split())
		elif mode == "lsf":
			bsub_prefix = "bsub -q cmscaf1nd "
			for command in self.commands:
				call( (bsub_prefix + command).split())
		elif mode == "condor":
			if not self.written: self.cf.write()
			retcode = call(["condor_submit", self.cf.condorFile])
			if retcode != 0:
				print "Error returned from condor_submit!"
				exit(retcode)
		else: return
