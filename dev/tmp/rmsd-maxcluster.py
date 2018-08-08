#!/usr/bin/python

USAGE  = "Calculate the RMSD value using the maxcluster tool (maxcluster64bit)\n"
USAGE += "USAGE: rmsd-maxcluster.py <PDB reference>  <PDB structure> [outputFilename]\n"

import os, sys

###############################################################################
# call a external programm that returns a value running on "workingDir"
###############################################################################
def runProgram (listOfParams, workingDir):
	import subprocess
	value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
	return value

############## CHECK ARGUMENTS ######################
if len (sys.argv) < 3:
	print USAGE
	sys.exit (0)
elif len (sys.argv) == 4:
	sys.stdout = open (sys.argv[3], "w")

############## MAIN #################################
pdbReference = sys.argv [1]
pdbFilename = sys.argv [2]

currentDir = os.getcwd()

strvalue = runProgram (["maxcluster64bit", pdbReference, pdbFilename, "-rmsd"], currentDir)

value = strvalue.split ()[1]

print value

