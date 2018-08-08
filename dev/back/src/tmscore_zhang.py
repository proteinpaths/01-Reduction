#!/usr/bin/python

USAGE  = "Calculate the TM-Score value using the Zhang&Skolnit tool\n"
USAGE += "USAGE: tmscore-zhang.py <PDB reference>  <PDB structure> [outputFilename]\n"

import os, sys

############## MAIN #################################
def main ():
    args = sys.argv
    ## CHECK ARGUMENTS 
    if len (args) < 3:
            print USAGE
            sys.exit (0)
    elif len (args) == 4:
            sys.stdout = open (sys.argv[3], "w")


    pdbReference = args [1]
    pdbTarget  = args [2]
    currentDir = os.getcwd()

    value = tmscore (pdbReference, pdbTarget, currentDir)
    print >> sys.stderr, os.path.basename (pdbReference), os.path.basename(pdbTarget), value

    print value

#-------------------------------------------------------------
# Return the tmscore value by calling a external program
#-------------------------------------------------------------
def tmscore (pdbReference, pdbTarget, currentDir=None):
    if currentDir == None:
        currentDir = os.getcwd ()

    value = runProgram (["TMscore", pdbReference, pdbTarget, "-rmsd"], currentDir)
    return value

#-------------------------------------------------------------
# Return the tmscore value by calling a external program
# It runs the command and processes the output.
# It get the value using the command line arguments and 
# running on the working dir
#-------------------------------------------------------------
def runProgram (listOfParams, workingDir):
	import subprocess
	#txtLines = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
        txtLines = subprocess.check_output(listOfParams, universal_newlines=True)
        value = getValueFromTxtLine (txtLines)

	return value

############## Get Value ARGUMENTS ######################
def getValueFromTxtLine (txtLines):
    for line in txtLines.split("\n"):
        fields = line.split ()
        if "TM-score" in line and fields[1]=="=":
            value = fields [2]
            return float (value)
    return -9999.999



#------------------------------------------------------------------
# Call main with input parameter
#------------------------------------------------------------------
if __name__ == "__main__":
	main ()
