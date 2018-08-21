#!/usr/bin/python

USAGE="\
Extracts from ini to end PDBs files from a directory \
USAGE: extract-pdbs.py <inputDir> <outputDir> <ini> <end>\n"

import os, sys, math
#SIZEBIN  = 100   # Number of files for each bin

#--------------------------------------------------------
# Main function to be called from command line
#--------------------------------------------------------
def main (args):
	if len (args) < 4:
		print USAGE
		sys.exit (1)

	inputDir  = "%s/%s" % (os.getcwd (), args [1])
	outputDir = "%s/%s" % (os.getcwd (), args [2])
	ini   = int (args [3]) 
	end   = int (args [4]) 

	createDir (outputDir)

	copyFiles (inputDir, outputDir, ini, end)

#--------------------------------------------------------
# Creates bins from files in an input dir 
# outputDir: destiny dir for bins
# binSize is the number of file by bin
# sizeFill is the prefix for each bin filename
#--------------------------------------------------------
def copyFiles (inputDir, outputDir, ini, end):
	inputFiles  = getSortedFilesDir (inputDir, ".pdb")
	n = len (inputFiles)
	subList = inputFiles [ini-1:end]

	for filename in subList:
		sourceFilename  = "%s/%s" % (inputDir, filename)
		destinyFilename = "%s/%s" % (outputDir, filename)
		#os.symlink (sourceFilename, destinyFilename)
		os.system ("cp %s %s" % (sourceFilename, destinyFilename))

#------------------------------------------------------------------
# Utility to create a directory safely.
# If it exists it is renamed as old-dir 
#------------------------------------------------------------------
def createDir (dir):
	def checkExistingDir (dir):
		if os.path.lexists (dir):
			headDir, tailDir = os.path.split (dir)
			oldDir = os.path.join (headDir, "old-" + tailDir)
			if os.path.lexists (oldDir):
					checkExistingDir (oldDir)

			os.rename (dir, oldDir)
	checkExistingDir (dir)
	os.system ("mkdir %s" % dir)
#--------------------------------------------------------------------
# Get the files containing the pattern from a inputDir 
#--------------------------------------------------------------------
def getSortedFilesDir (inputDir, pattern=""):
	files  = [x for x in os.listdir (inputDir) if pattern in x ]
	return sorted (files)
#--------------------------------------------------------------------
# Call main with input parameter
#--------------------------------------------------------------------
if __name__ == "__main__":
	main (sys.argv)
