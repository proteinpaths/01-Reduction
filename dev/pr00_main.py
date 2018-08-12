#!/usr/bin/python
import os, sys, os.path
"""
LOG: r1.2 (Aug3): Added K for global clustering
"""

"""
 Given a trayectory it uses a fast clustering algorith to
 reduce the trayectory to only the main representatives.
 INPUT:  
   <inputDir>	    An input directory with the trayectory files 
   <outputDir>	    An output directory with the results (representative files)
   <RMSD threshold> Threshold for local reduction comparisons that uses RMSD
   <Size Bin>	    Number of structures for each bin of the partiioned trajectory
   <N Cores>	    Number of cores to use for paralelizing the procedure
"""
USAGE  = "\nReduces a trayectory using a fast clustering"
USAGE += "\nUSAGE   : pr00_main.py <inputDir> <outputDir> <SizeBin> <TMSCORE> <K> <nCores>"
USAGE += "\nExample : pr00_main.py in out 40 0.5 5 4"

# Default values
threshold = 0.5   # TM-score threshold for comparisons between protein structures"
SIZEBIN   = 40	# Number of files for each bin
NCORES	  = 2	  # Number of cores for multiprocessing

def main (args):
	print (len (args))
	if len (args) < 6:
		print USAGE
		sys.exit (1)

	inputDir  = args [1]
	outputDir = args [2]
	SIZEBIN   = int (args [3])
	threshold = float (args [4])
	K         = int (args [5])
	NCORES	  = int (args [6])

	outputDirBins	= "%s/bins" % outputDir
	inputDirGlobal	= outputDir + "/binsLocal"
	outputDirGlobal = outputDir + "/pdbsGlobal"
	tmpDir          = outputDir 

	print "Parameters:"
	print "\t Input dir:          ", inputDir
	print "\t Output dir bins:    ", outputDirBins
	print "\t Size of bins: i     ", SIZEBIN
	print "\t TM-score threshold: ", threshold
	print "\t K:                  ", K
	print "\t NUM CORES :         ", NCORES
	print "\n"

	createDir (outputDir)
	writeLog (outputDir, inputDir, SIZEBIN, threshold, K)

	# Split full trajectory in bins (blocks of 1000 pdbs)
	cmm ="pr01_createBins.py %s %s %s" % (inputDir, outputDirBins, SIZEBIN)
	os.system (cmm) 

	# Get Representatives for each bin
	cmm = "pr02_localReduction.R %s %s %s %s" % (outputDirBins, outputDir, threshold, NCORES)
	os.system (cmm)

	# Uses K-Medoids to select K medoids for each local bin
	cmm = "pr03_globalReduction.R %s %s %s %s %s" % (inputDirGlobal, outputDirGlobal, tmpDir, K, NCORES)
	os.system (cmm)

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

#------------------------------------------------------------------
# Call main with input parameter
#------------------------------------------------------------------
#------------------------------------------------------------------
def writeLog (outputDir, proteinName, sizeOfBins, threshold, k):   
	inFile = open (outputDir+"/params.txt", "a")
	inFile.write ("Protein name: " + os.path.basename (proteinName) + "\n")
	inFile.write ("Size of bins : " + str (sizeOfBins) + "\n")
	inFile.write ("TM-score trhreshold : " + str (threshold) + "\n")
	inFile.write ("K representatives : " + str (k) + "\n")
	inFile.close ()
	
      
     
#------------------------------------------------------------------
if __name__ == "__main__":
	main (sys.argv)
