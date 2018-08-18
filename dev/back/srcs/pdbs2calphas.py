#!/usr/bin/python

import os, sys 

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


#---- Main ----
inputDir  = sys.argv[1]
outputDir = "shm/alphas"
createDir (outputDir)

pdbFilenames = os.listdir (inputDir)

for pdb in pdbFilenames:
    nombre = pdb.split (".")[0]
    cmm = 'grep "  CA  " %s/%s > %s/%s-A.pdb &' % (inputDir, pdb, outputDir, nombre)
    print (cmm)
    os.system (cmm)


