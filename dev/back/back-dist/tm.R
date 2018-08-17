#!/usr/bin/Rscript

listOfPDBPaths = list()

getDistMat <- function (inputBinDir, outputBinDir, listOfPDBPaths) {
	n = length (listOfPDBPaths)
	mat = matrix (seq (1,n))

	distMat = proxy::dist (mat, method=calculateTmscore)
	return (distMat)
}

#--------------------------------------------------------------
# Calculate pairwise using TM-score distance
#--------------------------------------------------------------
getTMDistanceMatrix <- function (inputBinDir, outputBinDir) {
	system (paste ("TMscore-distance-matrix.py", inputBinDir, outputBinDir))
	distanceFilename = paste (outputBinDir,"/", basename (inputBinDir), ".dist", sep="")
	matrixTM = read.table (distanceFilename, header=T)
	distanceMatrixTM =  as.dist (matrixTM)
	return (distanceMatrixTM)
}


	
calculateTmscore <- function (targetProtein, referenceProtein) {
	targetProtein    = listOfPDBPaths [[targetProtein]]
	referenceProtein = listOfPDBPaths [[referenceProtein]]

	allArgs = c (referenceProtein, targetProtein)
	output  = system2 ("TMscore", args=allArgs, stdout=T)

	lineTMscore = strsplit (output, "\n")[[17]]
	tmscore = strsplit (lineTMscore, "[ ]{1,}")
	return  (1-as.double (tmscore[[1]][[3]]))
}


args <- commandArgs (TRUE)
INPUTDIR          = args [1] 
inputBinDir = INPUTDIR
listOfPDBPaths    = list.files (INPUTDIR, pattern="p", full.names=T)
print (listOfPDBPaths)

#cat ("Calculating distance matrix using library proxy...\n")
#m = getDistMat (inputBinDir, ".", listOfPDBPaths)
#print (m)

cat ("\nCalculating distance matrix using python script\n")
p = getTMDistanceMatrix (inputBinDir, "bintmp")
print (p)
