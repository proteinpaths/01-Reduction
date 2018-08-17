#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

# LOG: 
#	r1.4 (Aug16): Modified clustering with initial medoids including the last pdb
#	r1.3 (Aug13): Fixed error when it gets number of pdbs < K
#	r1.2 (Aug3): Extracts K medoids and users TM-score instead RMSD

#----------------------------------------------------------
# Makes a detailed global clustering of protein conformations 
# from the representatives resulting from the local clustering 
# INPUT:  inputDir filename with the protein conformations
#         outputDir filename to write the results
#
# OUTPUT: Medoids for each bin cluster and their distance matrix
#----------------------------------------------------------
USAGE="USAGE: global.R <inputDir> <outputDir> <tmpDir> <K> <num cores>\n" 

#library (bio3d)
library (parallel)
library (cluster)
options (width=300)

#THRESHOLD = 1.3
#NCORES= 1

#----------------------------------------------------------
# Main function
#----------------------------------------------------------
main <- function () {
	args <- commandArgs (TRUE)
	#args = c("io/out1000/outbins", "io/out1000/outrepr", "1.5", "1")
	print (args)
	if (length (args) < 5){
		cat (USAGE)
		cat (length(args))
		quit (status=1)
	}

	INPUTDIR  = args [1] 
	OUTPUTDIR = args [2]
	tmpDir    = paste (args [3], "/tmp", sep="")
	K         = as.numeric (args [4])
	NCORES    = as.numeric (args [5])

	createDir (OUTPUTDIR)
	createDir (tmpDir)

	print ("\n FULL")
	print (OUTPUTDIR)
	listOfBinPaths    = list.files (INPUTDIR, pattern="bin", full.names=T)
	clusteringResults = mclapply (listOfBinPaths, reduceGlobal, 
							outputDir=OUTPUTDIR, tmpDir=tmpDir, K, mc.cores=NCORES)

	cat ("\nRESULTS MCLAPPLY\n")
	writeClusteringResults (clusteringResults, OUTPUTDIR)
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
# Clustering around medoids. Return k medoid for the bin
#----------------------------------------------------------
reduceGlobal <- function (inputBinPath, outputDir, tmpDir, K) {
	cat ("\n>>> Global Reducing ", inputBinPath )

	# Fast clustering for bin, writes representatives to clusDir
	listOfPDBPaths <- list.files (inputBinPath, pattern=".pdb", full.names=T)

  # Clustering around medoids. Return one medoid for all inputs
	nPdbs = length (listOfPDBPaths)
	if (nPdbs < 2)
		medoids = 1
	else if (nPdbs <= K)
		medoids = seq (nPdbs)
	else {
		binDir = inputBinPath
		print (binDir)
		distanceMatrix <- getTMDistanceMatrix (binDir, tmpDir)
		split          <- -1 * nPdbs / K
		initialMedoids <- round (seq (nPdbs, 1, split))
		pamPDBs        <- pam (distanceMatrix, k=K, diss=F, medoids=initialMedoids)
		medoids        <- pamPDBs$id.med
	}

	medoidName <- listOfPDBPaths [medoids]
	return (medoidName)
}

#----------------------------------------------------------
# Make links of the selected PDBs into the output dir
#----------------------------------------------------------
writeClusteringResults <- function (clusteringResults, outputDir) {
	for (binResults in clusteringResults) 
		for (pdbPath in binResults) {
			cmm <- sprintf ("ln -s %s/%s %s/%s", getwd(),
									 pdbPath, outputDir, basename (pdbPath))
			#cat (paste (">>> ", cmm, "\n"))
			system (cmm)
		}
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


#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}

			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}


getDistMat <- function (inputBinDir, outputBinDir) {
	n = length (listOfPDBPaths)
	mat = matrix (seq (1,n))

	distMat = proxy::dist (mat, method=calculateTmscore)
	return (distMat)
}
	


calculateTmscore <- function (targetProtein, referenceProtein) {
	cat (">>>", targetProtein, "\n")
	cat (">>>", referenceProtein,"\n")
	targetProtein    = listOfPDBPaths [[targetProtein]]
	referenceProtein = listOfPDBPaths [[referenceProtein]]
	cat (">>>", targetProtein, "\n")
	cat (">>>", referenceProtein,"\n")

	allArgs = c (referenceProtein, targetProtein)
	output  = system2 ("TMscore", args=allArgs, stdout=T)

	lineTMscore = strsplit (output, "\n")[[17]]
	tmscore = strsplit (lineTMscore, "[ ]{1,}")
	return  (as.double (tmscore[[1]][[3]]))
}



#--------------------------------------------------------------
#--------------------------------------------------------------
main () 
	
