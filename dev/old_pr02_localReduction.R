#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

# LOG: r1.2 (Aug3): Using TM-score instead RMSD

#----------------------------------------------------------
# Make a fast clustering of protein conformations from a 
# trayectory by doing first a fast local clustering 
# INPUT:  
#   <inputDir>       An input directory with the trayectory files 
#   <outputDir>      An output directory with the results (representative files)
#   <RMSD threshold> Threshold for local reduction comparisons that uses RMSD
#   <Size Bin>       Number of structures for each bin of the partiioned trajectory
#   <N Cores>        Number of cores to use for paralelizing the procedure
#
# OUTPUT:            An output dir with the clustering results for each bin
#----------------------------------------------------------
USAGE="USAGE: reduction.R <inputDir> <outputDir> <TM-score Threshold> <num cores>\n" 

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
	if (length (args) < 4){
		cat (USAGE)
		quit (status=1)
	}

	INPUTDIR  = args [1] 
	OUTPUTDIR = args [2]
	THRESHOLD = as.numeric (args [3])
	NCORES    = as.numeric (args [4])
	dirBins = paste (OUTPUTDIR,"/binsLocal",sep="")
	dirPdbs = paste (OUTPUTDIR,"/pdbsLocal",sep="")

	createDir (dirBins)
	createDir (dirPdbs)

	binPathLst         = list.files (INPUTDIR, pattern="bin", full.names=T)
	clusteringResults  = mclapply (binPathLst, reduceLocal, 
								   outputDir=dirBins, 
								   threshold=THRESHOLD, mc.cores=NCORES)
	writeClusteringResults (clusteringResults, dirPdbs)
}


#----------------------------------------------------------
# Reduction function to reduce a single bin
#----------------------------------------------------------
reduceLocal <- function (inputBinPath, outputDir, threshold) {
	cat ("\n>>> Local Reducing ", inputBinPath )
	# Create the output dir for representatives
	outputBinPath = (paste (getwd(), outputDir, basename (inputBinPath), sep="/"))
	createDir (outputBinPath)		

	inputProteinsLst <- list.files (inputBinPath, full.names=T)

	# Fast clustering for bin, writes representatives to outputBinPath
	listOfSelectedPdbs = fastClustering (inputBinPath, outputBinPath, threshold, inputProteinsLst)
	return (listOfSelectedPdbs)
}

#----------------------------------------------------------
# Fast clustering following hobbohm algorith
# The first protein in the bin is selected as representative, then
# if the others are differente from it, they are preserved.
# Write the links to the representatives in the output dir
# Return a list with the representative as the first pdb in the group
#----------------------------------------------------------
fastClustering <- function (inputBinPath, outputBinPath, threshold, inputProteinsLst) {
	n = length (inputProteinsLst)
	targetProteinPath  = inputProteinsLst [[n]]
	headProteinPath    = targetProteinPath   # Default head for the first group
	tmscoreValue       = -1  # To create the link for the first group
	listOfSelectedPdbs = list()
	for (k in n:1) {
		if (tmscoreValue < threshold) {
			listOfSelectedPdbs = append (listOfSelectedPdbs, targetProteinPath)
			cmm = sprintf ("ln -s %s/%s %s/%s", getwd(), targetProteinPath, 
							outputBinPath, basename (targetProteinPath))
			#cat ("\n>>>",cmm); 
			system (cmm)
			headProteinPath = targetProteinPath
		}

		targetProteinPath = inputProteinsLst[[k]] 
		#cat ("\n\n>>> Protein: ", basename (targetProteinPath), "<<<")
		tmscoreValue = runTmscore (targetProteinPath, headProteinPath)
	}
	return (listOfSelectedPdbs)
}	

#----------------------------------------------------------
# Make links of the selected PDBs into the output dir
#----------------------------------------------------------
writeClusteringResults <- function (clusteringResults, outputDir) {
	cat ("\nWriting results local reducing...\n")
	for (binResults in clusteringResults) 
		for (pdbPath in binResults) {
			cmm <- sprintf ("ln -s %s/%s %s/%s", getwd(),
									 pdbPath, outputDir, basename (pdbPath))
			#cat (paste (">>> ", cmm, "\n"))
			system (cmm)
		}
}

#----------------------------------------------------------
# Calculate the TM-scores using a external tool
#----------------------------------------------------------
runTmscore <- function (referenceProtein, targetProtein) {
		allArgs = c ("tmscore_zhang.py", referenceProtein, targetProtein)
		output  = system2 ("python", args=allArgs, stdout=T)
		return  (as.double (output))
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

#--------------------------------------------------------------
#--------------------------------------------------------------
main () 
	
