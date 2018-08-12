#!/usr/bin/Rscript

# Log: 11-08: Added Native Reference as parameter

library (bio3d)
library (parallel)
options (warn=0)

MC_CORES=1

# Creates a .pdf plot for a protein trajectory o pathway
# The input is either a compressed (.tgz) trajectory or
# a directory name with the PDBs files inside it.
USAGE="plot-pathway-tmscore.R <input pathway dir> [output dir]\n"
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function () {
	args = commandArgs (TRUE)
	if (length (args) < 0) {
		cat (USAGE)
		quit ()
	}
	pathname  = args [1]
	outputDir = getwd ()
	if (length (args) == 2)
		outputDir = args [2]

	filenames  = getInOutNames (pathname, outputDir)
	inputDir   = filenames$inputDir
	outputFile = filenames$outputFile

	# Extract or load filename to calculate RMSDs
	cat ("\nLoading PDBs from ", pathname, "...\n")
	files      = getPDBFiles (pathname)
	
	cat ("\nCalculating RMSDs...")	
	values     = parallelTmscorePathway (files$n, files$native, files$pdbs)

	cat ("\nWriting output file ", outputFile, "\n")
	plotPathway (values, outputFile)
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
parallelTmscorePathway <- function (n, native, pdbs) {
	values = mclapply (pdbs, calculateTmscore, native, 
						   mc.preschedule=T, mc.set.seed=T, mc.cores=30)
	return (values)
}
	
#----------------------------------------------------------
# Calculate the TM-scores using a external tool
#----------------------------------------------------------
calculateTmscore <- function (targetProtein, referenceProtein) {
	allArgs = c ("tmscore_zhang.py", referenceProtein, targetProtein)
	output  = system2 ("python", args=allArgs, stdout=T)
	return  (as.double (output))
}

#-------------------------------------------------------------
# Creata a XY plot from the RMSD values of each conformation
#-------------------------------------------------------------
plotPathway <- function (rmsdValues, outputFile) {
	pdf (outputFile, width=14)
		n = length(rmsdValues)
		rd = rmsdValues[1:(n-1)]
		time = 1:(n-1)
		plot(time, rd, typ = "l", ylab = "TM-score", xlab = "Frame No.")
		points (lowess(time,rd, f=2/10), typ="l", col="red", lty=2, lwd=2)
		#steps = n / 21
		#xPoints = seq (0,n, ceiling (steps))
		#axis (side=1, xPoints)
	dev.off ()
}

#--------------------------------------------------------------
# Get the PDB files from either a compressed file or a dir
#--------------------------------------------------------------
getPDBFiles <- function (pathname) {

	# Extracts to an inputDir 
	if (grepl ("gz", pathname)==T) {
		stemName = strsplit (pathname, split="[.]")[[1]][1]
		inputDir = stemName 
		untar (pathname, compressed=T, exdir=inputDir)
	} else 
		inputDir = pathname

	#inputFiles = list.files (inputDir)
	#inputFilesFull = sapply (inputFiles, function (x) paste (inputDir, x, sep="/"))
	inputFilesFull = list.files (inputDir, full.names=T)
	n = length (inputFilesFull)
	native = inputFilesFull [[n]]
	outfile = paste (inputDir, ".pdf", sep="")

	return (list (n=n, native=native, pdbs=inputFilesFull, outfile=outfile))
}

#--------------------------------------------------------------
#--------------------------------------------------------------
getInOutNames <- function (pathname, outputDir) {
	if (grepl ("gz", pathname)==T) 
		inputDir = strsplit (pathname, split="[.]")[[1]][1]
	else  
		inputDir = pathname

	outputFile = sprintf ("%s/%s.pdf", outputDir, basename (inputDir))
	return (list (inputDir=inputDir, outputFile=outputFile))
}
#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------
main ()
