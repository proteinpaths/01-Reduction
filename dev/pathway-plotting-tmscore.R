#!/usr/bin/Rscript

# Log: 11-08:  
	# r1.4(Aug13) Fixed ranges when get values for plotting
	# r1.3(Aug13) Added a fixed range for y-axis
	# r1.2: Added a parameter for the native protein reference

#library (bio3d)
library (parallel)
options (warn=0)

nCores=1

# Creates a .pdf plot for a protein trajectory o pathway
# The input is either a compressed (.tgz) trajectory or
# a directory name with the PDBs files inside it.
USAGE="plot-pathway-tmscore.R <Reference Native> <input pathway dir> [nCores=1]\n"
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function () {
	args = commandArgs (TRUE)
	if (length (args) < 3) {
		cat (USAGE)
		quit ()
	}
	referenceProtein = args [1]
	pathname         = args [2]
	outputDir        = getwd ()
	if (length (args) == 3)
		nCores = args [3]

	filenames  = getInOutNames (pathname, outputDir)
	inputDir   = filenames$inputDir
	outputFile = filenames$outputFile

	# Extract or load filename to calculate RMSDs
	cat ("\nLoading PDBs from ", pathname, "...\n")
	files      = getPDBFiles (pathname, referenceProtein)
	
	cat ("\nCalculating TMscores")	
	values = parallelTmscorePathway (files$n, files$native, files$pdbs, nCores)

	cat ("\nWriting output file ", outputFile, "\n")
	plotPathway (values, outputFile)
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
parallelTmscorePathway <- function (n, native, pdbs, nCores) {
	values = mclapply (pdbs, calculateTmscore, native, 
					  mc.preschedule=T, mc.set.seed=T, mc.cores=nCores)
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
	#pdf (outputFile, width=20)
	pdf (outputFile, width=14)
		n = length(rmsdValues)
		rd = rmsdValues[1:n]
		time = 0:(n-1)
		plot(time, rd, typ = "l", ylab = "TM-score", xlab = "Frame No.", 
			 cex.axis=1.5,cex.lab=1.5,
		     mar=c(5,4,2,2)+0.4,
			 axes=TRUE, ylim=range(c(0,1)))

		#x = c(0.1, 0.2,0.3, 0.4, 0.6, 0.8, 1,2)
		#axis (side=2, at = x, labels=x)
		#points (lowess(time,rd, f=2/10), typ="l", col="red", lty=2, lwd=2)
		#steps = n / 21
		#xPoints = seq (0,n, ceiling (steps))
		#axis (side=1, xPoints)
	dev.off ()
}

#--------------------------------------------------------------
# Get the PDB files from either a compressed file or a dir
#--------------------------------------------------------------
getPDBFiles <- function (pathname, referenceProtein) {

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
	#native = inputFilesFull [[n]]
	outfile = paste (inputDir, ".pdf", sep="")

	return (list (n=n, native=referenceProtein, pdbs=inputFilesFull, outfile=outfile))
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
log <- function (msgs) {
	cat ("\nLOG: ")
	for (i in msgs){
		cat (i)
		cat ("  ")
	}
}

main ()
