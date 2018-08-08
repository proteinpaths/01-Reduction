#----------------------------------------------------------
#----------------------------------------------------------
testTmscore <- function (inputProteinsLst) {
	inputProteinsLst = list.files ("io/in", pattern="pdb", full.names=T)
	n = length (inputProteinsLst)
	p0 = inputProteinsLst [[136]]
	for (k in 2:(n)) {
		p1 = inputProteinsLst [[k]]
		tmscoreValue = runTmscore (p0,p1)
		cat ("\n", p0, p1, tmscoreValue)
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


