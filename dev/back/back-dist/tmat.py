#!/usr/bin/python
import os, sys
import subprocess

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def printMatrix (matrix):
	n = len (matrix)
	for i in range (n):
		m = len (matrix[i])
		for j in range (m):
			print matrix [i][j],"\t",
		print ""

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def getTMScore (pdb1, pdb2):
	output = runProgram (["TMscore %s %s" % (pdb1, pdb2)])
	value  = preprocesOutputTM (output) 
	return  value
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def preprocesOutputTM (output):
	value = -1
	for line in output.split("\n"):
		if len (line.strip()) == 0:
			continue

		if line.split()[0] == "TM-score":
			value = float (line.split()[2])
			break
	return value

#--------------------------------------------------------------------
# call a external programm that returns a value running on "workingDir"
#--------------------------------------------------------------------
def runProgram (listOfParams):
	try:
		value = subprocess.Popen (listOfParams, stdout=subprocess.PIPE, shell=True).communicate()[0]
		return value
	except:
		print("runProgram Error:" + str(listOfParams))

