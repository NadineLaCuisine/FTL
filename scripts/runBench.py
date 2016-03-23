#!/usr/bin/python

import os, glob, time, sys
import subprocess

#Function that launch FTL and parse the result (it use a time out)
def runBench(readFile, referenceFile, anchorSize, maxMissmatch):
	output = subprocess.Popen(["./ftl", "-u", readFile, "-x", referenceFile, "-k", str(anchorSize), "-m", str(maxMissmatch)], stdout=subprocess.PIPE)

	#Change here the timeout ! (it's in seconds)
	sleepingSecond = 60
	while output.poll() == None:
		time.sleep(5)
		sleepingSecond = sleepingSecond -5	
		if (sleepingSecond <= 0):
			output.kill()
			print "Killing " + bench + " !"
			break

	if (sleepingSecond > 0):
		#The run was correct, we parse the results
		second = 0
		accuracy = 0
		for oneResultLine in output.stdout:

			if (oneResultLine.split()[-1] == "sec"):
				second = int(oneResultLine.split()[-2])
			elif (oneResultLine.split()[0] == "Percent"):
				accuracy = float(oneResultLine.split()[-1][:-1])

		return (second, accuracy)










benchmarks = ["minicoli.fa"]
benchmarks_ref = {"minicoli.fa":"ecoliref.fa"}

for readFile in benchmarks:
	referenceFile = benchmarks_ref[readFile]

	#We run FTL		
	(second, accuracy) = runBench(readFile, referenceFile, 31, 2)
	
	print "Result : " + str(second) + " sec for " + str(accuracy) + "% !"


