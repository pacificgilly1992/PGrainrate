############################################################################
# Project: The Lenard effect of preciptation at the RUAO,
# Title: Ensemble processing of the PG, Time and Rain Rate data,
# Author: James Gilmore,
# Email: james.gilmore@pgr.reading.ac.uk.
# Version: 1.0.8
# Date: 18/01/16
# Status: Operational
############################################################################

#Initialising the python script
from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import glob
import csv
import time
execfile("pgrrdiagnolim.py")

#Import the data files using a search criteria method

#pgfile=glob.glob('RUSOdata/*/*.rax') #UNIX system (Personal)
#pgfilelog=glob.glob('RUSOdata/*/*.rpt') #UNIX system (Personal)

pgfile=glob.glob('../../../../net/vina1s/vol/data1_s/meteorology_2/RUAOData/METFiDAS/data/raw/met/1sec/*/*.rax')
pgfilelog=glob.glob('../../../../net/vina1s/vol/data1_s/meteorology_2/RUAOData/METFiDAS/data/raw/met/1sec/*/*.rpt')

#pgfile=glob.glob('X:/RUAOData/METFiDAS/data/raw/met/1sec/*/*.rax')
#pgfile=glob.glob('C:/Users/james/Documents/Education/PhD/1sec/*/*.rax')

print(len(pgfile))
print(len(pgfilelog))

#Set Initial Conditions
iterator = len(pgfile)
bincount = 100

slope = intercept = r_value = std_err = p_value = np.zeros(iterator)
pearson_cor = np.zeros([iterator,2])
RainRateBin = np.zeros([iterator,bincount])
TimeTipBin = np.zeros([iterator,bincount])
PGtipBin = np.zeros([iterator,bincount])
timekeep = np.zeros([iterator,])
rain = np.zeros([iterator,])
TimeTip5mm = []
#RainRate5mm = np.zeros([iterator, 86400])
#PGtip = np.zeros([iterator, 86400])
#PGdata = np.zeros([iterator,])
ts = np.zeros(iterator)
PGdata = []
it=0

y=yy=n=nn=0

np.set_printoptions(threshold='nan')

tstart = time.time()

#try:
for i in range(iterator):
	it+=1
#		try:
	slope1, intercept1, r_value1, p_value1, std_err1, pearson_cor1, RainRateBin1, TimeTipBin1, PGTipBin1, day1 = PGrainrate(pgfile[i], bincount, pgfilelog[i])
		
	print(i)
	slope[i] = slope1
	intercept[i] = intercept1
	r_value[i] = r_value1
	p_value[i] = p_value1
	std_err[i] = std_err1
	pearson_cor[i] = pearson_cor1
	PGdata.append(day1)		
				
	ts[i] = time.time()
		
	print("Time Left (s): ", "%.0f" % ((iterator-it)*((ts[i]-tstart)/it)))
			
	if p_value[i]!=0:
		yy+=1
	else:	
		RainRateBin1
		nn+=1

	if (p_value[i] < 0.05) & (p_value[i]!=0):
		y+=1
	else:
		n+=1
	
	if np.sum(RainRateBin1) == 0:
		continue
	
	for j in range(bincount):

		TimeTipBin[i,j] = TimeTipBin1[j]
		PGtipBin[i,j] = PGTipBin1[j]
		
#		except:
#			print("########################### Error has on: ", os.path.basename(pgfile[i])[:-4])
#			continue

PGdata =  [val for sublist in PGdata for val in sublist]	

with open("processeddata/PGdata.csv", "wb") as output:
	writer = csv.writer(output, lineterminator='\n')
	writer.writerows(PGdata)

#print("Current Rainday: ", yy)
#print("Rain day ratio: ", yy/nn)
#print("Current Sig: ", y)
#print("Current Sig Ratio: ", y/n)
#print("Current Sig to rainday ratio: ", y/yy)

PGEnsemble(PGtipBin, bincount, iterator)

print("Sig: ", y, " Not Sig: ", n)
	
#except:
#	print("########################### MAJOR Error has on: ", os.path.basename(pgfile[i])[:-4], " calculation terminated. Fingers crossed and hope theres data in PGdatadump.csv")
#	PGdata =  [val for sublist in PGdata for val in sublist]
#	with open("processeddata/PGdatadump.csv", "wb") as output:
#		writer = csv.writer(output, lineterminator='\n')
#		writer.writerows(PGdata)

