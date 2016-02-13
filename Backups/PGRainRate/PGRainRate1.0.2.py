from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import glob
import csv
import time

execfile("pgrrdiag.py")

#Import the data files using a search criteria method

#pgfile=glob.glob('RUSOdata/*/*.rax') UNIX system
pgfile=glob.glob('X:/RUAOData/METFiDAS/data/raw/met/1sec/*/*.rax')
#pgfile=glob.glob('C:/Users/james/Documents/Education/PhD/1sec/*/*.rax')
print(len(pgfile))

slope = intercept = r_value = std_err = np.zeros(len(pgfile))

#Set Initial Conditions
iteratormin = 2400
interatormax = len(pgfile)
bincount = 100

p_value = np.zeros(interatormax)
pearson_cor = np.zeros([len(pgfile),2])
RainRateBin = np.zeros([interatormax,bincount])
TimeTipBin = np.zeros([interatormax,bincount])
PGtipBin = np.zeros([interatormax,bincount])
timekeep = np.zeros([interatormax,])
rain = np.zeros([interatormax,])
TimeTip5mm = []
#RainRate5mm = np.zeros([iterator, 86400])
#PGtip = np.zeros([iterator, 86400])
#PGdata = np.zeros([iterator,])
ts = np.zeros(interatormax)

y=0
yy=0
n=0
nn=0

np.set_printoptions(threshold='nan')


PGdata = []

try:
	for i in range(interatormin, interatormax):
		try:
			slope1, intercept1, r_value1, p_value1, std_err1, pearson_cor1, RainRateBin1, TimeTipBin1, PGTipBin1, day1 = PGrainrate(pgfile[i], bincount)
			print(i)
			slope[i] = slope1
			intercept[i] = intercept1
			r_value[i] = r_value1
			p_value[i] = p_value1
			std_err[i] = std_err1
			pearson_cor[i] = pearson_cor1
			
			PGdata.append(day1)
			
		
		
			ts[i] = time.time()
			print("Time Left: ", (iterator-i)*(ts[i]-ts[i-1]))
			
			if p_value[i]!=0:
				yy+=1
			else:	
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
				
	
			#print(pearson_cor[i])
			#print(RainRate5mm1)
		
		except:
			print("########################### Error has on: ", os.path.basename(pgfile[i])[:-4])
			continue

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
	
except:
	PGdata =  [val for sublist in PGdata for val in sublist]
	with open("processeddata/PGdatadump.csv", "wb") as output:
		writer = csv.writer(output, lineterminator='\n')
		writer.writerows(PGdata)

