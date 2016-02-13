from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import glob
import csv
from itertools import izip_longest
import time

execfile("pgrrdiag.py")

#Import the data files using a search criteria method

#pgfile=glob.glob('RUSOdata/2009/*.rax') UNIX system
pgfile=glob.glob('Z:/RUAOData/METFiDAS/data/raw/met/1sec/2009/*.rax')
print(len(pgfile))

slope = intercept = r_value = std_err = np.zeros(len(pgfile))

#Set Initial Conditions
iterator = len(pgfile)
bincount = 100

p_value = np.zeros(iterator)
pearson_cor = np.zeros([len(pgfile),2])
RainRateBin = np.zeros([iterator,bincount])
TimeTipBin = np.zeros([iterator,bincount])
PGtipBin = np.zeros([iterator,bincount])
timekeep = np.zeros([iterator,])
rain = np.zeros([iterator,])
TimeTip5mm = []
RainRate5mm = np.zeros([iterator, 86400])
PGtip = np.zeros([iterator, 86400])
PGdata = np.zeros([iterator,])
ts = np.zeros(iterator)

y=0
yy=0
n=0
nn=0

np.set_printoptions(threshold='nan')


for i in range(iterator):
	slope1, intercept1, r_value1, p_value1, std_err1, pearson_cor1, RainRateBin1, TimeTipBin1, PGTipBin1, day1 = PGrainrate(pgfile[i], bincount)
	print(i)
	slope[i] = slope1
	intercept[i] = intercept1
	r_value[i] = r_value1
	p_value[i] = p_value1
	std_err[i] = std_err1
	pearson_cor[i] = pearson_cor1
	
	ts[i] = time.time()
	print("Time Left: ", (iterator-i)*i/(ts[i]-ts[0]))
	
	if np.sum(RainRateBin1) == 0:
		continue
	
	for j in range(bincount):
	
		TimeTipBin[i,j] = TimeTipBin1[j]
		PGtipBin[i,j] = PGTipBin1[j]
		
	if p_value[i]!=0:
		yy+=1
	else:
		nn+=1

	if (p_value[i] < 0.05) & (p_value[i]!=0):
		y+=1
	else:
		n+=1	
	
	#print(pearson_cor[i])
	#print(RainRate5mm1)

#print("Current Rainday: ", yy)
#print("Rain day ratio: ", yy/nn)
#print("Current Sig: ", y)
#print("Current Sig Ratio: ", y/n)
#print("Current Sig to rainday ratio: ", y/yy)

PGEnsemble(PGtipBin, bincount, iterator)
	
print("Sig: ", y, " Not Sig: ", n)


