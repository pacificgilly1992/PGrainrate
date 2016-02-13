from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import glob

execfile("pgrrdiag.py")

#pgfile=glob.glob('RUSOdata/2009/*.rax') UNIX system
pgfile=glob.glob('Z:/RUAOData/METFiDAS/data/raw/met/1sec/2009/*.rax')
print(len(pgfile))

slope = intercept = r_value = std_err = np.zeros(len(pgfile))

p_value = np.zeros(len(pgfile))
pearson_cor = np.zeros([len(pgfile),2])
RainRateBin = np.zeros([len(pgfile),100])
TimeTipBin = np.zeros([len(pgfile),100])
PGtipBin = np.zeros([20,100])
timekeep = np.zeros([len(pgfile),])
rain = np.zeros([len(pgfile),])

y=0
yy=0
n=0
nn=0

bincount=100

for i in range(16,200):
	slope1, intercept1, r_value1, p_value1, std_err1, pearson_cor1, RainRateBin1, TimeTipBin1, PGTipBin1, timekeep1 = PGrainrate(pgfile[i], bincount)
	print(i)
	slope[i] = slope1
	intercept[i] = intercept1
	r_value[i] = r_value1
	p_value[i] = p_value1
	std_err[i] = std_err1
	pearson_cor[i] = pearson_cor1
	
	if np.sum(RainRateBin1) == 0:
		continue
	
	for j in range(bincount):
	
		RainRateBin[i,j] = RainRateBin1[j]
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

PGtotal = np.zeros(bincount)
total = np.zeros(bincount)

print("bincount: ", bincount)
print("totalbefore: ", total)

PGtotal = [sum(row[i] for row in PGtipBin) for i in range(len(PGtipBin[0]))]

print(PGtotal)

for j in xrange(bincount):
	for i in xrange(len(pgfile)):
		if PGtipBin[i, j] > 0:
			total[i] += 1	

print("total: ", total)
	
print("Sig: ", y, " Not Sig: ", n)


