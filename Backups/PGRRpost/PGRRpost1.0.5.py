############################################################################
# Project: The Lenard effect of preciptation at the RUAO,
# Title: Ensemble processing of the PG, Time and Rain Rate data,
# Author: James Gilmore,
# Email: james.gilmore@pgr.reading.ac.uk.
# Version: 1.0.5
# Date: 07/12/15
############################################################################

#Initialising the python script
from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate

execfile("externals.py")

np.set_printoptions(threshold='nan')

year, month, time, rainrate, pg = np.genfromtxt('processeddata/PGdata.csv', dtype=float, delimiter=',',  unpack=True)

Month = month.copy()[year.copy()!=0]
Time = time.copy()[year.copy()!=0]
Rainrate = rainrate.copy()[year.copy()!=0]
PG = pg.copy()[year.copy()!=0]
Year = year.copy()[year.copy()!=0]

PGRR = np.asarray(zip(Rainrate, PG))

PGRRsort = PGRR[np.lexsort((PGRR[:, 1], PGRR[:, 0]))]

PGsort = PGRRsort[:,1]
RRsort = PGRRsort[:,0]

bincount = 100

RainRateBin = np.zeros(bincount-1)
RainRateBinLimit = np.zeros(bincount)
TimeTipBin = np.zeros(bincount)
PGTipBin = np.zeros(bincount)
TotalBin = np.zeros(bincount)
PGTipBinMedian = np.zeros([bincount,150])
PGTipPosition = np.zeros(bincount)
PGTipBinMedianFinal = np.zeros(bincount)

#Define the Rain Rate for each bin with the centred values determined as well.
for i in range(bincount):
	RainRateBinLimit[i] = i*5/bincount
for i in range(bincount-1):
	RainRateBin[i] = 0.5*(RainRateBinLimit[i+1]-RainRateBinLimit[i])


############################################################################
#Define the mean (ensemble) PG and Tip Times for the statistically significant data.
for j in range(len(Year)):
	for i in range(1,bincount):
		if (Rainrate[j] < i*5/bincount and Rainrate[j] > (i-1)*5/bincount):
			PGTipBin[i] += PG[j]
			TimeTipBin[i] += Time[j]
			TotalBin[i] += 1	
PGTipBinned = PGTipBin.copy()/(TotalBin.copy())
TimeTipBinned = TimeTipBin.copy()/(TotalBin.copy())

#Removes NaN values
PGTipBinned = [0 if np.isnan(x) else x for x in PGTipBinned]
TimeTipBinned = [0 if np.isnan(x) else x for x in TimeTipBinned]

############################################################################

############################################################################
#Define the median PG and Tip Times for the statistically significant data.

print(PGTipBinMedian[i,:] )

for j in range(len(Year)):
	for i in range(bincount):
		if (Rainrate[j] < i*5/bincount and Rainrate[j] > (i-1)*5/bincount):
			PGTipBinMedian[i,PGTipPosition[i]] = PG[j]
			PGTipPosition[i]+=1

print("HHHHHHHHHHHH", PGTipBinMedian[3,:].copy()[PGTipBinMedian[3,:].copy()!=0])
			
for i in range(bincount):
	PGTipBinMedian[i:,]= PGTipBinMedian[i,:][PGTipBinMedian[i,:]!=0]
	PGTipBinMedianFinal[i] = np.median(PGTipBinMedian[i:,])		

print(PGTipBinMedianFinal)
print(PGTipPosition)

############################################################################

#print(TotalBin)
#print(PGTipBinned)

#Calculation of the linear regression model along with statistical parameters.
#slope, intercept, r_value, p_value, std_err = stats.linregress(RainRateBin, PGTipBinned)

#print("P-Value: ", p_value)
#print("R^2 Value: ", r_value**2)
#print("Standard Error: ", std_err)

#Plot the ensemble PG against Rain Rate. See external.py for the source function.
PGRainSlim(np.max(RainRateBin)+0.2, np.max(PGTipBinned)+0.2, "PGEnsemble" + str(bincount), "png", RainRateBin, PGTipBinned)
