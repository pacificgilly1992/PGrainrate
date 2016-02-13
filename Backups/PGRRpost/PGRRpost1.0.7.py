############################################################################
# Project: The Lenard effect of preciptation at the RUAO,
# Title: Ensemble processing of the PG, Time and Rain Rate data,
# Author: James Gilmore,
# Email: james.gilmore@pgr.reading.ac.uk.
# Version: 1.0.7
# Date: 18/01/16
# Status: Operational (Basic)
############################################################################

#Initialising the python script
from __future__ import absolute_import, division, print_function
from scipy import stats, interpolate
from lowess import lowess
from array import array
import sys
import numpy as np
import matplotlib.pyplot as plt
execfile("externals.py")
np.set_printoptions(threshold='nan')

y = "y"
n = "n"

#User input for the further processing of the PGRR data
print("####################################################################")
print("The Lenard effect of preciptation at the RUAO. Using the processed ")
print("data collected from the PGRainRate.py script the average for each ")
print("rain rate can be found.")
print("####################################################################\n")
selectcase = input("Please select the averaging method: Type '1' for Mean, Type '2' for Median: ")
loop = str(input('Do you want to loop over many bins? y/n: '))
if loop == "n":
	bincount = input("How many bins for the averaging would you like (recommended = 100): ")
	loop = 1
elif loop == "y":
	bincount = 30
	loop = 10


#Import the processed data for the significantly charged rain. See PGRainRate.py
year, month, time, rainrate, pg = np.genfromtxt('processeddata/PGdata.csv', dtype=float, delimiter=',',  unpack=True)

#Remove zero values from processed data
Month = month.copy()[year.copy()!=0]
Time = time.copy()[year.copy()!=0]
Rainrate = rainrate.copy()[year.copy()!=0]
PG = pg.copy()[year.copy()!=0]
Year = year.copy()[year.copy()!=0]

slope = intercept = r_value = p_value = std_err =  np.zeros(int(loop))
lowessval = np.zeros([int(loop),bincount+30*loop])

PGRR = np.asarray(zip(Rainrate, PG))

PGRRsort = PGRR[np.lexsort((PGRR[:, 1], PGRR[:, 0]))]

PGsort = PGRRsort[:,1]
RRsort = PGRRsort[:,0]

for k in xrange(loop):
	#Initalise the matrices and vectors
	RainRateBin = np.zeros((bincount+30*k)-1)
	RainRateBinLimit = np.zeros(bincount+30*(k-1))
	TimeTipBin = np.zeros(bincount+30*(k-1))
	PGTipBin = np.zeros(bincount+30*(k-1))
	TotalBin = np.zeros(bincount+30*(k-1))
	PGTipBinMedian = np.zeros([bincount+30*(k-1),len(Year)])
	PGTipPosition = np.zeros(bincount+30*(k-1))
	PGTipBinMedianFinal = np.zeros(bincount+30*(k-1))
	eps = sys.float_info.epsilon

	#Define the Rain Rate for each bin with the centred values determined as well.
	for i in range(bincount+30*(k-1)):
		RainRateBinLimit[i] = i*5/(bincount+30*(k-1))
	for i in range((bincount+30*(k-1))-1):
		RainRateBin[i] = 0.5*(RainRateBinLimit[i+1]-RainRateBinLimit[i])

	if selectcase == 1:
		############################################################################
		#Define the mean (ensemble) PG and Tip Times for the statistically significant data.
	
		for j in range(len(Year)):
			for i in range(1,bincount+30*(k-1)):
				if (Rainrate[j] < i*5/(bincount+30*(k-1)) and Rainrate[j] > (i-1)*5/(bincount+30*(k-1))):
					PGTipBin[i] += PG[j]
					TimeTipBin[i] += Time[j]
					TotalBin[i] += 1
		PGTipBinned = PGTipBin.copy()/(TotalBin.copy())
		TimeTipBinned = TimeTipBin.copy()/(TotalBin.copy())

		#Removes NaN values
		PGTipBinned = [0 if np.isnan(x) else x for x in PGTipBinned]
		TimeTipBinned = [0 if np.isnan(x) else x for x in TimeTipBinned]
	
		#Select values for plotting
		yvalue = PGTipBinned
		amethod = "Mean"
	
		############################################################################

	elif selectcase == 2:	
		############################################################################
		#Define the median PG and Tip Times for the statistically significant data.

		for j in range(len(Year)):
			for i in range(bincount+30*(k-1)):
				if (Rainrate[j] < i*5/(bincount+30*(k-1)) and Rainrate[j] > (i-1)*5/(bincount+30*(k-1))):
					PGTipBinMedian[i,PGTipPosition[i]] = PG[j]
					PGTipPosition[i]+=1

		for i in range(bincount+30*(k-1)):
			PGTipBinMedianFinal[i] = np.median(PGTipBinMedian[i,:].copy()[PGTipBinMedian[i,:].copy()!=0])		
			PGTipBinMedianFinal[np.isnan(PGTipBinMedianFinal)] = 0
		
		#Select values for plotting
		yvalue = PGTipBinMedianFinal
		amethod = "Median"

		############################################################################

	else:
		sys.exit("Please select either the Mean (1) or Median (2) case.")

	#Calculation of the linear regression model along with statistical parameters.
	slope[k], intercept[k], r_value[k], p_value[k], std_err[k] = stats.linregress(RainRateBinLimit, yvalue)
	for m in xrange(len(lowess(RainRateBinLimit+eps, yvalue+eps, 1/2))):
		lowessval[k,m] = lowess(RainRateBinLimit+eps, yvalue+eps, 1/2)[m]

if loop == 10:
	PGRainEnsembleMulti(np.max(RainRateBinLimit)+0.2, np.max(yvalue)+0.2, "PGEnsembleMulti" + str(amethod) + str(bincount), "png", RainRateBinLimit, yvalue, lowessval)
	
print(slope, intercept, r_value, p_value, std_err)

print(lowessval)

print("P-Value: ", p_value)
print("R^2 Value: ", r_value**2)
print("Standard Error: ", std_err)

#Plot the ensemble PG against Rain Rate. See external.py for the source function.
#PGRainSlim(np.max(RainRateBinLimit)+0.2, np.max(yvalue)+0.2, "PGEnsemble" + str(amethod) + str(bincount), "png", RainRateBinLimit, yvalue, slope, intercept)
