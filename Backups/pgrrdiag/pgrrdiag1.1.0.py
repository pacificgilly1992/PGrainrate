############################################################################
# Project: The Lenard effect of preciptation at the RUAO,
# Title: Ensemble processing of the PG, Time and Rain Rate data,
# Author: James Gilmore,
# Email: james.gilmore@pgr.reading.ac.uk.
# Version: 1.1.0
# Date: 10/02/16
# Status: Operational
# Change: Added all the error checking procedures! We've got it ALL cover now! Happy Days!
############################################################################

#Initialising the python script
from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate
import sys, os
import re

execfile("externals.py")

def PGrainrate(file="RUSOdata/raw/2009M062.rax", bin=100, filelog=None):

	#Load Data and Error Log files
	while True:
		try:
			time, rain, pg = np.genfromtxt(file, dtype=float, delimiter=',', usecols=(3, 26, 32),  unpack=True, skiprows=1)
			end=1
			break
		except:
			time, rain, pg = np.genfromtxt(file, dtype=float, delimiter=',', usecols=(3, 26, 32),  unpack=True, skiprows=2)
			end=2
			break
	while True:
		try:
			log = np.genfromtxt(filelog, dtype=str, delimiter=',',  unpack=True)
			break
	
		except:
			log=0
			break
	cal_year,cal_day,cal_Tipoffset,cal_Tipmultipler,cal_PGoffset,cal_PGmultipler = np.genfromtxt('RUAOdata/PGRRCalibration.ini', dtype=float, delimiter=',',  unpack=True, skiprows=1)

	######################################################################
	#Quality Control Section
	#PART 1: Fixes all the clock skips and jumps recoreded in the log files
	
	if log!=0: #for any day without an associated log file this will stop errors from cropping up and removing the particular data from processing.
		logged = np.zeros([len(log)-end,7])
		for i in xrange(len(log)-end):
			if log[i][:15] != "*** WARNING ***":
				bob = [int(s) for s in re.findall(r'[-\d]+', log[i])]
				logged[i,6]=len(log[i])
				for j in xrange(len(bob)):
					logged[i,j]=bob[j]
					bob[j]=0

		for i in xrange(len(logged)-1):
			time[logged[i,0]-end] += logged[i,1]*(1/3600)

	#PART 2: Find the correct calibration values for the data

	for i in xrange(len(cal_year)-1):
		if (float(os.path.basename(file)[:-8]+"."+os.path.basename(file)[5:-4])>(cal_year[i]+cal_day[i]/1000) and float(os.path.basename(file)[:-8]+"."+os.path.basename(file)[5:-4])<(cal_year[i+1]+cal_day[i+1]/1000)):
			cal_no = i
	
	######################################################################

	timekeep = time.copy()
	
	#Calibrate the PG data from V to V/m with correct values (hence cal_no)

	PG=(pg.copy()-cal_PGoffset[cal_no])/cal_PGmultipler[cal_no]
	
	PGTime = np.asarray(zip(timekeep, PG))


	#Calibrate the rain data from V to mm

	Rain = np.zeros_like(time)

	#time[0]=time[-1]=0	# Set the first and last value of the time series equal to 0	
	raintot=0   		# Number of tip buckets has occurred in a day

	#Buffer tip times

	for i in range(len(time)-1):
		if rain[i+1]-rain[i]<(0.2*cal_Tipmultipler[cal_no]):
			time[i]=0
			Rain[i]=0
		else:
			Rain[i]=0.2
			raintot+=1

	#Remove traces when tips has not occured

	timerain = time.copy()[time.copy()!=0]
	Rain = Rain[Rain!=0]

	if len(timerain) < 6:
		return 0,0,0,0,0,0,0,0,0,[(0,0,0,0,0)]

	RainRate = np.zeros_like(timerain)
	TimeTip = np.zeros_like(timerain)
	
	#Determine the rainfall rate (each bucket tip is 0.2mm)
	"Note that as the time is currently in units of hours that the"
	"derivatives will automatically be displayed in per unit hours"

	for i in range(len(timerain)-1):
		if timerain[i+1]>timerain[i]:
			RainRate[i]=0.2/(timerain[i+1]-timerain[i])
			TimeTip[i]=0.5*(timerain[i+1]+timerain[i])	
			
	#Determine rainfall rates less than 5mm/hr and relate to relevant time

	for i in range(len(RainRate)):
		if not RainRate[i]<=50:
			TimeTip[i]=-1

	RainRate5mm = RainRate[RainRate>0]
	TimeTip5mm = TimeTip[TimeTip>0]
	
	if len(RainRate5mm) < 1:
		return 0,0,0,0,0,0,0,0,0,[(0,0,0,0,0)]

	PGtip = np.zeros_like(TimeTip5mm)

	for i in range(len(timerain)-1):

		PGtemp = PG[(timerain[i] <= timekeep) & (timekeep < timerain[i+1])]
		PGtip[i] = np.median(PGtemp)

	#Fit a linear regression model for the PG against rain rate and return some statistics
	"This is assumed at the moment that the function is linearly dependent"
	"where more rigorous testing of the rain rate data needs to be done"
	
	eps = sys.float_info.epsilon
	#RRL = RainRate5mm[RainRate5mm<0.5*np.max(RainRate5mm)]/(PGtip[RainRate5mm<0.5*np.max(RainRate5mm)]+eps)
	#RRH = RainRate5mm[RainRate5mm>=0.5*np.max(RainRate5mm)]/(PGtip[RainRate5mm>=0.5*np.max(RainRate5mm)]+eps)

	#statistic, pvalue = stats.mannwhitneyu(RRL, RRH)
	pvalue=1
	#print("Mann: ",pvalue)

	slope, intercept, r_value, p_value, std_err = stats.linregress(RainRate5mm, PGtip)
	pearson_cor = stats.pearsonr(TimeTip5mm,RainRate5mm)

	#Print the results in a quadrant form displaying the variations in both variables

################REMOVED PVALUE CRITERIA###############	
	#if p_value >= 0.05:
	#	return 0,0,0,0,0,0,0,0,0,[(0,0,0,0,0)]
	
	#Plot the Rate Rate against PG once passed all basic quality checks
	print("Printing :)")	
	PGRainFull(np.max(RainRate5mm)+0.2, np.max(PGtip)+100, os.path.basename(file)[:-4], "png" , RainRate5mm, TimeTip5mm, timekeep, PG, PGtip, slope, intercept, p_value, r_value, pearson_cor, std_err, pvalue)
	
	#Now place each value of the PG into specific Rain rate ranges or 'bins'. This is
	#a form of discretisation as we are moving from the continuous to discrete regime.
	
	RainRateBin = np.zeros(bincount)
	TimeTipBin = np.zeros(bincount)
	PGTipBin = np.zeros(bincount)
	TotalBin = np.zeros(bincount)
	
	for i in range(bincount):
		RainRateBin[i] = i*5/bincount
		
	for j in range(len(RainRate5mm)):
		for i in range(1,bincount):
			if (RainRate5mm[j] < i*5/bincount and RainRate5mm[j] > (i-1)*5/bincount):
				PGTipBin[i] += PGtip[j]
				TimeTipBin[i] += TimeTip5mm[j]
				TotalBin[i] += 1

	PGTipBinned = PGTipBin.copy()/(TotalBin.copy()+eps)
	TimeTipBinned = TimeTipBin.copy()/(TotalBin.copy()+eps)
	
	#Removes NaN values
	PGTipBinned = [0 if np.isnan(x) else x for x in PGTipBinned]
	TimeTipBinned = [0 if np.isnan(x) else x for x in TimeTipBinned]	

	Year = np.zeros_like(TimeTip5mm)
	Day = np.zeros_like(TimeTip5mm)
	
	for i in xrange(len(Year)):
		Year[i] = os.path.basename(file)[:-8]
		Day[i] = os.path.basename(file)[5:-4]
	
	day = zip(Year, Day, TimeTip5mm, RainRate5mm, PGtip, timerain)
	#Return from the definition with a barrage of data

	return slope, intercept, r_value, p_value, std_err, pearson_cor, RainRateBin, TimeTipBinned, PGTipBin, day

def PGEnsemble(PGtipBin=None, bincount=None, iterator=None):
	"Creates an Ensemble of PG values for every statistically significant"
	"day that we have experience charged rain droplets via space charge."
	"Determines the PG in arbitrary set bins and thus averages the PG"
	"in each bin with the number of events that has happened in said bin."

	PGtotal = np.zeros(bincount)
	total = np.zeros(bincount)
	PGbinMean = np.zeros(bincount)
	RainRateBin = np.zeros(bincount)
	total = np.zeros_like(PGtipBin)

	#Print some information about the analysis
	print("bincount: ", bincount)

	#Separate out the PG data into regular time intervals 
	PGtotal = np.asarray([sum(row[i] for row in PGtipBin) for i in range(len(PGtipBin[0]))])

	#Convert PGtipBin into binary for when PG was recorded
	for i in xrange(iterator): #Usually the length of pgfile 
		for j in xrange(bincount):
			if PGtipBin[i, j] != 0:
				total[i, j] = 1
			
	#Then sum how many events happened on each bin			
	PGtotalcounts = np.asarray([sum(row[i] for row in total) for i in range(len(total[0]))])

	for i in range(bincount):
			RainRateBin[i] = i*5/bincount

	#Find the average PG for each bin and delete corresponding value in Rainratebin
	for i in xrange(bincount):
		if PGtotalcounts[i] != 0:
			PGbinMean[i] = PGtotal[i]/PGtotalcounts[i]
		else:
			RainRateBin[i] = 0

	#Remove zero entries from the list

	PGbinMean = PGbinMean[PGbinMean!=0]
	RainRateBin = RainRateBin[RainRateBin!=0]
	
	slope, intercept, r_value, p_value, std_err = stats.linregress(RainRateBin, PGbinMean)
	
	plt.scatter(RainRateBin, PGbinMean)
	plt.show()
	
	print("P-Value: ", p_value)
	print("R^2 Value: ", r_value**2)
	print("Standard Error: ", std_err)
	
	PGRainSlim(np.max(RainRateBin)+0.2, np.max(PGbinMean)+0.2, "PGEnsemble", "png", RainRateBin, PGbinMean, slope, intercept)
	
	return 
