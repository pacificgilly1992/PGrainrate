from __future__ import absolute_import, division, print_function
from array import array
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate
import sys, os

execfile("externals.py")

def PGrainrate(file="RUSOdata/raw/2009M062.rax", bin=100):

	time, rain, pg = np.genfromtxt(file, dtype=float, delimiter=',', usecols=(3, 26, 32),  unpack=True, skiprows=2)
 
	timekeep = time.copy()

	#Calibrate the PG data from V to V/m

	PG=(pg.copy()-0.00903)/0.00463

	PGTime = np.asarray(zip(timekeep, PG))


	#Calibrate the rain data from V to mm

	Rain = np.zeros_like(time)

	time[0]=time[-1]=0	# Set the first and last value of the time series equal to 0	
	raintot=0   		# Number of tip buckets has occurred in a day

	#Buffer tip times

	for i in range(len(time)-1):
		if rain[i+1]-rain[i]<=0.03:
			time[i]=0
			Rain[i]=0
		else:
			Rain[i]=0.2
			raintot+=1

	#Remove traces when tips has not occured

	timerain = time.copy()[time.copy()!=0]
	Rain = Rain[Rain!=0]

	if len(timerain) < 5:
		return 0,0,0,0,0,0,0,0,0,0

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
		if not RainRate[i]<=5:
			TimeTip[i]=-1

	RainRate5mm = RainRate[(RainRate<=5) & (RainRate>0)]
	TimeTip5mm = TimeTip[TimeTip>0]
	
	if len(RainRate5mm) < 1:
		return 0,0,0,0,0,0,0,0,0,0

	PGtip = np.zeros_like(TimeTip5mm)

	for i in range(len(TimeTip5mm)-1):

		PGtemp = PG[(TimeTip5mm[i] <= timekeep) & (timekeep < TimeTip5mm[i+1])]
		PGtip[i]=np.median(PGtemp)

	#Fit a linear regression model for the PG against rain rate and return some statistics
	"This is assumed at the moment that the function is linearly dependent"
	"where more rigorous testing of the rain rate data needs to be done"
	
	eps = sys.float_info.epsilon
	RRL = RainRate5mm[RainRate5mm<0.5*np.max(RainRate5mm)]/(PGtip[RainRate5mm<0.5*np.max(RainRate5mm)]+eps)
	RRH = RainRate5mm[RainRate5mm>=0.5*np.max(RainRate5mm)]/(PGtip[RainRate5mm>=0.5*np.max(RainRate5mm)]+eps)

	statistic, pvalue = stats.mannwhitneyu(RRL, RRH)

	print("Mann: ",pvalue)

	slope, intercept, r_value, p_value, std_err = stats.linregress(RainRate5mm, PGtip)
	pearson_cor = stats.pearsonr(TimeTip5mm,RainRate5mm)

	#Print the results in a quadrant form displaying the variations in both variables

	
	if p_value >= 0.05:
		return 0,0,0,0,0,0,0,0,0,0
		
	PGRainFull(np.max(RainRate5mm)+0.2, np.max(PGtip)+100, os.path.basename(file)[:-4], "png" , RainRate5mm, TimeTip5mm, timekeep, PG, PGtip, slope, intercept, p_value, r_value, pearson_cor, std_err, pvalue)
	
	RainRateBin = np.zeros(bincount)
	TimeTipBin = np.zeros(bincount)
	PGTipBin = np.zeros(bincount)
	
	for i in range(bincount):
		RainRateBin[i] = i*5/bincount
		
	for j in range(len(RainRate5mm)):
		for i in range(1,bincount):
			if (RainRate5mm[j] < i*5/bincount and RainRate5mm[j] > (i-1)*5/bincount):
				PGTipBin[i] = PGtip[j]
				TimeTipBin[i] = TimeTip5mm[j]
	
	#Return from the definition with a barrage of data

	return slope, intercept, r_value, p_value, std_err, pearson_cor, RainRateBin, TimeTipBin, PGTipBin, timekeep
