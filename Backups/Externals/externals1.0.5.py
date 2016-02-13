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
import matplotlib.pyplot as plt
import numpy as np
from lowess import lowess
import sys

############################################################################
#Plotting functions for PG and Rain Rate

def PGRainFull(xlimmax=None, ylimmax=None, outFile=None, fileformat=None, RainRate5mm=None, TimeTip5mm=None, timekeep=None, PG=None, PGtip=None, slope=None, intercept=None, p_value=None, r_value=None, pearson_cor=None, std_err=None, mann_wht=None):
	"Plot 3 subplots all of which completment the main focus, i.e. (1) PG vs." 
	"Rain Rate along with side plots for (2) Rain Rate and (3) PG between the"
	"times that charged rain was detected. Statistical information was also "
	"added in the remaining quadrant to fill the white space but can easily "
	"be removed if neseccary."
	
	plt.clf()
	fig = plt.figure()
	#plt.suptitle("Raindrop Charge: " + outFile)

	pgrain = fig.add_subplot(222)
	pgrain.scatter(RainRate5mm, PGtip)
	pgrain.set_xlabel("Rain Rate (mm/hr)")
	pgrain.set_ylabel("Potential Gradient (V/m)")
	pgrain.grid()
	pgrain.set_xlim(-.1,xlimmax)
	pgrain.set_ylim(-1050, ylimmax)
	pgrain.invert_yaxis()
	
	pgrain.plot(np.arange(-.1, xlimmax+0.3, 0.2),np.arange(-.1, xlimmax+0.3, 0.2)*slope+intercept)

	PGRainsort = np.array(sorted(zip(RainRate5mm, PGtip)))
 
	eps = sys.float_info.epsilon
	pgrain.plot(PGRainsort[:,0], lowess(PGRainsort[:,0]+eps, PGRainsort[:,1]+eps, 1/2))

	x0, x1 = pgrain.get_xlim()
	y0, y1 = pgrain.get_ylim()
	pgrain.set_aspect(np.abs((x1-x0)/(y1-y0)))
	
	#PG Plot

	pg = fig.add_subplot(221)
	pg.plot(timekeep,PG)
	pg.set_xlabel("Time (hrs)")
	pg.set_xlim(np.min(TimeTip5mm),np.max(TimeTip5mm))
	pg.set_ylim(-1050, ylimmax)
	pg.invert_yaxis()
	#pg.axes.get_yaxis().set_visible(False)
	pg.grid()
	
	x0, x1 = pg.get_xlim()
	y0, y1 = pg.get_ylim()
	pg.set_aspect(np.abs((x1-x0)/(y1-y0)))

	#Rain plot

	rain = fig.add_subplot(224)
	rain.plot(RainRate5mm,TimeTip5mm)
	rain.set_ylabel("Time (hrs)")
	rain.set_ylim(np.min(TimeTip5mm),np.max(TimeTip5mm))
	rain.set_xlim(-.1,xlimmax)
	rain.grid()

	x0, x1 = rain.get_xlim()
	y0, y1 = rain.get_ylim()
	rain.set_aspect(np.abs((x1-x0)/(y1-y0)))

	#Info Plot

	info = fig.add_subplot(223)
	info.axis('off')
	info.text(-0.1, .9, '$Year and Day$', fontsize=15)
	info.text(-0.1, .75, '$P-Value$: ', fontsize=15)
	info.text(-0.1, .6, '$R^2$: ', fontsize=15)
	info.text(-0.1, .45, "$Pearson's Cor$: ", fontsize=15)
	info.text(-0.1, .3, "$Standard Error$: ", fontsize=15)
	info.text(-0.1, .15, "$Mann-Whitney$: ", fontsize=15)

	info.text(0.6, .9, outFile, fontsize=15)
	info.text(0.6, .75, round(p_value,7), fontsize=15)
	info.text(0.6, .6, round(r_value**2,5), fontsize=15)
	info.text(0.6, .45, round(pearson_cor[1],5), fontsize=15)
	info.text(0.6, .3, round(std_err,5), fontsize=15)
	info.text(0.6, .15, round(mann_wht,5), fontsize=15)
	

	x0, x1 = info.get_xlim()
	y0, y1 = info.get_ylim()
	info.set_aspect(np.abs((x1-x0)/(y1-y0)))

	plt.tight_layout(pad=0.4, w_pad=-0.5, h_pad=0.5)


	plt.savefig('plots/new/' + outFile + "." + fileformat)
	plt.close(fig)

def PGRainSlim(xlimmax=None, ylimmax=None, outFile=None, fileformat=None, RainRate5mm=None, PGtip=None, slope=None, intercept=None):

	plt.clf()
	fig = plt.figure()
	#plt.suptitle("Raindrop Charge: " + outFile)

	pgrain = fig.add_subplot(111)
	pgrain.scatter(RainRate5mm, PGtip)
	pgrain.set_xlabel("Rain Rate (mm/hr)")
	pgrain.set_ylabel("Potential Gradient (V/m)")
	pgrain.grid()
	pgrain.set_xlim(-.1,xlimmax)
	pgrain.set_ylim(-200, ylimmax)
	pgrain.invert_yaxis()
	
	#pgrain.plot(np.arange(-.1, xlimmax+0.3, 0.2),np.arange(-.1, xlimmax+0.3, 0.2)*slope+intercept)

	PGRainsort = np.array(sorted(zip(RainRate5mm, PGtip)))
 
	eps = sys.float_info.epsilon
	pgrain.plot(PGRainsort[:,0], lowess(PGRainsort[:,0]+eps, PGRainsort[:,1]+eps, 1/2))

	x0, x1 = pgrain.get_xlim()
	y0, y1 = pgrain.get_ylim()
	pgrain.set_aspect(np.abs((x1-x0)/(y1-y0)))	

	plt.savefig('plots/new/' + outFile + "." + fileformat)
	plt.close(fig)
