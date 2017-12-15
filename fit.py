import json
import numpy
import codecs
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#def myfunction(xdata,a,b):
#	return b*numpy.sqrt((xdata**2/a**2)-1)

#aug_high=[1981,1982,1986,1989,1998,2003,2014]
#aug_low =[1979,1984,1985,1995]
#sep_high=[1982,1988,1992,1998,2002,2006,2014]
#sep_low =[1984,1985,1995,2001,2009,2012]
#patrik_high=[1994,1996,2001,2013,2014]
#patrik_low=[1981,1989,1993,2003,2007,2008,2012]


xdata = numpy.arange(1979,2015)
data = codecs.open('patrickEOF1_1979_2014','r').read()
data2 = json.loads(data)
data3 = numpy.array(data2)
#xdata = numpy.arange(10,46)
ydata = data3
zdata = numpy.polyfit(xdata,ydata,2)
poly  = numpy.poly1d(zdata)
#popt, pcov = curve_fit(myfunction,xdata,ydata)
#print myfunction(xdata,*popt)
residu = ydata - poly(xdata)
#residu = ydata - myfunction(xdata,*popt)
std = numpy.repeat(numpy.std(residu),xdata.size)
average = numpy.repeat(numpy.average(residu),xdata.size)
stdUp = average + std
stdDown = average - std

##################################
'''
ydata = codecs.open('DivQIndex_1979_2014_Sep_10W40E_35N75N','r').read()
#ydata = codecs.open('DivDIndex_1979_2014_Sep_10W40E_35N75N','r').read()
ydata2 = json.loads(ydata)
ydata3 = numpy.array(ydata2)
ydata = ydata3
zdata = numpy.polyfit(xdata,ydata,1)
poly  = numpy.poly1d(zdata)
residu = ydata - poly(xdata)
#residu = ydata - myfunction(xdata,*popt)
std = numpy.repeat(numpy.std(residu),xdata.size)
#average = numpy.average(ydata)
residu = ydata - poly(xdata)
average = numpy.repeat(numpy.average(residu),xdata.size)
std = numpy.repeat(numpy.std(residu),xdata.size)
stdUp = average + std
stdDown = average - std
'''
##################################

fig, ax = plt.subplots(1,figsize = (13,6))
fitData = ax.plot(xdata,ydata,'-b',linewidth = 3, label = 'data')
fitPlot = ax.plot(xdata,poly(xdata),'-r',linewidth = 3, label = 'fit')
resiplot = ax.plot(xdata,residu,'-g',linewidth = 3, label = 'residual')
stdUpPlot = ax.plot(xdata,stdUp,'--m',linewidth = 3, label = 'std')
stdDownPlot = ax.plot(xdata,stdDown,'--m',linewidth = 3)

plt.xlabel('Year')
#plt.ylabel('Latent energy index (Wm-2)')
plt.ylabel('Ice Index')
plt.title("Aug ice index anomalies (1979,2014)")
plt.grid(True)
plt.legend(loc=1)
plt.xlim(1979,2014)

#plt.savefig('/home/nba035/plot/iceIndexAnomalyFitSepEOF1_79_14.eps', format = 'eps')
plt.show()


