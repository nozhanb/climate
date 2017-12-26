import json
import numpy
import codecs
import matplotlib
import matplotlib.pyplot as plt

month = 'Sep'
highYear1 = str(1982)
highYear2 = str(1988)
lowYear1 = str(1984)
lowYear2 = str(1985)

'''
highPlot = ['DivQIndexLand_'+highY1+'_'+highY2+'_Mar_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Apr_10W40E_35N75N',\
'DivQIndexLand_'+highY1+'_'+highY2+'_May_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Jun_10W40E_35N75N',\
'DivQIndexLand_'+highY1+'_'+highY2+'_Jul_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Aug_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Sep_10W40E_35N75N',\
'DivQIndexLand_'+highY1+'_'+highY2+'_Oct_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Nov_10W40E_35N75N','DivQIndexLand_'+highY1+'_'+highY2+'_Dec_10W40E_35N75N',\
'DivQIndexLand_'+highY1+'_'+highY2+'_Jan_10W40E_35N75N']

avePlot = ['DivQIndexLand_1979_2014_Mar_10W40E_35N75N','DivQIndexLand_1979_2014_Apr_10W40E_35N75N',\
'DivQIndexLand_1979_2014_May_10W40E_35N75N','DivQIndexLand_1979_2014_Jun_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jul_10W40E_35N75N','DivQIndexLand_1979_2014_Aug_10W40E_35N75N','DivQIndexLand_1979_2014_Sep_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Oct_10W40E_35N75N','DivQIndexLand_1979_2014_Nov_10W40E_35N75N','DivQIndexLand_1979_2014_Dec_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jan_10W40E_35N75N']

lowPlot = ['DivQIndexLand_'+lowY1+'_'+lowY2+'_Mar_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Apr_10W40E_35N75N',\
'DivQIndexLand_'+lowY1+'_'+lowY2+'_May_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Jun_10W40E_35N75N',\
'DivQIndexLand_'+lowY1+'_'+lowY2+'_Jul_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Aug_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Sep_10W40E_35N75N',\
'DivQIndexLand_'+lowY1+'_'+lowY2+'_Oct_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Nov_10W40E_35N75N','DivQIndexLand_'+lowY1+'_'+lowY2+'_Dec_10W40E_35N75N',\
'DivQIndexLand_'+lowY1+'_'+lowY2+'_Jan_10W40E_35N75N']


highPlot = ['DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Mar_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Apr_10W40E_35N75N',\
'DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_May_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Jun_10W40E_35N75N',\
'DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Jul_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Aug_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Sep_10W40E_35N75N',\
'DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Oct_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Nov_10W40E_35N75N','DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Dec_10W40E_35N75N',\
'DivQIndexEOF1HighSepLand_'+highY1+'_'+highY2+'_Jan_10W40E_35N75N']

avePlot = ['DivQIndexLand_1979_2014_Mar_10W40E_35N75N','DivQIndexLand_1979_2014_Apr_10W40E_35N75N',\
'DivQIndexLand_1979_2014_May_10W40E_35N75N','DivQIndexLand_1979_2014_Jun_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jul_10W40E_35N75N','DivQIndexLand_1979_2014_Aug_10W40E_35N75N','DivQIndexLand_1979_2014_Sep_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Oct_10W40E_35N75N','DivQIndexLand_1979_2014_Nov_10W40E_35N75N','DivQIndexLand_1979_2014_Dec_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jan_10W40E_35N75N']

lowPlot = ['DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Mar_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Apr_10W40E_35N75N',\
'DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_May_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Jun_10W40E_35N75N',\
'DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Jul_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Aug_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Sep_10W40E_35N75N',\
'DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Oct_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Nov_10W40E_35N75N','DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Dec_10W40E_35N75N',\
'DivQIndexEOF1LowSepLand_'+lowY1+'_'+lowY2+'_Jan_10W40E_35N75N']

'''

months = ['Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan']

xdata = numpy.arange(len(months))

highAv = []
aveAv  = []
lowAv  = []

for aveMonth in months:
	highData3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexHigh'+month+'BarentsLand_'+highYear1+'_'+highYear2+'_'+aveMonth+'_10W40E_35N75N','r').read()
	highData2 = json.loads(highData3)
	highData = numpy.array(highData2)
	highAv.append(numpy.average(highData))

for aveMonth in months:
	average3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexLand_1979_2014_'+aveMonth+'_10W40E_35N75N','r').read()
	average2 = json.loads(average3)
	average = numpy.array(average2)
	aveAv.append(numpy.average(average))

for aveMonth in months:
	lowData3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexLow'+month+'BarentsLand_'+lowYear1+'_'+lowYear2+'_'+aveMonth+'_10W40E_35N75N','r').read()
	lowData2 = json.loads(lowData3)
	lowData = numpy.array(lowData2)
	lowAv.append(numpy.average(lowData))

highMinusAve = numpy.array(highAv) - numpy.array(aveAv)

lowMinusAve  = numpy.array(lowAv)  - numpy.array(aveAv)

fig, ax = plt.subplots(1,1,figsize = (13,6))
ax.plot(xdata, highMinusAve, '--or', linewidth = 3., label = 'high')
ax.plot(xdata, lowMinusAve, '-ob', linewidth = 3., label = 'low')
#ax.plot(xdata, numpy.array(aveAv), '-og', linewidth = 3., label = 'ave')

plt.xlabel('Month')
plt.ylabel('Latent energy (Wm-2)')
plt.xticks(xdata,month)
plt.title("Latent energy anomolies of the years with highest and lowest SIC in Aug.")
plt.grid(True)
plt.legend(loc=4)
#plt.xlim(1979,2014)

#plt.savefig('/home/nba035/plot/sicSelectedLatentTrend'+saveMonth+'.eps', format = 'eps')
plt.show()


