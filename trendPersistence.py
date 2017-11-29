import json
import numpy
import codecs
import matplotlib
import matplotlib.pyplot as plt


highPlot = ['DivQIndexLand_1982_1988_Mar_10W40E_35N75N','DivQIndexLand_1982_1988_Apr_10W40E_35N75N',\
'DivQIndexLand_1982_1988_May_10W40E_35N75N','DivQIndexLand_1982_1988_Jun_10W40E_35N75N',\
'DivQIndexLand_1982_1988_Jul_10W40E_35N75N','DivQIndexLand_1982_1988_Aug_10W40E_35N75N','DivQIndexLand_1982_1988_Sep_10W40E_35N75N',\
'DivQIndexLand_1982_1988_Oct_10W40E_35N75N','DivQIndexLand_1982_1988_Nov_10W40E_35N75N','DivQIndexLand_1982_1988_Dec_10W40E_35N75N',\
'DivQIndexLand_1982_1988_Jan_10W40E_35N75N']

avePlot = ['DivQIndexLand_1979_2014_Mar_10W40E_35N75N','DivQIndexLand_1979_2014_Apr_10W40E_35N75N',\
'DivQIndexLand_1979_2014_May_10W40E_35N75N','DivQIndexLand_1979_2014_Jun_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jul_10W40E_35N75N','DivQIndexLand_1979_2014_Aug_10W40E_35N75N','DivQIndexLand_1979_2014_Sep_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Oct_10W40E_35N75N','DivQIndexLand_1979_2014_Nov_10W40E_35N75N','DivQIndexLand_1979_2014_Dec_10W40E_35N75N',\
'DivQIndexLand_1979_2014_Jan_10W40E_35N75N']

lowPlot = ['DivQIndexLand_1984_1985_Mar_10W40E_35N75N','DivQIndexLand_1984_1985_Apr_10W40E_35N75N',\
'DivQIndexLand_1984_1985_May_10W40E_35N75N','DivQIndexLand_1984_1985_Jun_10W40E_35N75N',\
'DivQIndexLand_1984_1985_Jul_10W40E_35N75N','DivQIndexLand_1984_1985_Aug_10W40E_35N75N','DivQIndexLand_1984_1985_Sep_10W40E_35N75N',\
'DivQIndexLand_1984_1985_Oct_10W40E_35N75N','DivQIndexLand_1984_1985_Nov_10W40E_35N75N','DivQIndexLand_1984_1985_Dec_10W40E_35N75N',\
'DivQIndexLand_1984_1985_Jan_10W40E_35N75N']

xdata = numpy.arange(len(avePlot))

month = ['Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan']

highAv = []
aveAv  = []
lowAv  = []

for high in highPlot:
	highData3 = codecs.open(high,'r').read()
	highData2 = json.loads(highData3)
	highData = numpy.array(highData2)
	highAv.append(numpy.average(highData))

for avi in avePlot:
	average3 = codecs.open(avi,'r').read()
	average2 = json.loads(average3)
	average = numpy.array(average2)
	aveAv.append(numpy.average(average))

for low in lowPlot:
	lowData3 = codecs.open(low,'r').read()
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
plt.title("Latent energy anomolies of the years with highest and lowest SIC")
plt.grid(True)
plt.legend(loc=4)
#plt.xlim(1979,2014)

#plt.savefig('/home/nba035/plot/latentEnergyAnomalies.eps', format = 'eps')
plt.show()


