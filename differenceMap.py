import json
import numpy
import codecs
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

month = 'Sep'
highYear1 = str(1982)
highYear2 = str(1988)
lowYear1 = str(1984)
lowYear2 = str(1985)

'''
high3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageEOF1High'+month+'Land_'+highYear1+'_'+highYear2+'_'+month+'_10W40E_35N75N','r').read()
high2 = json.loads(high3)
high = numpy.array(high2)

average3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageLand_1979_2014_'+month+'_10W40E_35N75N','r').read()
average2 = json.loads(average3)
average = numpy.array(average2)

low3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageEOF1Low'+month+'Land_'+lowYear1+'_'+lowYear2+'_'+month+'_10W40E_35N75N','r').read()
low2 = json.loads(low3)
low = numpy.array(low2)

highIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexEOF1High'+month+'Land_'+highYear1+'_'+highYear2+'_'+month+'_10W40E_35N75N','r').read()
highIndex2 = json.loads(highIndex3)
highIndex = numpy.array(highIndex2)

averageIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexLand_1979_2014_'+month+'_10W40E_35N75N','r').read()
averageIndex2 = json.loads(averageIndex3)
averageIndex = numpy.array(averageIndex2)

lowIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexEOF1Low'+month+'Land_'+lowYear1+'_'+lowYear2+'_'+month+'_10W40E_35N75N','r').read()
lowIndex2 = json.loads(lowIndex3)
lowIndex = numpy.array(lowIndex2)

contour1 = numpy.fromfile('DivQMontecarloHigh'+month+'EOF1Land_1981_1989_'+month+'_10W40E_35N75N', float, -1, ",")

'''

high3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageHigh'+month+'BarentsLand_'+highYear1+'_'+highYear2+'_'+month+'_10W40E_35N75N','r').read()
high2 = json.loads(high3)
high = numpy.array(high2)

average3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageLand_1979_2014_'+month+'_10W40E_35N75N','r').read()
average2 = json.loads(average3)
average = numpy.array(average2)

low3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQAverageLow'+month+'BarentsLand_'+lowYear1+'_'+lowYear2+'_'+month+'_10W40E_35N75N','r').read()
low2 = json.loads(low3)
low = numpy.array(low2)

highIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexHigh'+month+'BarentsLand_'+highYear1+'_'+highYear2+'_'+month+'_10W40E_35N75N','r').read()
highIndex2 = json.loads(highIndex3)
highIndex = numpy.array(highIndex2)

averageIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexLand_1979_2014_'+month+'_10W40E_35N75N','r').read()
averageIndex2 = json.loads(averageIndex3)
averageIndex = numpy.array(averageIndex2)

lowIndex3 = codecs.open('/global/work/nba035/WRF/thesis/data/test/backUp/DivQIndexLow'+month+'BarentsLand_'+lowYear1+'_'+lowYear2+'_'+month+'_10W40E_35N75N','r').read()
lowIndex2 = json.loads(lowIndex3)
lowIndex = numpy.array(lowIndex2)

contour1 = numpy.fromfile('DivQMontecarloLow'+month+'BarentsLand_1984_1985_'+month+'_10W40E_35N75N', float, -1, ",")


highDiff = high - average
lowDiff  = low  - average

highMinusAveIndex = numpy.average(highIndex) - numpy.average(averageIndex)

lowMinusAveIndex = numpy.average(lowIndex) - numpy.average(averageIndex)

lat1, lat2 = (35,75)
latitude123 = numpy.arange(lat1,lat2+0.5,0.5)
lon1, lon2 = (-10,40)
longitude123= numpy.arange(lon1, lon2+0.5,0.5)

fig1,ax1 = plt.subplots(1,1,figsize=(13,13))

map1 = Basemap(projection='cyl',llcrnrlon=lon1,llcrnrlat=lat1, \
    urcrnrlon=lon2,urcrnrlat=lat2, lon_0 = 30, lat_0 = 30,\
    rsphere=6371200., resolution = 'l', area_thresh=10000)

lons,lats = numpy.meshgrid(longitude123,latitude123)
x, y = map1(lons,lats)
#myData = highDiff.reshape((latitude123.size,longitude123.size))
myData = lowDiff.reshape((latitude123.size,longitude123.size))
finalCon = contour1.reshape((latitude123.size,longitude123.size))
finalCon = numpy.flipud(finalCon)
myData = numpy.flipud(myData)   #   This line is important since data are read from north pole to equator!!!
cmap = plt.get_cmap('RdBu_r')	#   RdBu_r
plotType = 'contourf'
cmap = plt.cm.get_cmap('RdBu_r')
cmap.set_under('c')
cmap.set_over('m')
bounds = numpy.arange(-40, 41, 5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
colorf = map1.contourf(x, y, myData,bounds,cmap=cmap,norm=norm,extend = 'both')
cset = map1.contour(x, y, finalCon,levels=[0.90,0.95],colors = ['g','y'], linewidths=[2.,2.],extent='both')
plt.clabel(cset, fmt='%1.2f', inline=1, fontsize=10)
cb = map1.colorbar(colorf, location='bottom', pad="5%")
cb.set_label(r'Latent energy divergence (Wm-2)')

map1.drawparallels(numpy.arange(lat1,lat2 + 1,10),labels=[1,0,0,0], linewidth=0.0)
map1.drawmeridians(numpy.arange(lon1,lon2 + 1,10),labels=[0,0,0,1], linewidth=0.0)

map1.drawcoastlines()
map1.drawcountries()
#ax1.set_title('High SIC anomalies in '+month )
ax1.set_title('Low SIC anomalies in '+month )

#ax1.annotate('ave: {:.2f}, high-ave: {:.2f}'.format(numpy.average(averageIndex), highMinusAveIndex), xy = (10,40), xytext=(35,75.5))

ax1.annotate('ave: {:.2f}, low-ave: {:.2f}'.format(numpy.average(averageIndex), lowMinusAveIndex), xy = (10,40), xytext=(35,75.5))

#plt.savefig('/home/nba035/plot/eof1'+month+'SicHighMinusAverageAnomalies'+month+'.eps',dpi = 80, format = 'eps')

plt.show()

