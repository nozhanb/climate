import json
import numpy
import codecs
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


high3 = codecs.open('DivQAverage_1982_1988_Sep_10W40E_35N75N','r').read()
high2 = json.loads(high3)
high = numpy.array(high2)

average3 = codecs.open('DivQAverage_1979_2014_Sep_10W40E_35N75N','r').read()
average2 = json.loads(average3)
average = numpy.array(average2)

low3 = codecs.open('DivQAverage_1984_1985_Sep_10W40E_35N75N','r').read()
low2 = json.loads(low3)
low = numpy.array(low2)


highIndex3 = codecs.open('DivQIndex_1982_1988_Sep_10W40E_35N75N','r').read()
highIndex2 = json.loads(highIndex3)
highIndex = numpy.array(highIndex2)

averageIndex3 = codecs.open('DivQIndex_1979_2014_Sep_10W40E_35N75N','r').read()
averageIndex2 = json.loads(averageIndex3)
averageIndex = numpy.array(averageIndex2)

lowIndex3 = codecs.open('DivQIndex_1984_1985_Sep_10W40E_35N75N','r').read()
lowIndex2 = json.loads(lowIndex3)
lowIndex = numpy.array(lowIndex2)


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
myData = highDiff.reshape((latitude123.size,longitude123.size))
#myData = lowDiff.reshape((latitude123.size,longitude123.size))
myData = numpy.flipud(myData)   #   This line is important!!!
cmap = plt.get_cmap('RdBu_r')
plotType = 'contourf'
v = numpy.arange(-200, 201, 20)
cmap = plt.cm.get_cmap('RdBu_r')
bounds = v
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
colorf = map1.contourf(x, y, myData,v,cmap=cmap,norm=norm,extend = 'both')
cb = map1.colorbar(colorf, location='bottom', pad="5%")
cb.cmap.set_under('m')
cb.cmap.set_over('g')
cb.set_label(r'Latent energy divergence (Wm-2)')

map1.drawparallels(numpy.arange(lat1,lat2 + 1,10),labels=[1,0,0,0], linewidth=0.0)
map1.drawmeridians(numpy.arange(lon1,lon2 + 1,10),labels=[0,0,0,1], linewidth=0.0)

map1.drawcoastlines()
map1.drawcountries()
ax1.set_title('high - average')
ax1.annotate('high {:.2f}, ave {:.2f},low {:.2f}, high-av {:.2f} and low-av {:.2f}'.format(numpy.average(highIndex),\
numpy.average(averageIndex), numpy.average(lowIndex), highMinusAveIndex, lowMinusAveIndex), xy = (10,40), xytext=(-10,77))
#plt.savefig('/home/nba035/plot/'+name+'.eps',dpi = 80, format = 'eps')

plt.show()

