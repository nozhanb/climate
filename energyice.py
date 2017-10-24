import json
import numpy
import netCDF4
import calendar
import collections  #   to make an ordered dictionary.
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import numpy.ma as ma  #   I do not remember what it does! must be checked.
#numpy.set_printoptions(threshold=numpy.nan) #   if activated it prints out the numpy arrays in their entierty 
                                                #   regrdless of how big they are.

#ice = energyice.EnergyIce('SIC',(2010,2014),['01','02','03'],(-10,40),(35,75))
#energy = energyice.EnergyIce('DivQ',(2010,2012),['09'],(-10,40),(35,75))    for Energy files

class EnergyIce:

    def __init__(self, fileContent, year, month, longitude, latitude):
        
        #   fileContent: it has to be a either 'SIC', 'DivQ', 'DivQ'.
        #   year: must be a tuple of int. format. e.g.  (1979,2014). If you intend to do it for only one year repeat the 
        #       the year. e.g. (2010,2010)
        #   month: has to be given as a list of strings with a preceding zero. e.g. ['01','05'] which stands for Jan. 
        #       and May.
        #   longitude and latitude: must be tuples, integers and in degree (0,90), from 0-degree to 90-degree. for latitude
        #       and (0,360) for longitude. 
        #   !!!!!!!!!
        #   NOTE:(in the line totalAreaInMeter you will face problem if longitude or latitude are given in fraction 
        #       because of the multipication by 2)!!!
        #   !!!!!!!!!
        
        #!!! IMPORTANT !!!
        #   testsL Check for the last latitude at 90 and also at 0 and make sure the right value is generated.        
        #   Watch the 30 that is hard coded in averageEnergyPerCellPerYear array. It will produce wrong average for 31-day months !!!
        #   in iceReader where the average is taken you have to account for the zeros you put in there as a replace for
        #   ---> regions with ice value ---. First, check if --- realy exist in your file. If so, find a solution to find the 
        #   correct aerage value.
        #   Notice the use of ordered dictionary.
        
        self.month = month
        self.lat = ((90-latitude[0])*2,(90-latitude[1])*2)
        self.lon = ((abs(180)+longitude[0])*2,(abs(180)+longitude[1])*2)
        self.lat5 = latitude
        self.lon5 = longitude
        self.totalSize = int(((latitude[1]-latitude[0])*2 + 1)*((longitude[1]-longitude[0])*2 + 1))
        self.fileContent = fileContent
        self.filex = []
        
        self.year1, self.year2 = year[0],year[1]
        if fileContent  == 'SIC':
            for years in range(self.year1, self.year2+1):
                self.filex.append('SIC.'+str(years)+'.nc')
        elif fileContent == 'DivQ':
            for years in range(self.year1, self.year2+1):
                for months in self.month:
                    self.filex.append('DivQ.'+str(years)+'.'+months+'.nc')
        elif fileContent == 'DivD':
            for years in range(self.year1, self.year2+1):
                for months in self.month:
                    self.filex.append('DivD.'+str(years)+'.'+str(months)+'.nc')



    def sicReader(self,save = True, outputName = None):
        ''' This function reads in a .nc file and returns a dictionary of differnt indices for indicated months.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        totalIceFraction = []
        for fileName in self.filex:
            for monthNo in self.month:
                printName = fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'
                loaded_file = netCDF4.Dataset(fileName)
                latiInverse = loaded_file.variables['g0_lat_1'][int(self.lat[1]):int(self.lat[0]) + 1]
                lati = numpy.flipud(latiInverse.reshape((latiInverse.size,1))).flatten('C')
                longi = loaded_file.variables['g0_lon_2'][int(self.lon[0]):int(self.lon[1]) + 1]
                numerator = numpy.array([])
                denominator = numpy.array([])
                latCounter = numpy.arange(self.lat[0] + 1, self.lat[1], -1)
                lonCounter = numpy.arange(self.lon[0], self.lon[1] + 1, 1)

                counter1 = 0                    
                for i in lati:      #   lati is a list of all latitude between an initial and final latitudes given by the user.
                    lat1 = i
                    counter2 = 0
                    for j in longi: #   similar to lati but it stands for longitudes
                        surat = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][monthNo,latCounter[counter1],lonCounter[counter2]]*(numpy.cos(lat1*(3.14/180))*0.5)*0.5)
                        makh = numpy.array((numpy.cos(lat1*(3.14/180))*0.5)*0.5)    #   area of each cell (including the conversion between radian and degree).
                        numerator = numpy.append(numerator,numpy.array(surat))      #   stores the value of the ice fraction at each cell at each step.
                        denominator = numpy.append(denominator,makh)                #   stores the value of the area of each cell at each step.
                        if int(lat1) == 90:     # if the upper edge (upper latitude) is 90, equate the are of that region to zero because the cos(90) = 0.
                            numerator = numpy.append(numerator,0)       #   store the values similar to the line above.
                            denominator = numpy.append(denominator,0)   #   similar to above.
                        counter2 += 1   #   at each step increment by 1.
                    counter1 += 1       #   increment by 1.
                
                sumnum = numpy.sum(numerator)       #   This is where we run a sum function over the series of ice fraction at each cell to find the total ice fraction.
                sumdenom = numpy.sum(denominator)   #   We run a sum over a series of area of each cell to find the total area of the region.
                print 'Dic name --->', printName
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = float(sumnum)/float(sumdenom)
                totalIceFraction.append(float(sumnum)/float(sumdenom))
        
                if save == True:
                    if self.lat5[0] < 0.0:
                        latFir = 'S'
                    elif self.lat5[0] >= 0.0:
                        latFir = 'N'
                    if self.lat5[1] < 0.0:
                        latSec = 'S'
                    elif self.lat5[1] >= 0.0:
                        latSec = 'N'
                    if self.lon5[0] < 0.0:
                        lonFir = 'W'
                    elif self.lon5[0] >= 0.0:
                        lonFir = 'E'
                    if self.lon5[1] < 0.0:
                        lonSec = 'W'
                    elif self.lon5[1] >= 0.0:
                        lonSec = 'E'
                    if outputName == None:
                        nameToWrite = 'iceFraction'
                    elif outputName is not None:
                        nameToWrite = outputName
                    ff = open(nameToWrite+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[monthNo][0:3]+'_'\
                    +str(self.lon5[0])+lonFir+str(self.lon5[1])+lonFir+'_'+str(self.lat5[0])+latFir+str(self.lat5[1])+\
                    latFir,'w')
                    ff.write(json.dumps(totalIceFraction))
                    ff.close()

        return self.output_dic


        
    def energyReader(self,save = True, outputName = None):
        self.averagePerYearDictionary = collections.OrderedDict()   #   This dictionary collects the average energy during each month for each montha for each year for each cell.
        self.energyIndexDictionary = collections.OrderedDict()
        totalAverage = numpy.zeros((self.totalSize,))
        newArray = numpy.zeros((self.totalSize,))
        totalIndexList = []

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(fileName)
            latiInverse = loaded_file.variables['lat_NH'][int(self.lat[1]):int(self.lat[0]) + 1]
            lati = numpy.flipud(latiInverse.reshape((latiInverse.size,1))).flatten('C')
            longi = loaded_file.variables['lon'][int(self.lon[0]):int(self.lon[1]) + 1]
            averageEnergyPerDayOverEurope = []
            cosValueArray = numpy.array([])
            cellAreaEnergyArray = numpy.array([])
            energyPerCellArrayList = []
            latCounter = numpy.arange(self.lat[0] + 1, self.lat[1], -1)
            lonCounter = numpy.arange(self.lon[0], self.lon[1] + 1, 1)

            for time in range(0,30):
                energyPerCellArray = numpy.array([])
                counter1 = 0
                for i in lati:      #   lati is a list of all latitude between an initial and final latitudes given by the user.
                    lat2 = i + 0.5  
                    counter2 = 0
                    for j in longi: #   similar to lati but it stands for longitudes
                        cellEnergyWattPerMeter = numpy.array(loaded_file.variables['Divergence'][time,latCounter[counter1],lonCounter[counter2]])
                        if int(lat2) == 90:
                            cellEnergyWattPerMeter = 0      # I am not sure about this line. Needs to be tested !!!
                        cellAreaEnergy = cellEnergyWattPerMeter*numpy.cos(lat2*(3.14/180))
                        cellAreaEnergyArray = numpy.append(cellAreaEnergyArray,cellAreaEnergy)
                        energyPerCellArray = numpy.append(energyPerCellArray,cellEnergyWattPerMeter)
                        counter2 += 1
                    cosValue = numpy.cos(lat2*(3.14/180.))
                    cosValueArray = numpy.append(cosValueArray,cosValue)
                    counter1 += 1
                energyPerCellArrayList.append(energyPerCellArray)
                averageEnergyPerDayOverEurope.append(float(numpy.sum(cellAreaEnergyArray))/float((self.lon5[1]-self.lon5[0])\
                *2*numpy.sum(cosValueArray)))
            counter5 = 0
            for counter3 in energyPerCellArrayList:
                newArray = numpy.add(newArray,counter3)
                newArray2 = newArray
                newArray = newArray2
                counter5 += 1
            finalIndex = sum(averageEnergyPerDayOverEurope)/30.
            self.energyIndexDictionary['energyIndexPerYear'+fileName[5:12]] = finalIndex
            totalIndexList.append(finalIndex)
            self.averagePerYearDictionary['averageEnergyPerCellPerYear'+fileName[5:12]] = newArray2/30.
            totalAverage = numpy.add(totalAverage,self.averagePerYearDictionary['averageEnergyPerCellPerYear'+fileName[5:12]])
            print fileName[5:12], '-------------> DONE!!!' 

        if save == True:
            if self.lat5[0] < 0.0:
                latFir = 'S'
            elif self.lat5[0] >= 0.0:
                latFir = 'N'
            if self.lat5[1] < 0.0:
                latSec = 'S'
            elif self.lat5[1] >= 0.0:
                latSec = 'N'
            if self.lon5[0] < 0.0:
                lonFir = 'W'
            elif self.lon5[0] >= 0.0:
                lonFir = 'E'
            if self.lon5[1] < 0.0:
                lonSec = 'W'
            elif self.lon5[1] >= 0.0:
                lonSec = 'E'
            if outputName == None:
                nameToWrite1 = self.fileContent+'Index'
                nameToWrite2 = self.fileContent+'Average'
                
            if len(self.filex) == 1:
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            ff1.write(json.dumps(totalIndexList))
            totalAverageList = totalAverage.tolist()
            ff2.write(json.dumps(totalAverageList))
            ff1.close()
            ff2.close()
        return totalAverage, numpy.array(totalIndexList)
        
        
        
def europeMap(data,numberOfYears,longitude,latitude, plotType = 'mesh',figTitle = None):
    lat1, lat2 = latitude
    latitude123 = numpy.arange(lat1,lat2+0.5,0.5)
    lon1, lon2 = longitude
    longitude123= numpy.arange(lon1, lon2+0.5,0.5)
    years = int(numberOfYears)

    fig1,ax1 = plt.subplots(1,1,figsize=(13,13))

    ax1.set_title(str(figTitle))

    map1 = Basemap(projection='cyl',llcrnrlon=-10.0,llcrnrlat=35.0, \
    urcrnrlon=40.0,urcrnrlat=75.0, lon_0 = 30, lat_0 = 30,\
    rsphere=6371200., resolution = 'l', area_thresh=10000)

    lons,lats = numpy.meshgrid(longitude123,latitude123)
    x, y = map1(lons,lats)
    data = data/years
    myData = data.reshape((81,101))
    cmap = plt.get_cmap('jet')
    if plotType == 'mesh':
        colormesh = map1.pcolormesh(x,y,myData,cmap = cmap)
        cb = map1.colorbar(colormesh, location='bottom', pad="5%")
    elif plotType == 'contourf':
        colorf = map1.contourf(x, y, myData, cmap=cmap)
        cb = map1.colorbar(colorf, location='bottom', pad="5%")   
#    map1.colorbar(colormesh)
    cb.set_label(r'Latent Energy Flux (W/m^2)')

    map1.drawparallels(numpy.arange(35,76,10),labels=[1,0,0,0], linewidth=0.0)
    map1.drawmeridians(numpy.arange(-10,41,10),labels=[0,0,0,1], linewidth=0.0)

    map1.drawcoastlines()
    map1.bluemarble()
    map1.etopo()
    map1.drawcountries()

#    cbar = map1.colorbar(ice,location='right',pad="3%")

    plt.show()    

        
def plotIndex(xdata, ydata, numberofmonth, overplot = False, inputClass = 'array'):

    #   xdata: give the initial and final year in int and tuple form.
		    #   ydata: is the output from the ice_reader function. It is a dictionary.
    #   numberofmonth: is the number of month that is covered in the ydata. It is of the int format.
    #   overplot: if it is True it will the function will generate a plot with all the given month overplotted.
    #   if it is False the function will generate a seperate graph for each month (no overplot).
    #   inputClass: can be array for a simple numpy array or dictionary. It must be given in string format.

    #   Tests: check for the order of data both in year and in ydata to see if the right year is associated 
    #   with the right year.
    #   Also fix the labels and tickmarks (there are many tick marks).


    colour = ['g','r','k','m','b']
    x_axis = xdata
    y_axis = ydata
    valueArray = numpy.array([])
    for key, value in y_axis.iteritems():
	valueArray = numpy.append(valueArray, value)

    arrayLength = valueArray.size
    interval = int(arrayLength/int(numberofmonth))
    intervalArray = numpy.arange(0,arrayLength+1,interval)  #   +1 so the end point is included.
    
    if overplot == True:
        fig1, ax1 = plt.subplots(1, figsize=(15,6))
        ini = intervalArray[0]
        month = 1
        for counter1 in range(numberofmonth):
            end = intervalArray[counter1+1]
            yplot = valueArray[ini:end]
            yearArray = numpy.arange(x_axis[0], x_axis[1]+1)
            ax1.plot(yearArray,yplot,color = colour[counter1], linestyle ='-', linewidth = 3.0, label = calendar.month_name[month])
            ini = end
            ax1.set_xlim(x_axis[0], x_axis[1])
            ax1.grid('on')
            month += 1
        plt.legend(bbox_to_anchor=(0.15, 0.25))
    elif overplot == False:
        if int(numberofmonth%2) == 0:
            fig1, ax1 = plt.subplots(int(numberofmonth/2),2, figsize=(15,6))
        elif int(numberofmonth%2) != 0:
            fig1, ax1 = plt.subplots(int(numberofmonth),1, figsize=(15,6), sharex = True)
        
        ini = intervalArray[0]
        row = 0
        for counter2 in range(numberofmonth):
            end = intervalArray[counter2+1]
            yplot = valueArray[ini:end]
            yearArray = numpy.arange(x_axis[0], x_axis[1]+1)
            ax1[row].plot(yearArray, yplot, color = colour[counter2], linestyle ='-', linewidth = 3.0)
            ax1[row].set_xlim(x_axis[0], x_axis[1])
            ax1[row].grid('on')
            ini = end
            row +=1
            
            #   Tick-adjusment section:
            if row == 1:
                ytick = ax1[1].yaxis.get_major_locator
#                   if len(ytick)%2 == 0:
                       
            
        fig1.subplots_adjust(hspace=0)

    plt.xlabel('Year', fontsize = 17)
    plt.ylabel('SSI Index', fontsize = 17)
    plt.title('Sea Surface Ice Index for Jan. Feb. and Mar. 1979-2016 for lat (75,90) and lon (0,360)')
    
    plt.show()















