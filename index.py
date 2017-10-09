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

#data = index.index((2010,2014),['01','02','03'],(40,70),(10,40),fileContent = 'ice',all_years = True)

class EnergyIce:

    def __init__(self, year, month, latitude, longitude, fileContent = 'energy', file_name = None, all_years = False):
        
        #  file_name must be a list of file names in string format. e.g. ['file1','file2','file3']
        #  month has to be given as a list of integers e.g. [1,5] which stands for Jan. and May.
        #  latitude and longitude must be tuples, float and in degree (0,90), from 0-degree to 90-degree.
        #  all_years: if True it will find the ice cover over all the years between two given years in the file_name.
        
        #   testsL Check for the last latitude at 90 and also at 0 and make sure the right value is generated.        
        #   Watch the 30 that is hard coded in averageEnergyPerCellPerYear array. It will produce wrong average for 31-day months !!!
        #   in iceReader where the average is taken you have to account for the zeros you put in there as a replace for
        #   ---> regions with ice value ---. First, check if --- realy exist in your file. If so, find a solution to find the 
        #   correct aerage value.
                
        
        self.month = month  #   must be a list of number in string format e.g. ['01','09','12'] ---> Jan. Sep. and Dec.
        self.lat = ((90-latitude[0])*2,(90-latitude[1])*2)
        self.lon = ((abs(180)+longitude[0])*2,(abs(180)+longitude[1])*2)
        self.lat5 = latitude
        self.lon5 = longitude
        self.totalSize = int(((latitude[1]-latitude[0])*2 + 1)*((longitude[1]-longitude[0])*2 + 1))
        self.year1, self.year2 = year
        self.fileContent = fileContent
        if all_years == False:
            self.filex = file_name
        elif all_years == True:
            if fileContent  == 'ice':
                self.filex = []
                for years in range(self.year1, self.year2+1):
                    self.filex.append('SIC.'+str(years)+'.nc')
            elif fileContent == 'energy':
                self.filex = []
                for years in range(self.year1, self.year2+1):
                    for months in self.month:
                        self.filex.append('DivD.'+str(years)+'.'+months+'.nc')
                
                
    def iceReader(self):
        ''' This function reads in a .nc file and returns a dictionary of differnt indices for indicated months.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        for fileName in self.filex:
            for monthNo in self.month:
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
                    lat2 = i + 0.5  
                    counter2 = 0
                    for j in longi: #   similar to lati but it stands for longitudes
                        surat = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][monthNo,latCounter[counter1],lonCounter[counter2]]*(numpy.cos(lat2*(3.14/180))*0.5)*0.5)
                        makh = numpy.array((numpy.cos(lat2*(3.14/180))*0.5)*0.5)    #   area of each cell (including the conversion between radian and degree).
                        numerator = numpy.append(numerator,numpy.array(surat))      #   stores the value of the ice fraction at each cell at each step.
                        denominator = numpy.append(denominator,makh)                #   stores the value of the area of each cell at each step.
                        if int(lat2) == 90:     # if the upper edge (upper latitude) is 90, equate the are of that region to zero because the cos(90) = 0.
                            numerator = numpy.append(numerator,0)       #   store the values similar to the line above.
                            denominator = numpy.append(denominator,0)   #   similar to above.
                        counter2 += 1   #   at each step increment by 1.
                    counter1 += 1       #   increment by 1.
                
                sumnum = numpy.sum(numerator)       #   This is where we run a sum function over the series of ice fraction at each cell to find the total ice fraction.
                sumdenom = numpy.sum(denominator)   #   We run a sum over a series of area of each cell to find the total area of the region.
                print 'Dic name --->', fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = float(sumnum/sumdenom)
        return self.output_dic
        
        
    def energyReader(self,index = False,save = True):
        self.averagePerYearDictionary = collections.OrderedDict()   #   This dictionary collects the average energy during each month for each montha for each year for each cell.
	self.energyIndexDictionary = collections.OrderedDict()
        totalAverage = numpy.zeros((self.totalSize,))
	totalIndexList = []
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(fileName)
            latiInverse = loaded_file.variables['lat_NH'][int(self.lat[1]):int(self.lat[0]) + 1]
            lati = numpy.flipud(latiInverse.reshape((latiInverse.size,1))).flatten('C')
            longi = loaded_file.variables['lon'][int(self.lon[0]):int(self.lon[1]) + 1]
            averageEnergyPerCellPerYear = numpy.array([])
            latCounter = numpy.arange(self.lat[0] + 1, self.lat[1], -1)
            lonCounter = numpy.arange(self.lon[0], self.lon[1] + 1, 1)

            counter1 = 0                    
            for i in lati:      #   lati is a list of all latitude between an initial and final latitudes given by the user.
                lat2 = i + 0.5  
                counter2 = 0
                for j in longi: #   similar to lati but it stands for longitudes
                    totalEnergy = []
                    for time in range(0,30):    #   this way we store the total energy for each cell over 30 days in the energyPerCEll array
                            #   This is where we find the average energy per cell per year.
                        cellEnergy = numpy.array(loaded_file.variables['Divergence'][time,latCounter[counter1],lonCounter[counter2]]*(numpy.cos(lat2*(3.14/180))*0.5)*0.5)
			print 'energy per cell --->', cellEnergy
                        if int(lat2) == 90:
                            cellEnergy = 0      #   I am not sure about this line. Needs to be tested !!!
                        totalEnergy.append(cellEnergy)
                    averageEnergyPerCellPerYear = numpy.append(averageEnergyPerCellPerYear, sum(totalEnergy)/30)  #   We divide by 30 days to find the average for September.
                    counter2 += 1   #   at each step increment by 1.
                counter1 += 1       #   increment by 1.
            if index == True:
		print 'Sum over all cells and totalSize: --------->', numpy.sum(averageEnergyPerCellPerYear), self.totalSize
		finalIndex = numpy.sum(averageEnergyPerCellPerYear)/self.totalSize
		self.energyIndexDictionary['energyIndexPerYear'+fileName[5:12]] = finalIndex
		print 'energyIndexPerYear'+fileName[5:12], '----->', 'Done!!!'
		totalIndexList.append(finalIndex)
	    elif index == False:
		self.averagePerYearDictionary['averageEnergyPerCellPerYear'+fileName[5:12]] = averageEnergyPerCellPerYear
        	print 'averageEnergyPerCellPerYear'+fileName[5:12], '----->', 'Done!!!'
	        totalAverage = numpy.add(totalAverage,self.averagePerYearDictionary['averageEnergyPerCellPerYear'+fileName[5:12]])
        if save == True:
	    if index == True:
		ff = open('energyIndex','w')
		ff.write(json.dumps(totalIndexList))
	    elif index == False:
		totalAverageList = totalAverage.tolist()
		ff = open('energyAverage','w')
		ff.write(json.dumps(totalAverageList))
            ff.close()
	if index == True:
		return numpy.array(totalIndexList)
        elif index == False:
		return totalAverage
        
        
def europeMap(data,numberOfYears,longitude,latitude):
    lat1, lat2 = latitude
    latitude123 = numpy.arange(lat1,lat2+0.5,0.5)
    lon1, lon2 = longitude
    longitude123= numpy.arange(lon1, lon2+0.5, 0.5)
    years = int(numberOfYears)

    fig1,ax1 = plt.subplots(1,1,figsize=(13,13))

    ax1.set_title(r'September mean of Div(vD) over Europe (1979-2014)')

    map1 = Basemap(projection='cyl',llcrnrlon=-10.0,llcrnrlat=35.0, \
    urcrnrlon=40.0,urcrnrlat=75.0, lon_0 = 30, lat_0 = 30,\
    rsphere=6371200., resolution = 'l', area_thresh=10000)

    lons,lats = numpy.meshgrid(longitude123,latitude123)
    x, y = map1(lons,lats)
    data = data/years
    myData = data.reshape((81,101))
    cmap = plt.get_cmap('jet')
    colormesh = map1.pcolormesh(x,y,myData,cmap = cmap)
#    ice = map1.contourf(x, y, myData, cmap=plt.cm.jet)
   
#    map1.colorbar(colormesh)
    cb = map1.colorbar(colormesh, location='bottom', pad="3%")
    cb.set_label(r'Dry Static Energy (Watt)')

    map1.drawcoastlines()
    map1.bluemarble()
    map1.etopo()
    map1.drawcountries()

#    cbar = map1.colorbar(ice,location='right',pad="3%")

    plt.show()    

        
def plotIndex(xdata, ydata, numberofmonth, overplot = False):

    #   xdata: give the initial and final year in int and tuple form.
    #   ydata: is the output from the ice_reader function. It is a dictionary.
    #   numberofmonth: is the number of month that is covered in the ydata. It is of the int format.
    #   overplot: if it is True it will the function will generate a plot with all the given month overplotted.
    #   if it is False the function will generate a seperate graph for each month (no overplot).

    #   Tesets: check for the order of data both in year and in ydata to see if the right year is associated 
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















