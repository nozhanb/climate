import json
import numpy
import codecs
import netCDF4
import calendar
import collections  #   to make an ordered dictionary.
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import numpy.ma as ma  #   It is used to work with MASKED arrays. See numpy manual for more details on masked arrays.
#numpy.set_printoptions(threshold=numpy.nan) #   if activated it prints out the numpy arrays in their entierty 
                                                #   regrdless of how big they are.

#ice = energyice.EnergyIce('SIC',(2010,2014),['01','02','03'],(-10,40),(35,75))
#energy = energyice.EnergyIce('DivQ',(2010,2012),['09'],(-10,40),(35,75))    for Energy files

class EnergyIce:

    def __init__(self, fileContent, year, month, longitude, latitude, allYears = True):
        
        #   fileContent: it has to be a either 'SIC', 'DivQ', 'DivQ'.
        #   year: must be a tuple of int. format. e.g.  (1979,2014). If you intend to do it for only one year repeat the 
        #       the year. e.g. (2010,2010)
        #   month: has to be given as a list of strings with a preceding zero. e.g. ['01','05'] which stands for Jan. 
        #       and May.
        #   longitudes and latitudes must be given in tuple format and be aither integer or with 0.5 decimal. Also, the 
        #   sign convention has to be followed. Meaning that for western and southern hemisphere the sign must be negative
        #   and positive for northern and eastern hemisphere. e.g. longitude = (-180.,47.5), latitude=(-43.5,69.5).
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
        totalLon = numpy.arange(-180,180,0.5)
        totalLat = numpy.arange(90,29.5,-0.5)
        lonInitialIndex = numpy.where(totalLon==float(self.lon5[0]))
        lonEndIndex = numpy.where(totalLon==float(self.lon5[1]))
        latInitialIndex = numpy.where(totalLat==float(self.lat5[0]))
        latEndIndex = numpy.where(totalLat==float(self.lat5[1]))
        self.finalLon = numpy.arange(lonInitialIndex[0][0], lonEndIndex[0][0] + 1)
        self.finalLat = numpy.arange(latEndIndex[0][0], latInitialIndex[0][0] + 1)  #   Notice that latitude is ordered differently.
        
        if allYears == True:
            self.yearRange = range(self.year1, self.year2+1)        
        elif allYears == False:
            self.yearRange = year
        
        if fileContent  == 'SIC':
            for years in self.yearRange:
                    self.filex.append('SIC.'+str(years)+'.nc')
        elif fileContent == 'DivQ':
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivQ.'+str(years)+'.'+months+'.nc')
        elif fileContent == 'DivD':
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivD.'+str(years)+'.'+str(months)+'.nc')
            self.loaded_file = netCDF4.Dataset('DivD.'+str(years)+'.'+str(months)+'.nc')
            self.fastLat = numpy.array(self.loaded_file.variables['lat_NH'][:])
            self.fastLon = numpy.array(self.loaded_file.variables['lon'][:])



    def sicReader(self,save = True, outputName = None):
        ''' This function reads in a .nc file and returns a dictionary of differnt indices for indicated months.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        totalIceFraction = []
        for fileName in self.filex:
            for monthNo in self.month:
                printName = fileName+'_'+str(monthNo)+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'
                loaded_file = netCDF4.Dataset(fileName)     #   Notice the following arrays are MASKED arrays!!!
                fastLat = numpy.ma.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
                fastLon = numpy.ma.array(loaded_file.variables['g0_lon_2'][self.finalLon[0]:self.finalLon[-1]+1])
                cosValue = numpy.cos(fastLat*(3.14/180.))
                cosValue = cosValue.reshape(self.finalLat.size, 1)
                maskedIceFraction = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][int(monthNo),self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])*cosValue
                #   There are some cells in which the value of ice fraction is e+29 (i.e. regions covered with land) so we set them to zero.
#                landCoor = numpy.ma.masked_less(maskedIceFraction,1.5).filled(0.)
#                print 'This is the dimension of the landcoor: ', landCoor.shape               
                iceFractionMaskedWithZero = numpy.ma.masked_greater(maskedIceFraction,10).filled(0.)
                numerator = numpy.sum(numpy.sum(iceFractionMaskedWithZero,0))
                denominator = fastLon.size*numpy.sum(numpy.sum(cosValue,0))
                iceIndex = numerator/denominator
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = iceIndex
                totalIceFraction.append(iceIndex)
                print printName, '------> DONE!!!'
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
                    ff = open(nameToWrite+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(monthNo)][0:3]+'_'\
                    +str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+\
                    latSec,'w')
                    ff.write(json.dumps(totalIceFraction))
                    ff.close()
        return numpy.array(totalIceFraction), landCoor



##########################
#   REMEMBER TO CHECK FOR LEAP YEARS!!!
##########################
        
    def energyReader(self,save = True, outputName = None):
        self.averagePerYearDictionary = collections.OrderedDict()   #   This dictionary collects the average energy during each month for each montha for each year for each cell.
        self.energyIndexDictionary = collections.OrderedDict()
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        totalFinalIndex = numpy.array([])

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(fileName)
            days = calendar.monthrange(int(fileName[5:9]),int(fileName[10:12]))[1]
            fastLat = numpy.array(loaded_file.variables['lat_NH'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(1,self.finalLat.size,1)
            fastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            EnergyAverage = EnergyAverage + numpy.sum(fastEnergy,0)/days
            numerator = numpy.sum(numpy.sum(fastEnergy*cosValue,0)/days)
            denominator = self.finalLon.size*numpy.sum(cosValue)
            finalIndex = numerator/denominator
            totalFinalIndex = numpy.append(totalFinalIndex,finalIndex)
            print fileName[5:12], '-----> DONE!!!'
#            print 'Average index (Wm-2)----->', finalIndex
        print 'Average index (Wm-2)----->', numpy.average(totalFinalIndex)
        if save == True:
            if self.lat5[0] < 0.0:
                latFir = 'S'
            elif self.lat5[0] >= 0.0:
                latFir = 'N'
            if self.lat5[1] < 0.0:
                latSec = 'S'
            elif self.lat5[1] >= 0.0:
                latSec = 'N'
            if self.lon5[0] <   0.0:
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

            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((EnergyAverage/float(len(self.yearRange))).tolist()))
            ff1.close()
            ff2.close()
        return EnergyAverage/float(len(self.yearRange)), totalFinalIndex



    def landEnergyReader(self,save = True, outputName = None):
        self.averagePerYearDictionary = collections.OrderedDict()   #   This dictionary collects the average energy during each month for each montha for each year for each cell.
        self.energyIndexDictionary = collections.OrderedDict()
        totalFinalIndex = numpy.array([])
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        europeCoor = numpy.fromfile('europeCoordinate',float,-1,",")
        europeCoor = europeCoor.reshape(30,81,101)   #   The number of days in this line does not count for LEAP years.
        maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(europeCoor, 1))
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(fileName)
            days = calendar.monthrange(int(fileName[5:9]),int(fileName[10:12]))[1]
            fastLat = numpy.array(loaded_file.variables['lat_NH'][self.finalLat[0]:self.finalLat[-1]+1])
            fastLon = numpy.array(loaded_file.variables['lon'][self.finalLon[0]:self.finalLon[-1]+1])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(1,self.finalLat.size,1)
            fastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastEnergy = numpy.ma.masked_where(maskOfEurope,fastEnergy).filled(0.)
            EnergyAverage = EnergyAverage + numpy.sum(fastEnergy,0)/days
            localIndex = numpy.array([])
            for i in range(self.finalLat.size):
                numerator = numpy.sum(numpy.sum(fastEnergy[:,i:i+1,:]*cosValue[:,i:i+1,:],0)/days)
                denominator = numpy.count_nonzero(fastEnergy[:,i:i+1,:])*numpy.sum(cosValue[:,i:i+1,:])
                if denominator == 0.:
                    localIndex = numpy.append(localIndex,0.)
                elif denominator != 0.:
                    localIndex = numpy.append(localIndex,numerator/denominator)
#                print i,'----->',fastLat[i],fastEnergy[0:1,i:i+1,0:1], denominator, localIndex
#                for index in numpy.ndindex(fastEnergy[0:1,i:i+1,0:1].shape):
#                    print index,'--->',i,fastLat[i],fastLon[index[2]] ,numpy.count_nonzero(fastEnergy[:,i:i+1,:])#'----->',denominator,' = ',fastEnergy[index[0],i,index[2]], '*', cosValue[:,i,:]

#                print fastEnergy[:,i:i+1,:],'-----------' ,cosValue[:,i:i+1,:],fastEnergy[:,i:i+1,:].count(),numerator/denominator
#            numerator = numpy.sum(numpy.sum(fastEnergy*cosValue,0)/days)
#            denominator = 124560.*numpy.sum(cosValue)   #   124560  the  number of unmasked cells
#            finalIndex = numerator/denominator
#            finalIndex.append(numpy.sum(localIndex))
#            totalFinalIndex = numpy.append(totalFinalIndex,finalIndex)
            totalFinalIndex = numpy.append(totalFinalIndex,numpy.sum(localIndex))
            print fileName[5:12], '-----> DONE!!!'
#            print 'Average index (Wm-2)----->', finalIndex
        print totalFinalIndex
        print 'Average index (Wm-2)----->', numpy.average(totalFinalIndex)
        
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
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec+'land','w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec+'land','w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(self.month[0])][0:3]+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((EnergyAverage/float(len(self.yearRange))).tolist()))
            ff1.close()
            ff2.close()
        return EnergyAverage/float(len(self.yearRange)), totalFinalIndex




    def fitting(self, xData, yData, fitFunction = None, **params):
        import scipy
        if fitFunction is None:
            pars = numpy.polyfit(xData,yData,2)
            poly = numpy.poly1d(pars)
            return poly(xData), pars
        elif fitFunction is not None:
            def defined(xData,params):
                return fitFunction
            popt, pcov = scipy.optimize.curve_fit(defined, xData, yData)
            return defined(xData,*popt)
        
        
        

    def smoothing(self, inputData,longitude,latitude,smoothingLength = 1):
        "This function applies a smoothing technique to a given region and makes the features of that region more clear."
        #   A square is drwan around each point. To do so, we need to introduce 9 data points (including the central).
        #   This function smooths out longitudinally. It starts from the starting point and smooths out along the starting 
        #   latitude. Once it reaches the last longitude (indicated by longitude in endingPoint) it moves to the upper 
        #   latitude and the process is repeated until it reached the very last point (endingPoint).
        #   To make things easier it is best to think of the startingPoint and endingPoint as the lower-left corner and 
        #   the upper-right corner of the smoothing frame/region.
        #   smoothingLength: 1 means the immediate neighbours (i.e. half a degree distance between the central point and 
        #   the neighbouring points). If the user is intrested in one degree distance then number 2 must be placed. 
        #   Generally, the following formula must be followed to set the desired distance/radius around the central point
        #   desired distance (in degree) = smoothingLength*0.5
        
        '''
        self.sLength = smoothingLength
        self.startLat = self.totalLat[0]
        self.startLon = self.totalLon[0]
        self.endLat = self.totalLat[-1]
        self.endLon = self.totalLon[-1]

        # The coordinate of the neighbouring points:
        self.p0 = (self.startLat-self.sLength,self.startLon-self.sLength)
        self.p1 = (self.startLat-self.sLength,self.startLon)
        self.p2 = (self.startLat-self.sLength,self.startLon+self.sLength)
        self.p3 = (self.startLat,self.startLon-self.sLength)
        self.pCenter = (self.startLat,self.startLon)     #   This is the central/original point!!!
        self.p5 = (self.startLat,self.startLon+self.sLength)
        self.p6 = (self.startLat+self.sLength,self.startLon-self.sLength)
        self.p7 = (self.startLat+self.sLength,self.startLon)
        self.p8 = (self.startLat+self.sLength,self.startLon+self.sLength)

        '''
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(fileName)
            days = calendar.monthrange(int(fileName[5:9]),int(fileName[10:12]))[1]
            fastEnergy = numpy.array(loaded_file.variables['Divergence'][:,:,:])
            fastLat = numpy.array(loaded_file.variables['lon'][:])
            fastLon = numpy.array(loaded_file.variables['lat_NH'][:])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(1,fastLat.size,1)
#        newLon = numpy.arange(0,720)
#        newLon = numpy.append(newLon,0)
#        newLon = numpy.insert(newLon,0,719,0)
        avEnergy = numpy.sum(fastEnergy,0)/days
        cosValue = numpy.cos(fastLat*3.14/180.).reshape(fastLat.size,1)
        cosEnergy = avEnergy*cosValue
        
        for lati in self.finalLat:
                upArea = cosEnergy[lati,:]
                downArea = cosEnergy[lati + 1,:]
                twoArea = (upArea + downArea)/2.
                for loni in self.finalLon:
                    finalArea = numpy.sum(twoArea[loni - 1:loni + 1])/2.
                
        '''     
           0.5     0.5
        6-------7-------8
        |       |       |
        |       |       |
        3-------C-------5
        |       |       |
        |       |       |
        0-------1-------2

        '''

        
def europeMap(data,longitude,latitude,sphereProjection= False, figTitle = None,\
     plotType = 'contourf', store = False, name = None):

    lat1, lat2 = latitude
    latitude123 = numpy.arange(lat1,lat2+0.5,0.5)
    lon1, lon2 = longitude
    longitude123= numpy.arange(lon1, lon2+0.5,0.5)

    fig1,ax1 = plt.subplots(1,1,figsize=(13,13))

    map1 = Basemap(projection='cyl',llcrnrlon=lon1,llcrnrlat=lat1, \
    urcrnrlon=lon2,urcrnrlat=lat2, lon_0 = 30, lat_0 = 30,\
    rsphere=6371200., resolution = 'l', area_thresh=10000)

    if sphereProjection == True:
        map1 = Basemap(projection='npstere', boundinglat=60, llcrnrlon=lon1,llcrnrlat=lat1, \
        urcrnrlon=lon2, urcrnrlat=lat2, lon_0 = 270, round=True)
        
    lons,lats = numpy.meshgrid(longitude123,latitude123)
    x, y = map1(lons,lats)
    myData = data.reshape((latitude123.size,longitude123.size))
    myData = numpy.flipud(myData)   #   This line is important!!!
    cmap = plt.get_cmap('RdBu_r')
    if plotType == 'mesh':
        colormesh = map1.pcolormesh(x,y,myData,cmap = cmap)
        cb = map1.colorbar(colormesh, location='bottom', pad="5%")
    elif plotType == 'contourf':
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

    if sphereProjection == True:
        map1.drawmeridians([-150,-120,-90,-60,-30,0,30,60,90,120,150,180], labels=[1,1,1,1], linewidth=1.0)
        map1.drawparallels(numpy.arange(60,81,10),labels=[1,1,1,1], linewidth=1.0)

    map1.drawcoastlines()
#    map1.bluemarble()
#    map1.etopo()
    map1.drawcountries()
    ax1.set_title(str(figTitle))
    ax1.annotate('Ave = 5.01', xy = (10,40), xytext=(40,76))
    if store == True:
        plt.savefig('/home/nba035/plot/'+name+'.eps',dpi = 80, format = 'eps')

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















