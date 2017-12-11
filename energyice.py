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
numpy.set_printoptions(threshold=numpy.nan) #   if activated it prints out the numpy arrays in their entierty 
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
        self.year = year
        self.days = calendar.monthrange(int(self.year[0]),int(self.month[0]))[1]
        self.monthName = calendar.month_name[int(self.month[0])][0:3]    # this line creates the output month to be used in output file's name.

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
            self.path = './dmDivQ_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivQ.'+str(years)+'.'+months+'.nc')
        elif fileContent == 'DivD':
            self.path = './dmDivD_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivD.'+str(years)+'.'+str(months)+'.nc')
#            self.loaded_file = netCDF4.Dataset('DivD.'+str(years)+'.'+str(months)+'.nc')
#            self.fastLat = numpy.array(self.loaded_file.variables['lat_NH'][:])
#            self.fastLon = numpy.array(self.loaded_file.variables['lon'][:])



    def sicReader(self,save = True, outputName = None, landMask = False):
        ''' This function reads in a .nc file and returns a dictionary of differnt indices for indicated months.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        totalIceFraction = []
        path = './sic/'
        for fileName in self.filex:
            for monthNo in self.month:
                printName = fileName+'_'+str(monthNo)
                loaded_file = netCDF4.Dataset(path+fileName)     #   Notice the following arrays are MASKED arrays!!!
                fastLat = numpy.ma.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
                fastLon = numpy.ma.array(loaded_file.variables['g0_lon_2'][self.finalLon[0]:self.finalLon[-1]+1])
                cosValue = numpy.cos(fastLat*(3.14/180.))
                cosValue = cosValue.reshape(self.finalLat.size, 1)
                maskedIceFraction = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][int(monthNo),self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])*cosValue
                #   There are some cells in which the value of ice fraction is e+29 (i.e. regions covered with land) so we set them to zero.
                if landMask == True:
                    landMa = numpy.ma.masked_less(maskedIceFraction,2.).filled(0.)

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
                    
                    if landMask == True:
                        days = calendar.monthrange(int(fileName[4:8]),int(monthNo))[1]
                        ff2 = open('landMask'+str(days)+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+\
                    latSec,'w')
                        landMa.tofile(ff2,",")
                        ff2.close()
                    ff1 = open(nameToWrite+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(monthNo)][0:3]+'_'\
                    +str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+\
                    latSec,'w')
                    ff1.write(json.dumps(totalIceFraction))
                    ff1.close()
        return numpy.array(totalIceFraction)



##########################
#   REMEMBER TO CHECK FOR LEAP YEARS!!!
##########################
        
    def energyReader(self,save = True, outputName = None):
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        totalFinalIndex = numpy.array([])

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
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
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((EnergyAverage/float(len(self.yearRange))).tolist()))
            ff1.close()
            ff2.close()
        return EnergyAverage/float(len(self.yearRange)), totalFinalIndex



    def landEnergyReader(self, fName, save = True, highLowAndMonth = None, outputName = None):
        totalFinalIndex = numpy.array([])
        day1 = fName[8:10]
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor = numpy.fromfile(fName,float,-1,",")
        landCoor = landCoor.reshape(1,81,101).repeat(int(day1),0)   #   The number of days in this line does not count for LEAP years.
        maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 1))
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
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
            totalFinalIndex = numpy.append(totalFinalIndex,numpy.sum(localIndex))
            print fileName[5:12], '-----> DONE!!!'
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
            
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) == 1:
                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((EnergyAverage/float(len(self.yearRange))).tolist()))
            ff1.close()
            ff2.close()
        return EnergyAverage/float(len(self.yearRange)), totalFinalIndex


    def montecarlo(self, maskName = 'landMask30_10W40E_35N75N', content = 'DivQ', \
                    counts = 100, save = True, highLowAndMonth = None, outputName = None):
        #   years has to be a list of years in integer format.
        #   The self.days is correct for now but in the future if you are planning to change the code watch out self.days!!!
        import random
        yArray = numpy.arange(1979,2015)
        totalFinalIndex = numpy.array([])
        day1 = maskName[8:10]
        totalAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor = numpy.fromfile(maskName,float,-1,",")
        landCoor = landCoor.reshape(1,81,101).repeat(int(day1),0)   #   The number of days in this line does not count for LEAP years.
        maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 1))
#        fastLat = numpy.array(loaded_file.variables['lat_NH'][self.finalLat[0]:self.finalLat[-1]+1])
#        fastLon = numpy.array(loaded_file.variables['lon'][self.finalLon[0]:self.finalLon[-1]+1])
        fastLat = numpy.arange(75,34.5,-0.5)
        cosValue = numpy.cos(fastLat*(3.14/180.))
        cosValue = cosValue.reshape(fastLat.size,1)
        totalFiles = []        
        for totalYear in numpy.arange(1979,2015):
                totalFiles.append(content+'.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
        for totalName in totalFiles:
            loaded_file = netCDF4.Dataset(self.path+totalName)
            totalFastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            totalFastEnergy = numpy.ma.masked_where(maskOfEurope,totalFastEnergy).filled(0.)
            totalAverage = totalAverage + numpy.sum(totalFastEnergy,0)/self.days
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastEnergy = numpy.ma.masked_where(maskOfEurope,fastEnergy).filled(0.)
            EnergyAverage = EnergyAverage + numpy.sum(fastEnergy,0)/self.days
        
        subtractedSpecific = EnergyAverage - totalAverage
        for i in range(0,counts):
            print "Monte Carlo cycle: ", i
            mcFiles = []           
            years = random.sample(yArray,len(self.filex))
            mcEnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for ycount in years:
                mcFiles.append(content+'.'+str(ycount)+'.'+str(self.month[0])+'.nc')
            for j in mcFiles:
                loaded_file = netCDF4.Dataset(self.path+j)
                mcFastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                mcFastEnergy = numpy.ma.masked_where(maskOfEurope,mcFastEnergy).filled(0.)
                mcEnergyAverage = mcEnergyAverage + numpy.sum(mcFastEnergy,0)/self.days
            subtractedMc = mcEnergyAverage - totalAverage
            cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedSpecific),numpy.absolute(subtractedMc)).astype(float)
        finalDensity = cellDensity/counts

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
                nameToWrite2 = self.fileContent+'Montecarlo'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) == 1:
#                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
#                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            finalDensity.tofile(ff2,sep=",")
#            ff1.close()
            ff2.close()
        return finalDensity


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
        
        
    

def smoothing(inArray, length = 7, mode = 'wrap'):
    """This function applies a smoothing technique to a given region."""
    from scipy.ndimage import uniform_filter
    filt = numpy.zeros((61,720))                          #((81,101))
    fastLat = numpy.arange(89.5,29.5,-0.5)          #(90,59.5,-0.5)
#    fastLat = numpy.arange(75,34.5,-0.5)
    print '----->', inArray.shape, fastLat.shape
    cosValue = numpy.cos(fastLat*(3.14/180.))
    cosValue = cosValue.reshape(fastLat.size,1)
    numerator = uniform_filter((inArray*cosValue),size=length,origin=0,mode = mode)
    denominator = numerator.reshape(inArray.shape)/cosValue
    for ilat in range(len(fastLat)):
        lonStep = int(round(min([720.,14.]/numpy.cos(fastLat[ilat]*(3.14/180.)))))
        filt[ilat:ilat+1,:] = uniform_filter(denominator[ilat,:],size=lonStep,origin=0,mode = mode)
    return numerator

        
        
        
def europeMap(data,longitude,latitude,sphereProjection= False, figTitle = None,\
     plotType = 'contourf', montecarlo = False, store = False, name = None):

    lat1, lat2 = latitude
    latitude123 = numpy.arange(lat1,lat2+0.5,0.5)
    lon1, lon2 = longitude
    longitude123= numpy.arange(lon1,lon2+0.5,0.5)

    fig1,ax1 = plt.subplots(1,1,figsize=(13,13))

    map1 = Basemap(projection='cyl',llcrnrlon=lon1,llcrnrlat=lat1, \
    urcrnrlon=lon2,urcrnrlat=lat2, lon_0 = 30, lat_0 = 30,\
    rsphere=6371200., resolution = 'l', area_thresh=10000)

    if sphereProjection == True:
        map1 = Basemap(projection='npstere', boundinglat=60, llcrnrlon=lon1,llcrnrlat=lat1, \
        urcrnrlon=lon2, urcrnrlat=lat2, lon_0 = 270, round=True)
    print "Thi sis the shapes: ----->", data.shape, latitude123.shape, longitude123.shape
    lons,lats = numpy.meshgrid(longitude123,latitude123)
    x, y = map1(lons,lats)
    myData = data.reshape((latitude123.size,longitude123.size))
    myData = numpy.flipud(myData)   #   This line is important!!!
    cmap = plt.get_cmap('RdBu_r')
    if plotType == 'mesh':
        colormesh = map1.pcolormesh(x,y,myData,cmap = cmap)
        cb = map1.colorbar(colormesh, location='bottom', pad="5%")
    elif plotType == 'contourf':
        if montecarlo == True:
            v = numpy.arange(0, 1.05, 0.05)
            cmap = plt.cm.get_cmap('Blues')
            bounds = v
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            colorf = map1.contourf(x, y, myData,v,cmap=cmap,norm=norm)
            cset = map1.contour(x, y, myData,levels=[0.95,0.99],colors = ['y','r'], linewidths=[2.,2.],extent='both')
            plt.clabel(cset, fmt='%1.2f', inline=1, fontsize=10)
            cb = map1.colorbar(colorf,ticks = numpy.arange(0,1.2,0.2),location='bottom', pad="5%")
            cb.ax.set_xticklabels(['0.0','20%','40%','60%','80%','100%'])
            cb.cmap.set_under('w')
            cb.cmap.set_over('g')
            cb.set_label(r'Latent energy divergence (Wm-2)')
        elif montecarlo == False:
            v = numpy.arange(-200, 201, 20)
            cmap = plt.cm.get_cmap('RdBu_r')
            bounds = v
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            colorf = map1.contourf(x, y, myData,v,cmap=cmap,norm=norm,extend = 'both')
            cb = map1.colorbar(colorf, location='bottom', pad="5%")
            cb.cmap.set_under('m')
            cb.cmap.set_over('g')
            cb.set_label(r'Latent energy divergence (Wm-2)')

    if sphereProjection == True:
        map1.drawmeridians([-150,-120,-90,-60,-30,0,30,60,90,120,150,180], labels=[1,1,1,1], linewidth=1.0)
        map1.drawparallels(numpy.arange(60,81,10),labels=[1,1,1,1], linewidth=1.0)
    elif sphereProjection == False:
        map1.drawparallels(numpy.arange(lat1,lat2 + 1,10),labels=[1,0,0,0], linewidth=0.0)
        map1.drawmeridians(numpy.arange(lon1,lon2 + 1,10),labels=[0,0,0,1], linewidth=0.0)
    
    map1.drawcoastlines()
#    map1.bluemarble()
#    map1.etopo()
    map1.drawcountries()
    ax1.set_title(str(figTitle))
#    ax1.annotate('Ave = 5.01', xy = (10,40), xytext=(40,76))
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





