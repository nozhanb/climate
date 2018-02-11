import json
import numpy
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
        self.lon = ((180+longitude[0])*2,(180+longitude[1])*2)
        self.lat5 = latitude
        self.lon5 = longitude
        self.totalSize = int(((latitude[1]-latitude[0])*2 + 1)*((longitude[1]-longitude[0])*2 + 1))
        self.fileContent = fileContent
        self.filex = []
        self.filexWater = []
        self.filexIce = []
        self.year1, self.year2 = year[0],year[1]
        self.year = year
        self.days = calendar.monthrange(int(self.year[0]),int(self.month[0]))[1]
        self.monthName = calendar.month_name[int(self.month[0])][0:3]    # this line creates the output month to be used in output file's name.

        totalLon = numpy.arange(-180,180,0.5)
        totalLat = numpy.arange(90,29.5,-0.5)
        self.lonInitialIndex = numpy.where(totalLon==float(self.lon5[0]))
        self.lonEndIndex = numpy.where(totalLon==float(self.lon5[1]))
        self.latInitialIndex = numpy.where(totalLat==float(self.lat5[0]))
        self.latEndIndex = numpy.where(totalLat==float(self.lat5[1]))
        self.finalLon = numpy.arange(self.lonInitialIndex[0][0], self.lonEndIndex[0][0] + 1)
        self.finalLat = numpy.arange(self.latEndIndex[0][0], self.latInitialIndex[0][0] + 1)  #   Notice that latitude is ordered differently.
        
        if allYears == True:
            self.yearRange = range(self.year1, self.year2+1)        
        elif allYears == False:
            self.yearRange = year
        
        if self.fileContent  == 'SIC':
            for years in self.yearRange:
                self.filex.append('SIC.'+str(years)+'.nc')
        elif self.fileContent == 'DivQ':
            self.path = './dmDivQ_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivQ.'+str(years)+'.'+months+'.nc')
        elif self.fileContent == 'DivD':
            self.path = './dmDivD_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivD.'+str(years)+'.'+str(months)+'.nc')
        elif self.fileContent == 'Geopo':
            self.path = './data/geoPo/'
            for years in self.yearRange:
                self.filex.append('Z.'+str(years)+'.nc')
        elif self.fileContent == 'LiqWat' or self.fileContent == 'FrozWat':
            self.pathWater = './data/liqWat/'
            for years in self.yearRange:
                self.filexWater.append('LiqWatPath.'+str(years)+'.nc')
            self.pathIce = './data/frozWat/'
            for years in self.yearRange:
                self.filexIce.append('FrozWatPath.'+str(years)+'.nc')
#        elif self.fileContent == 'FrozWat':
#            self.path = './data/frozWat'
#            for years in self.yearRange:
#                self.filexIce.append('FrozWaterPath.'+str(years)+'.nc')





    def sicReader(self,save = True, outputName = None, landMask = False):
        ''' Something needs to be written here.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        totalIceFraction = numpy.array([])
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
                #   There are some cells in which the value of ice fraction is e+19 (i.e. regions covered with land) so we set them to zero.
                if landMask == True:
                    landMa = numpy.ma.masked_less(maskedIceFraction,1.1).filled(0.)

                iceFractionMaskedWithZero = numpy.ma.masked_greater(maskedIceFraction,1.1).filled(0.)
                numerator = numpy.sum(numpy.sum(iceFractionMaskedWithZero,0))
                denominator = fastLon.size*numpy.sum(cosValue,0)
                iceIndex = numerator/denominator
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = iceIndex
                totalIceFraction = numpy.append(totalIceFraction,iceIndex)
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
            ff1.write(json.dumps(totalIceFraction.tolist()))
            ff1.close()
        return numpy.array(totalIceFraction)

###############################################################################
###############################################################################

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

###############################################################################
###############################################################################

    def landEnergyReader(self, fName = 'totalLandMask', save = True, highLowAndMonth = None,\
                        geoLevel = None, outputName = None):
        totalFinalIndex = numpy.array([])
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor2 = numpy.fromfile(fName,float,-1,",").reshape(121,720)
        landCoor1 = landCoor2[self.latEndIndex[0][0]:self.latInitialIndex[0][0]+1,self.lonInitialIndex[0][0]:self.lonEndIndex[0][0]+1]
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            days = calendar.monthrange(int(fileName[5:9]),int(fileName[10:12]))[1]
            landCoor = landCoor1.reshape(1,self.finalLat.size,self.finalLon.size).repeat(days,0)
            maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 1))
            fastLat = numpy.array(loaded_file.variables['lat_NH'][self.finalLat[0]:self.finalLat[-1]+1])
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


###############################################################################
###############################################################################


    def geoPoReader(self, fName = 'totalLandMask', geoLevel = 0, save = True, \
                    montecarlo=False, counts=None, highLowAndMonth = None, outputName = None):

        geoPoAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastGeoPo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            geoPoAverage = geoPoAverage + fastGeoPo/9.8
            print fileName[2:6], '-----> DONE!!!'            
            
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalGeopo = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('Z.'+str(totalYear)+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastGeopo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalGeopo = totalGeopo + totalFastGeopo/9.8
            
            subtractedGeopo = geoPoAverage - totalGeopo
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcGeopoAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append('Z.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcFastGeopo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcGeopoAverage = mcGeopoAverage + mcFastGeopo/9.8
                subtractedMc = mcGeopoAverage - totalGeopo
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedGeopo),numpy.absolute(subtractedMc)).astype(float)
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
                nameToWrite2 = self.fileContent+'Average'
                nameToWrite3 = self.fileContent+'Montecarlo'   
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) == 1:
#                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
#               ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            finalDensity.tofile(ff3,sep=",")
#            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((geoPoAverage/float(len(self.yearRange))).tolist()))
#            ff1.close()
            ff2.close()
            ff3.close()
            if montecarlo == True:
                return geoPoAverage/float(len(self.yearRange)), finalDensity
            elif montecarlo == False:
                return geoPoAverage/float(len(self.yearRange))

###############################################################################
###############################################################################

    def cloudReader(self, fName = 'totalLandMask', waterIce = 'water', save = True, \
                    montecarlo=False, counts=None, highLowAndMonth = None, outputName = None):
        
        totalFinalIndex = numpy.array([])
        totalFinalIndexWater = numpy.array([])
        totalFinalIndexIce = numpy.array([])
        waterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        iceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor2 = numpy.fromfile(fName,float,-1,",").reshape(121,720)
        landCoor1 = landCoor2[self.latEndIndex[0][0]:self.latInitialIndex[0][0]+1,self.lonInitialIndex[0][0]:self.lonEndIndex[0][0]+1]
        
        for fileNameWater, fileNameIce in zip(self.filexWater,self.filexIce):
            loaded_file_water = netCDF4.Dataset(self.pathWater+fileNameWater)
            loaded_file_ice = netCDF4.Dataset(self.pathIce+fileNameIce)
            landCoor = landCoor1.reshape(self.finalLat.size,self.finalLon.size)
            maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 1))
            fastLat = numpy.array(loaded_file_water.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(self.finalLat.size,1)
            fastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastWater = numpy.ma.masked_where(maskOfEurope,fastWater).filled(0.)
            fastIce = numpy.ma.masked_where(maskOfEurope,fastIce).filled(0.)
            waterAverage = waterAverage + fastWater
            iceAverage = iceAverage + fastIce
            cloudAverage = waterAverage + iceAverage
            localIndexWater = numpy.array([])
            localIndexIce = numpy.array([])
            for i in range(self.finalLat.size):
                numeratorWater = numpy.sum(numpy.sum(fastWater[i:i+1,:]*cosValue[i:i+1,:],0))
                numeratorIce = numpy.sum(numpy.sum(fastIce[i:i+1,:]*cosValue[i:i+1,:],0))
                denominatorWater = numpy.count_nonzero(fastWater[i:i+1,:])*numpy.sum(cosValue[i:i+1,:])
                denominatorIce = numpy.count_nonzero(fastWater[i:i+1,:])*numpy.sum(cosValue[i:i+1,:])
                if denominatorWater == 0.:
                    localIndexWater = numpy.append(localIndexWater,0.)
                elif denominatorWater != 0.:
                    localIndexWater = numpy.append(localIndexWater,numeratorWater/denominatorWater)
                if denominatorIce == 0.:
                    localIndexIce = numpy.append(localIndexIce,0.)
                elif denominatorIce != 0.:
                    localIndexIce = numpy.append(localIndexIce,numeratorIce/denominatorIce)
            totalFinalIndexWater = numpy.append(totalFinalIndexWater,numpy.sum(localIndexWater))
            totalFinalIndexIce = numpy.append(totalFinalIndexIce,numpy.sum(localIndexIce))
            totalFinalIndex = totalFinalIndexWater + totalFinalIndexIce
            print fileNameIce[12:16], '-----> DONE!!!'
#            print totalFinalIndex
        print 'Average index (Kgm-2)----->', numpy.average(totalFinalIndex)
        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalWaterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalIceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFilesWater = []
            totalFilesIce = []
            for totalYear in numpy.arange(1979,2015):
                totalFilesWater.append('LiqWatPath.'+str(totalYear)+'.nc')
                totalFilesIce.append('FrozWatPath.'+str(totalYear)+'.nc')
            for totalNameWater, totalNameIce in zip(totalFilesWater, totalFilesIce):
                loaded_file_water = netCDF4.Dataset(self.pathWater+totalNameWater)
                loaded_file_ice = netCDF4.Dataset(self.pathIce+totalNameIce)
                totalFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalFastWater = numpy.ma.masked_where(maskOfEurope,totalFastWater).filled(0.)
                totalFastIce = numpy.ma.masked_where(maskOfEurope,totalFastIce).filled(0.)
                totalWaterAverage = totalWaterAverage + totalFastWater
                totalIceAverage = totalIceAverage + totalFastIce
                totalCloudAverage = totalWaterAverage + totalIceAverage
            
            subtractedCloud = cloudAverage - totalCloudAverage
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFilesWater = []
                mcFilesIce = []
                years = random.sample(yArray,len(self.filexWater))
                mcWaterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcIceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFilesWater.append('LiqWatPath.'+str(ycount)+'.nc')
                    mcFilesIce.append('FrozWatPath.'+str(ycount)+'.nc')
                for i,j in zip(mcFilesWater,mcFilesIce):
                    loaded_file_water = netCDF4.Dataset(self.pathWater+i)
                    loaded_file_ice = netCDF4.Dataset(self.pathIce+j)
                    mcFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcFastWater = numpy.ma.masked_where(maskOfEurope,mcFastWater).filled(0.)
                    mcFastIce = numpy.ma.masked_where(maskOfEurope,mcFastIce).filled(0.)
                    mcWaterAverage = mcWaterAverage + mcFastWater
                    mcIceAverage = mcIceAverage + mcFastIce
                    mcCloudAverage = mcWaterAverage + mcIceAverage
                
                subtractedMc = mcCloudAverage - totalCloudAverage
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedCloud),numpy.absolute(subtractedMc)).astype(float)
            finalDensity = cellDensity/counts        

        if montecarlo == False:
            return cloudAverage/float(len(self.yearRange))
        elif montecarlo == True:
            return cloudAverage/float(len(self.yearRange)), finalDensity
            
###############################################################################
###############################################################################

    def montecarlo(self, maskName = 'totalLandMask', content = 'DivQ', \
                    counts = 100, indexSig = False, save = True, \
                    highLowAndMonth = None, outputName = None):
        ''' This method calls for a montecarlo comparison simulation. '''
        #   The self.days is correct for now but in the future if you are planning to change the code watch out self.days!!!
        import random
        yArray = numpy.arange(1979,2015)
        monteCount = 0.
        totalFinalIndex = numpy.array([])
        normalFinalIndex= numpy.array([])
        totalAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        EnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor2 = numpy.fromfile(maskName,float,-1,",").reshape(121,720)  # 121 and 720 are the total grid point in latitude and longitude direction respectively.
        landCoor1 =landCoor2[self.latEndIndex[0][0]:self.latInitialIndex[0][0]+1,self.lonInitialIndex[0][0]:self.lonEndIndex[0][0]+1]
        fastLat = numpy.arange(self.lat5[1],self.lat5[0]-0.5,-0.5)
        cosValue = numpy.cos(fastLat*(3.14/180.))
        cosValue = cosValue.reshape(fastLat.size,1)
        totalFiles = []
        for totalYear in numpy.arange(1979,2015):
            totalFiles.append(content+'.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
        for totalName in totalFiles:
            loaded_file = netCDF4.Dataset(self.path+totalName)
            days = calendar.monthrange(int(totalName[5:9]),int(totalName[10:12]))[1]
            landCoor = landCoor1.reshape(1,self.finalLat.size,self.finalLon.size).repeat(days,0)
            maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 2.))  # It does not have to be 2 just any number bigger than 1
            totalFastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            totalFastEnergy = numpy.ma.masked_where(maskOfEurope,totalFastEnergy).filled(0.)
            totalAverage = totalAverage + numpy.sum(totalFastEnergy,0)/days
            if indexSig == True:
                cosValue = cosValue.reshape(1,fastLat.size,1)
                totalLocalIndex = numpy.array([])
                for i in range(self.finalLat.size):
                    numerator = numpy.sum(numpy.sum(totalFastEnergy[:,i:i+1,:]*cosValue[:,i:i+1,:],0)/days)
                    denominator = numpy.count_nonzero(totalFastEnergy[:,i:i+1,:])*cosValue[:,i:i+1,:]
                    if denominator == 0.:
                        totalLocalIndex = numpy.append(totalLocalIndex,0.)
                    elif denominator != 0.:
                        totalLocalIndex = numpy.append(totalLocalIndex,numerator/denominator)
                totalFinalIndex = numpy.append(totalFinalIndex,numpy.sum(totalLocalIndex))
        
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            days = calendar.monthrange(int(fileName[5:9]),int(fileName[10:12]))[1]
            landCoor = landCoor1.reshape(1,self.finalLat.size,self.finalLon.size).repeat(days,0)
            maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 2.))
            fastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastEnergy = numpy.ma.masked_where(maskOfEurope,fastEnergy).filled(0.)
            EnergyAverage = EnergyAverage + numpy.sum(fastEnergy,0)/days
            if indexSig == True:
                cosValue = cosValue.reshape(1,fastLat.size,1)
                normalLocalIndex = numpy.array([])
                for i in range(self.finalLat.size):
                    numerator = numpy.sum(numpy.sum(fastEnergy[:,i:i+1,:]*cosValue[:,i:i+1,:],0)/days)
                    denominator = numpy.count_nonzero(fastEnergy[:,i:i+1,:])*cosValue[:,i:i+1,:]
                    if denominator == 0.:
                        normalLocalIndex = numpy.append(normalLocalIndex,0.)
                    elif denominator != 0.:
                        normalLocalIndex = numpy.append(normalLocalIndex,numerator/denominator)
                normalFinalIndex = numpy.append(normalFinalIndex,numpy.sum(normalLocalIndex))

        if indexSig == True:
            totalFinalIndexAve = numpy.average(totalFinalIndex)
            normalFinalIndex = numpy.average(normalFinalIndex)
            finalValue = normalFinalIndex - totalFinalIndexAve
            
        subtractedSpecific = EnergyAverage - totalAverage
        
        yearCollector = numpy.array([])
        for i in range(0,counts):
            print "Monte Carlo cycle: ", i
            mcFiles = []
            monteFinalIndex = numpy.array([])
            years = random.sample(yArray,len(self.filex))
            mcEnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for ycount in years:
                mcFiles.append(content+'.'+str(ycount)+'.'+str(self.month[0])+'.nc')
            for j in mcFiles:
                loaded_file = netCDF4.Dataset(self.path+j)
                days = calendar.monthrange(int(j[5:9]),int(j[10:12]))[1]
                landCoor = landCoor1.reshape(1,self.finalLat.size,self.finalLon.size).repeat(days,0)
                maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 2.))
                mcFastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                mcFastEnergy = numpy.ma.masked_where(maskOfEurope,mcFastEnergy).filled(0.)
                mcEnergyAverage = mcEnergyAverage + numpy.sum(mcFastEnergy,0)/days
                if indexSig == True:
                    cosValue = cosValue.reshape(1,fastLat.size,1)
                    monteLocalIndex = numpy.array([])
                    for i in range(self.finalLat.size):
                        numerator = numpy.sum(numpy.sum(mcFastEnergy[:,i:i+1,:]*cosValue[:,i:i+1,:],0)/days)
                        denominator = numpy.count_nonzero(mcFastEnergy[:,i:i+1,:])*cosValue[:,i:i+1,:]
                        if denominator == 0.:
                            monteLocalIndex = numpy.append(monteLocalIndex,0.)
                        elif denominator != 0.:
                            monteLocalIndex = numpy.append(monteLocalIndex,numerator/denominator)
                    monteFinalIndex = numpy.append(monteFinalIndex,numpy.sum(monteLocalIndex))
            if indexSig == True:
                monteValue = numpy.average(monteFinalIndex) - totalFinalIndexAve
                if abs(finalValue) > abs(monteValue):
                    monteCount += 1.
                elif abs(finalValue) < abs(monteValue):
                   yearCollector = numpy.append(yearCollector,years)
            subtractedMc = mcEnergyAverage - totalAverage
            cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedSpecific),numpy.absolute(subtractedMc)).astype(float)
        finalDensity = cellDensity/counts
        if indexSig == True:
            indexSignificance = monteCount/counts
        
        uni, cou = numpy.unique(yearCollector, return_counts = True)
        
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
                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'Land'+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            finalDensity.tofile(ff2,sep=",")
            ff1.close()
            ff2.close()
        if indexSig == True:
            return finalDensity , indexSignificance, zip(uni,cou)
        elif indexSig == False:
            return finalDensity
            
###############################################################################
###############################################################################

def smoothing(inArray, length = 7, mode = 'wrap'):
    """This function applies a smoothing technique to a given region."""
    from scipy.ndimage import uniform_filter1d
    filt = numpy.zeros((101,720))
    fastLat = numpy.arange(80,29.5,-0.5)
    cosValue = numpy.cos(fastLat*(3.14/180.))
    cosValue = cosValue.reshape(fastLat.size,1)
    numerator = uniform_filter1d((inArray*cosValue),axis=0,size=length,origin=0,mode = mode)
    denominator = numerator.reshape(inArray.shape)/cosValue
    for ilat in range(len(fastLat)):
        lonStep = int(min([720.,(length*2.)/numpy.cos(fastLat[ilat]*(3.14/180.))]))
        filt[ilat:ilat+1,:] = uniform_filter1d(denominator[ilat,:],size=lonStep,origin=0,mode = mode)
    
    return filt

###############################################################################
###############################################################################

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
    print "This is the shapes: ----->", data.shape, latitude123.shape, longitude123.shape
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
            cset = map1.contour(x, y, myData,levels=[0.90,0.95],colors = ['y','r'], linewidths=[2.,2.],extent='both')
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
