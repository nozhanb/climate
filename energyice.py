import json
import time
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
        self.filexDry = []
        self.filexLate = []
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
            self.path = './data/sic/'
            for years in self.yearRange:
                self.filex.append('SIC.'+str(years)+'.nc')
        elif self.fileContent == 'DivQ':
            self.path = './data/energy/latent/dmDivQ_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivQ.'+str(years)+'.'+months+'.nc')
        elif self.fileContent == 'DivD':
            self.path = './data/energy/static/dmDivD_NH/'
            for years in self.yearRange:
                for months in self.month:
                    self.filex.append('DivD.'+str(years)+'.'+str(months)+'.nc')
        elif self.fileContent == 'Energies':
            self.pathDry = './data/energy/static/dmDivD_NH/'
            for years in self.yearRange:
                self.filexDry.append('DivD.'+str(years)+'.'+str(self.month[0])+'.nc')
            self.pathLate = './data/energy/latent/dmDivQ_NH/'
            for years in self.yearRange:
                self.filexLate.append('DivQ.'+str(years)+'.'+str(self.month[0])+'.nc')
        elif self.fileContent == 'Geopo':
            self.path = './data/geoPo/'
            for years in self.yearRange:
                self.filex.append('Z.'+str(years)+'.nc')
        elif self.fileContent == 'Temperature':
            self.path = './data/2mTemp/'
            for years in self.yearRange:
                self.filex.append('SAT.'+str(years)+'.nc')
        elif self.fileContent == 'Seatemp':
            self.path = './data/seaTemp/'
            for years in self.yearRange:
                self.filex.append('SST.'+str(years)+'.nc')
        elif self.fileContent == 'Seapressure':
            self.path = './data/seaPress/'
            for years in self.yearRange:
                self.filex.append('SP.'+str(years)+'.nc')
        elif self.fileContent == 'Cloud':
            self.pathWater = './data/liqWat/'
            for years in self.yearRange:
                self.filexWater.append('LiqWatPath.'+str(years)+'.nc')
            self.pathIce = './data/frozWat/'
            for years in self.yearRange:
                self.filexIce.append('FrozWatPath.'+str(years)+'.nc')
        elif self.fileContent == 'Flux':
            self.path = './data/'
        elif self.fileContent == 'Precipitation':
            self.path = './data/precip/'
            for years in self.yearRange:
                self.filex.append('Z.'+str(years)+'.nc')
        elif self.fileContent == 'Moisture':
            self.path = './data/sh/'
            for years in self.yearRange:
                self.filex.append('SH.'+str(years)+'.nc')

###############################################################################
###############################################################################


    def sicReader(self,save = True, outputName = None, landMask = False, highLowAndMonth=None):
        ''' Something needs to be written here.'''
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        totalIceFraction = numpy.array([])
        iceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        for fileName in self.filex:
            printName = fileName+'_'+str(self.month[0])
            loaded_file = netCDF4.Dataset(self.path+fileName)     #   Notice the following arrays are MASKED arrays!!!
            fastLat = numpy.ma.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(self.finalLat.size, 1)
            fastIce = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastIce = numpy.ma.masked_greater(fastIce,1.1).filled(0.)
            #   There are some cells in which the value of ice fraction is e+19 (i.e. regions covered with land) so we set them to zero.
            if landMask == True:
                landMa = numpy.ma.masked_less(fastIce,1.1).filled(0.)
                
            iceAverage = iceAverage + fastIce
#            print fileName[4:8], '-----> DONE!!!'            
            numerator = numpy.sum(numpy.sum(fastIce*cosValue,0))
            denominator = self.finalLon.size*numpy.sum(cosValue)
            iceIndex = numerator/denominator
#            iceFractionMaskedWithZero = numpy.ma.masked_greater(fastIce,1.1).filled(0.)
#==============================================================================
#             numerator = numpy.sum(numpy.sum(iceFractionMaskedWithZero,0))
#             denominator = self.finalLon.size*numpy.sum(cosValue)
#             iceIndex = numerator/denominator
#==============================================================================
#            self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = iceIndex
            totalIceFraction = numpy.append(totalIceFraction,iceIndex)
            print printName, '------> DONE!!!'
        iceAverage=iceAverage/float(len(self.yearRange))

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
            
            if landMask == True:
                days = calendar.monthrange(int(fileName[4:8]),int(self.month))[1]
                ff3 = open('landMask'+str(days)+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+\
            latSec,'w')
                landMa.tofile(ff3,",")
                ff3.close()
                
            if len(self.filex) == 1:
                ff1 = open(nameToWrite1+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
            elif len(self.filex) > 1:
                ff1 = open(nameToWrite1+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')

            ff1.write(json.dumps(totalIceFraction.tolist()))
            ff2.write(json.dumps((iceAverage.tolist())))
            ff1.close()
            ff2.close()
        return iceAverage, totalIceFraction
                
                
#==============================================================================
#             ff1 = open(nameToWrite+'_'+str(self.year1)+'_'+str(self.year2)+'_'+calendar.month_name[int(monthNo)][0:3]+'_'\
#             +str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+\
#             latSec,'w')
#             ff1.write(json.dumps(totalIceFraction.tolist()))
#==============================================================================

###############################################################################
###############################################################################

    def energyReader(self,montecarlo=False,counts=None, smoothing=False, length=None,\
                        smoothingCount=None,mode = 'wrap',save = True, outputName = None,\
                        highLowAndMonth = None):
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
        print 'Average index (Wm-2)-----> ', numpy.average(totalFinalIndex)
        EnergyAverage=EnergyAverage/float(len(self.yearRange))
    
            
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalEnergy = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append(self.fileContent+'.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalEnergy = totalEnergy + numpy.sum(totalFastEnergy,0)/days
            totalEnergy = totalEnergy/float(len(yArray))
            
            subtractedEnergy = EnergyAverage - totalEnergy
        
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcEnergyAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append(self.fileContent+'.'+str(ycount)+'.'+str(self.month[0])+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcfastEnergy = numpy.array(loaded_file.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcEnergyAverage = mcEnergyAverage + numpy.sum(mcfastEnergy,0)/days
                mcEnergyAverage=mcEnergyAverage/float(len(self.yearRange))
    
                subtractedMc = mcEnergyAverage - totalEnergy

                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedEnergy),numpy.absolute(subtractedMc)).astype(float)
            finalDensity = cellDensity/counts


        if smoothing == True:
            from scipy.ndimage import uniform_filter1d
            fastLat = numpy.arange(self.lat5[1],self.lat5[0]-0.5,-0.5)
            filt = numpy.zeros((EnergyAverage.shape))#(101,720)
            densityFilt = numpy.zeros((EnergyAverage.shape))
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(fastLat.size,1)
            theArea = EnergyAverage
            count1 = 0
            if montecarlo==True:
                count2 = 0
                theDensity = finalDensity
                while count2 < smoothingCount:
                    densityNumerator = uniform_filter1d((theDensity*cosValue),axis=0,size=length,origin=0,mode=mode)
#==============================================================================
#                     densityNumerator = uniform_filter1d(theDensity,axis=0,size=length,origin=0,mode=mode)
#==============================================================================
                    densityDenominator = densityNumerator/cosValue
#==============================================================================
#                     densityDenominator = densityNumerator
#==============================================================================
                    for ilat in range(len(fastLat)):
                        lonStep = int(min([720.,(length*2.)/cosValue[ilat]]))  #   numpy.cos(fastLat[ilat]*(3.14/180.))
                        densityFilt[ilat:ilat+1,:] = uniform_filter1d(densityDenominator[ilat,:],size=lonStep,origin=0,mode = mode)
#==============================================================================
#                         densityFilt[ilat:ilat+1,:] = uniform_filter1d(densityDenominator[ilat,:],size=length,origin=0,mode = mode)
#==============================================================================
                    theDensity = densityFilt
                    count2 += 1
            while count1 < smoothingCount:
                numerator = uniform_filter1d((theArea*cosValue),axis=0,size=length,origin=0,mode=mode)
                denominator = numerator/cosValue    #   numerator.reshape(theArea.shape)/cosValue                
                for ilat in range(len(fastLat)):
                    lonStep = int(min([720.,(length*2.)/cosValue[ilat]]))  #   numpy.cos(fastLat[ilat]*(3.14/180.))
                    filt[ilat:ilat+1,:] = uniform_filter1d(denominator[ilat,:],size=lonStep,origin=0,mode = mode)
                theArea = filt
                count1 += 1

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
                nameToWrite3 = self.fileContent+'Montecarlo'
                nameToWrite4 = self.fileContent+'Smoothed'
                nameToWrite5 = self.fileContent+'SmoothedMontecarlo'
            
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) > 1:
                ff1 = open(nameToWrite1+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
                if smoothing == True:
                    if montecarlo==True:
                        ff5 = open(nameToWrite5+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                        densityFilt.tofile(ff5,sep=",")
                        ff5.close()
                    elif montecarlo==False:
                        ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                        filt.tofile(ff4,sep=",")
                        ff4.close()


            ff1.write(json.dumps(totalFinalIndex.tolist()))
            ff2.write(json.dumps((EnergyAverage.tolist())))
            ff1.close()
            ff2.close()

        if montecarlo == True:
            return EnergyAverage, totalFinalIndex, finalDensity
        elif montecarlo == False:
            return EnergyAverage, totalFinalIndex
        

###############################################################################
###############################################################################

    def landEnergyReader(self, fName = 'totalLandMask', save = True, highLowAndMonth = None,\
                    outputName = None):
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
            totalFinalIndex = numpy.append(totalFinalIndex,numpy.sum(localIndex)) # Make sure this has to be sum not MEAN!!!
            print fileName[5:12], '-----> DONE!!!'
        EnergyAverage=EnergyAverage/float(len(self.yearRange))
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
            ff2.write(json.dumps((EnergyAverage.tolist())))
            ff1.close()
            ff2.close()
        return EnergyAverage, totalFinalIndex


###############################################################################
###############################################################################


    def energiesReader(self, save = True, montecarlo=False, counts=100,\
                    smoothing=True,smoothingCount=3,length=3, mode='wrap',\
                    highLowAndMonth = None, outputName = None):
        
        dryAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        lateAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        energiesAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        
        
        for fileNameDry, fileNameLate in zip(self.filexDry,self.filexLate):
            loaded_file_dry = netCDF4.Dataset(self.pathDry+fileNameDry)
            days = calendar.monthrange(int(fileNameDry[5:9]),int(fileNameDry[10:12]))[1]
            loaded_file_late = netCDF4.Dataset(self.pathLate+fileNameLate)
            fastDry  = numpy.array(loaded_file_dry.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastLate = numpy.array(loaded_file_late.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            dryAverage = numpy.sum(fastDry,0)/days
            lateAverage   = numpy.sum(fastLate,0)/days
            energiesAverage = energiesAverage + dryAverage + lateAverage
            print fileNameDry[5:12], '-----> DONE!!!'
        energiesAverage = energiesAverage/float(len(self.yearRange))    # this is the final average value!!!
        

        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalDryAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalLateAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalEnergiesAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFilesDry = []
            totalFilesLate = []
            for totalYear in numpy.arange(1979,2015):
                totalFilesDry.append('DivD.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
                totalFilesLate.append('DivQ.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
            for totalNameDry, totalNameLate in zip(totalFilesDry, totalFilesLate):
                loaded_file_dry = netCDF4.Dataset(self.pathDry+totalNameDry)
                loaded_file_late = netCDF4.Dataset(self.pathLate+totalNameLate)
                totalFastDry = numpy.array(loaded_file_dry.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalFastLate = numpy.array(loaded_file_late.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalDryAverage = numpy.sum(totalFastDry,0)/days
                totalLateAverage   = numpy.sum(totalFastLate,0)/days
                totalEnergiesAverage = totalEnergiesAverage + totalDryAverage + totalLateAverage
            totalEnergiesAverage=totalEnergiesAverage/float(len(yArray))
            subtractedEnergies = energiesAverage - totalEnergiesAverage
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFilesDry = []
                mcFilesLate = []
                years = random.sample(yArray,len(self.filexDry))
                mcDryAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcLateAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcEnergiesAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFilesDry.append('DivD.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
                    mcFilesLate.append('DivQ.'+str(totalYear)+'.'+str(self.month[0])+'.nc')
                for i,j in zip(mcFilesDry,mcFilesLate):
                    loaded_file_dry = netCDF4.Dataset(self.pathDry+i)
                    loaded_file_late = netCDF4.Dataset(self.pathLate+j)
                    mcFastDry = numpy.array(loaded_file_dry.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcFastLate = numpy.array(loaded_file_late.variables['Divergence'][:,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcDryAverage = numpy.sum(mcFastDry,0)/days
                    mcLateAverage   = numpy.sum(mcFastLate,0)/days
                    mcEnergiesAverage = mcEnergiesAverage + mcDryAverage + mcLateAverage
                mcEnergiesAverage = mcEnergiesAverage/float(len(self.yearRange))
                subtractedMc = mcEnergiesAverage - totalEnergiesAverage
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedEnergies),numpy.absolute(subtractedMc)).astype(float)
            finalDensity = cellDensity/counts        
        
        
        if smoothing == True:
            from scipy.ndimage import uniform_filter1d
            fastLat = numpy.arange(self.lat5[1],self.lat5[0]-0.5,-0.5)
            filt = numpy.zeros((energiesAverage.shape))#(101,720)
            densityFilt = numpy.zeros((energiesAverage.shape))
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(fastLat.size,1)
            theArea = energiesAverage
            count1 = 0
            if montecarlo==True:
                count2 = 0
                theDensity = finalDensity
                while count2 < smoothingCount:
                    densityNumerator = uniform_filter1d((theDensity*cosValue),axis=0,size=length,origin=0,mode=mode)
                    densityDenominator = densityNumerator/cosValue
                    for ilat in range(len(fastLat)):
                        lonStep = int(min([720.,(length*2.)/cosValue[ilat]]))  #   numpy.cos(fastLat[ilat]*(3.14/180.))
                        densityFilt[ilat:ilat+1,:] = uniform_filter1d(densityDenominator[ilat,:],size=lonStep,origin=0,mode = mode)
                    theDensity = densityFilt
                    count2 += 1        
            while count1 < smoothingCount:
                numerator = uniform_filter1d((theArea*cosValue),axis=0,size=length,origin=0,mode=mode)
                denominator = numerator/cosValue    #   numerator.reshape(theArea.shape)/cosValue                
                for ilat in range(len(fastLat)):
                    lonStep = int(min([720.,(length*2.)/cosValue[ilat]]))  #   numpy.cos(fastLat[ilat]*(3.14/180.))
                    filt[ilat:ilat+1,:] = uniform_filter1d(denominator[ilat,:],size=lonStep,origin=0,mode = mode)
                theArea = filt
                count1 += 1        

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
                nameToWrite4 = self.fileContent+'Smoothed'
                nameToWrite5 = self.fileContent+'SmoothedMontecarlo'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth

            if len(self.filexDry) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filexDry) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')    
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
                if smoothing == True:
                    if montecarlo==True:
                        ff5 = open(nameToWrite5+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                        densityFilt.tofile(ff5,sep=",")
                        ff5.close()
                    elif montecarlo==False:
                        ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                        filt.tofile(ff4,sep=",")
                        ff4.close()
                        
                        
            ff2.write(json.dumps((energiesAverage).tolist()))
            ff2.close()
            
        if montecarlo == True:
            return energiesAverage, finalDensity
        elif montecarlo == False:
            return energiesAverage
            

###############################################################################
###############################################################################


    def geoPoReader(self, fName = 'totalLandMask', geoLevel = 0, save = True, \
                    montecarlo=False, counts=None, highLowAndMonth = None, \
                    combine=False, outputName = None):

        geoPoAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        
        if combine == False:
            for fileName in self.filex:
                loaded_file = netCDF4.Dataset(self.path+fileName)
                fastGeoPo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                geoPoAverage = geoPoAverage + fastGeoPo/9.8
                print fileName[2:6], '-----> DONE!!!'            
            geoPoAverage=geoPoAverage/float(len(self.yearRange))
        elif combine == True:
            for fileName in self.filex:
                geoPoMonthAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                loaded_file = netCDF4.Dataset(self.path+fileName)
                for mons in self.month:
                    fastGeoPo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(mons)-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    geoPoMonthAverage = geoPoMonthAverage + fastGeoPo/9.8
                geoPoMonthAverage = geoPoMonthAverage/float(len(self.month))
                geoPoAverage = geoPoAverage + geoPoMonthAverage
                print fileName[2:6], '-----> DONE!!!'
            geoPoAverage=geoPoAverage/float(len(self.yearRange))
        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalGeopo = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('Z.'+str(totalYear)+'.nc')
            if combine == False:
                for totalName in totalFiles:
                    loaded_file = netCDF4.Dataset(self.path+totalName)
                    totalFastGeopo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalGeopo = totalGeopo + totalFastGeopo/9.8
                totalGeopo = totalGeopo/float(len(yArray))
            elif combine == True:
                for totalName in totalFiles:
                    totalMonthGeopo = numpy.zeros((self.finalLat.size,self.finalLon.size))
                    loaded_file = netCDF4.Dataset(self.path+totalName)
                    for mons in self.month:
                        totalFastGeopo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        totalMonthGeopo = totalMonthGeopo + totalFastGeopo/9.8
                    totalMonthGeopo = totalMonthGeopo/float(len(self.month))
                    totalGeopo = totalGeopo + totalMonthGeopo
                totalGeopo = totalGeopo/float(len(yArray))
            
            subtractedGeopo = geoPoAverage - totalGeopo
            
            
            if combine == False:
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
                    mcGeopoAverage=mcGeopoAverage/float(len(self.yearRange))
                    subtractedMc = mcGeopoAverage - totalGeopo
                    cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedGeopo),numpy.absolute(subtractedMc)).astype(float)
                finalDensity = cellDensity/counts
            elif combine == True:
                for i in range(0,counts):
                    print "Monte Carlo cycle: ", i
                    mcFiles = []
                    years = random.sample(yArray,len(self.filex))
                    mcGeopoAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                    for ycount in years:
                        mcFiles.append('Z.'+str(ycount)+'.nc')
                    for j in mcFiles:
                        mcGeopoMonthAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                        loaded_file = netCDF4.Dataset(self.path+j)
                        for mons in self.month:
                            mcFastGeopo = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(mons)-1,geoLevel,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                            mcGeopoMonthAverage = mcGeopoMonthAverage + mcFastGeopo/9.8
                        mcGeopoMonthAverage = mcGeopoMonthAverage/float(len(self.month))
                        mcGeopoAverage = mcGeopoAverage + mcGeopoMonthAverage
                    mcGeopoAverage=mcGeopoAverage/float(len(self.yearRange))
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
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filex) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            
            
            ff2.write(json.dumps((geoPoAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return geoPoAverage, finalDensity
        elif montecarlo == False:
            return geoPoAverage



###############################################################################
###############################################################################

    def geoThicknessReader(self, thickLevels=[11,26], montecarlo=False, counts=None):

        #   Remember the leves: 500 hPa --- 11, 850 hPa --- 20, 1000 hPa --- 26
        geoPoAverage1 = numpy.zeros((self.finalLat.size,self.finalLon.size))
        geoPoAverage2 = numpy.zeros((self.finalLat.size,self.finalLon.size))

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastGeoPo1 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[0],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            fastGeoPo2 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[1],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            geoPoAverage1 = geoPoAverage1 + fastGeoPo1/9.8
            geoPoAverage2 = geoPoAverage2 + fastGeoPo2/9.8
        geoPoAverage1=geoPoAverage1/float(len(self.yearRange))
        geoPoAverage2=geoPoAverage2/float(len(self.yearRange))
        fastLat = numpy.array(loaded_file.variables['g0_lat_2'][self.finalLat[0]:self.finalLat[-1]+1])
        cosValue = numpy.cos(fastLat*(3.14/180.))
        cosValue = cosValue.reshape(self.finalLat.size,1)
        thickIndex1 = numpy.array([])
        thickIndex2 = numpy.array([])
        numerator1   = numpy.sum(numpy.sum(geoPoAverage1*cosValue,0))
        denominator1 = geoPoAverage1.shape[1]*numpy.sum(cosValue)
        thickIndex1  = numpy.append(thickIndex1,numerator1/denominator1)
        numerator2   = numpy.sum(numpy.sum(geoPoAverage2*cosValue,0))
        denominator2 = geoPoAverage2.shape[1]*numpy.sum(cosValue)
        thickIndex2  = numpy.append(thickIndex2,numerator2/denominator2)
#        thickFinalIndex1  = numpy.append(thickFinalIndex1,numpy.sum(thickIndex1))
#        thickFinalIndex2  = numpy.append(thickFinalIndex2,numpy.sum(thickIndex2))
        thickness = thickIndex1 - thickIndex2

        totalGeoPoAverage1 = numpy.zeros((self.finalLat.size,self.finalLon.size))
        totalGeoPoAverage2 = numpy.zeros((self.finalLat.size,self.finalLon.size))
        totalThickFiles = []
        totalThickIndex1 = numpy.array([])
        totalThickIndex2 = numpy.array([])
        for totalThickYear in numpy.arange(1979,2015):
            totalThickFiles.append('Z.'+str(totalThickYear)+'.nc')
        for totalThickName in totalThickFiles:
            loaded_file = netCDF4.Dataset(self.path+totalThickName)
            totalThickFastGeopo1 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[0],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            totalThickFastGeopo2 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[1],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            totalGeoPoAverage1 = totalGeoPoAverage1 + totalThickFastGeopo1/9.8
            totalGeoPoAverage2 = totalGeoPoAverage2 + totalThickFastGeopo2/9.8
        totalGeoPoAverage1=totalGeoPoAverage1/float(len(totalThickFiles))
        totalGeoPoAverage2=totalGeoPoAverage2/float(len(totalThickFiles))

        totalNumerator1   = numpy.sum(numpy.sum(totalGeoPoAverage1*cosValue,0))
        totalDenominator1 = totalGeoPoAverage1.shape[1]*numpy.sum(cosValue)
        totalThickIndex1  = numpy.append(totalThickIndex1,totalNumerator1/totalDenominator1)
        totalNumerator2   = numpy.sum(numpy.sum(totalGeoPoAverage2*cosValue,0))
        totalDenominator2 = totalGeoPoAverage2.shape[1]*numpy.sum(cosValue)
        totalThickIndex2  = numpy.append(totalThickIndex2,totalNumerator2/totalDenominator2)
#        totalThickFinalIndex1  = numpy.append(totalThickFinalIndex1,numpy.sum(totalThickIndex1))
#        totalThickFinalIndex2  = numpy.append(totalThickFinalIndex2,numpy.sum(totalThickIndex2))
        
        totalThickness = totalThickIndex1 - totalThickIndex2        

        thicknessAnomaly = thickness - totalThickness

        if montecarlo==True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((1))
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcGeopoAverage1 = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcGeopoAverage2 = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcThickIndex1 = numpy.array([])
                mcThickIndex2 = numpy.array([])
                for ycount in years:
                    mcFiles.append('Z.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcFastGeopo1 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[0],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcFastGeopo2 = numpy.array(loaded_file.variables['Z_GDS0_ISBL_123'][int(self.month[0])-1,thickLevels[1],self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1]) 
                    mcGeopoAverage1 = mcGeopoAverage1 + mcFastGeopo1/9.8
                    mcGeopoAverage2 = mcGeopoAverage2 + mcFastGeopo2/9.8
                mcGeopoAverage1=mcGeopoAverage1/float(len(self.yearRange))
                mcGeopoAverage2=mcGeopoAverage2/float(len(self.yearRange))
                mcNumerator1   = numpy.sum(numpy.sum(mcGeopoAverage1*cosValue,0))    #cosValue[i:i+1,:]
                mcDenominator1 = mcGeopoAverage1.shape[1]*numpy.sum(cosValue)    #numpy.sum(cosValue[i:i+1,0]
                mcThickIndex1  = numpy.append(mcThickIndex1,mcNumerator1/mcDenominator1)
                mcNumerator2   = numpy.sum(numpy.sum(mcGeopoAverage2*cosValue,0))
                mcDenominator2 = mcGeopoAverage2.shape[1]*numpy.sum(cosValue)
                mcThickIndex2  = numpy.append(mcThickIndex2,mcNumerator2/mcDenominator2)
#                mcThickFinalIndex1  = numpy.append(mcThickFinalIndex1,numpy.sum(mcThickIndex1))
#                mcThickFinalIndex2  = numpy.append(mcThickFinalIndex2,numpy.sum(mcThickIndex2))
            
                mcThickness = mcThickIndex1 - mcThickIndex2
                
                subtractedMc = mcThickness - totalThickness
                cellDensity = cellDensity + numpy.greater(numpy.absolute(thicknessAnomaly),numpy.absolute(subtractedMc)).astype(float)
            finalDensity = cellDensity/counts
        
        if montecarlo == True:
            return thicknessAnomaly, finalDensity
        elif montecarlo == False:
            return thicknessAnomaly


###############################################################################
###############################################################################

    def temperatureReader(self, fName = 'totalLandMask', save = True, accumulation=False,\
                    montecarlo=False, counts=None, highLowAndMonth = None, outputName = None):

        tempAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        landCoor2 = numpy.fromfile(fName,float,-1,",").reshape(121,720)
        landCoor1 = landCoor2[self.latEndIndex[0][0]:self.latInitialIndex[0][0]+1,self.lonInitialIndex[0][0]:self.lonEndIndex[0][0]+1]
        landCoor = landCoor1.reshape(self.finalLat.size,self.finalLon.size)
        maskOfEurope = numpy.ma.getmask(numpy.ma.masked_less(landCoor, 1))

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastTemp = numpy.array(loaded_file.variables['2T_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            #fastTemp = numpy.ma.masked_where(maskOfEurope,fastTemp).filled(0.)
            tempAverage = tempAverage + fastTemp
            print fileName[4:8], '-----> DONE!!!'            
        tempAverage=tempAverage/float(len(self.yearRange))
        
        
        if accumulation==True:
            totalFiles2 = []
            accuFinalIndex=numpy.array([])
            for totalYear in self.yearRange:
                totalFiles2.append('SAT.'+str(totalYear)+'.nc')
            for totalName2 in totalFiles2:
                loaded_file = netCDF4.Dataset(self.path+totalName2)
                totalFastTemp2 = numpy.array(loaded_file.variables['2T_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                fastLat = numpy.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
                #totalFastTemp2 = numpy.ma.masked_where(maskOfEurope,totalFastTemp2).filled(0.)
                cosValue = numpy.cos(fastLat*(3.14/180.))
                cosValue = cosValue.reshape(self.finalLat.size,1)
                accuIndex = numpy.array([])
                for i in range(self.finalLat.size):
                    numerator = numpy.sum(numpy.sum(totalFastTemp2[i:i+1,:]*cosValue[i:i+1,:],0))
                    denominator = numpy.count_nonzero(totalFastTemp2[i:i+1,:])*numpy.sum(cosValue[i:i+1,0])
                    if denominator == 0.:
                        accuIndex = numpy.append(accuIndex,0.)
                    elif denominator != 0.:
                        accuIndex = numpy.append(accuIndex,numerator/denominator)
                accuFinalIndex = numpy.append(accuFinalIndex,numpy.mean(numpy.ma.masked_equal(accuIndex,0.))-273.15)
                print totalName2[4:8], '-----> Temperature accumulation DONE!!!'

        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalTemp = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('SAT.'+str(totalYear)+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastTemp = numpy.array(loaded_file.variables['2T_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                #totalFastTemp = numpy.ma.masked_where(maskOfEurope,totalFastTemp).filled(0.)
                totalTemp = totalTemp + totalFastTemp
            totalTemp = totalTemp/float(len(yArray))
            
            subtractedTemp = tempAverage - totalTemp
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mctempAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append('SAT.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcfastTemp = numpy.array(loaded_file.variables['2T_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mctempAverage = mctempAverage + mcfastTemp
                mctempAverage=mctempAverage/float(len(self.yearRange))    
                subtractedMc = mctempAverage - totalTemp
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedTemp),numpy.absolute(subtractedMc)).astype(float)
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
                nameToWrite4 = self.fileContent+'Index'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filex) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close() 
                if accumulation==True:
                    ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    accuFinalIndex.tofile(ff4,sep=",")
                    ff4.close()
                    
            
            ff2.write(json.dumps((tempAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return tempAverage, finalDensity
        elif montecarlo==False and accumulation==True:
            return tempAverage,accuFinalIndex
        elif montecarlo==False and accumulation==False:
            return tempAverage


###############################################################################
###############################################################################


    def seatempReader(self, fName = 'totalLandMask', save = True, montecarlo=False,\
                        accumulation=False,counts=None, highLowAndMonth = None,\
                        outputName = None):

        seatempAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastSeatemp = numpy.array(loaded_file.variables['SSTK_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            seatempAverage = seatempAverage + fastSeatemp
            print fileName[4:8], '-----> DONE!!!'            
        seatempAverage=seatempAverage/float(len(self.yearRange))
           
           
        if accumulation==True:
            totalFiles2 = []
            accuSeatempAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for totalYear in self.yearRange:
                totalFiles2.append('SST.'+str(totalYear)+'.nc')
            for totalName2 in totalFiles2:
                loaded_file = netCDF4.Dataset(self.path+totalName2)
                accuSeatemp = numpy.array(loaded_file.variables['SSTK_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                accuSeatempAverage = accuSeatempAverage + accuSeatemp
                accuSeatempAverage = numpy.ma.masked_greater(accuSeatempAverage,999.)
            accuSeatempAverage=accuSeatempAverage/float(len(self.yearRange))
            fastLat = numpy.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(fastLat*(3.14/180.))
            cosValue = cosValue.reshape(self.finalLat.size,1)
            accuIndex = numpy.array([])
            
            accNumerator1 = numpy.sum(numpy.sum(accuSeatempAverage*cosValue,0))
            accuDenominator1 = accuSeatempAverage.shape[1]*numpy.sum(cosValue)
            accuIndex = numpy.append(accuIndex,accNumerator1/accuDenominator1)

            print 'This is the average value of the SST: ', accuIndex           
           
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalSeatemp = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('SST.'+str(totalYear)+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastSeatemp = numpy.array(loaded_file.variables['SSTK_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalSeatemp = totalSeatemp + totalFastSeatemp
            totalSeatemp = totalSeatemp/float(len(yArray))
            
            subtractedTemp = seatempAverage - totalSeatemp
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcSeatempAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append('SST.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcfastSeatemp = numpy.array(loaded_file.variables['SSTK_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcSeatempAverage = mcSeatempAverage + mcfastSeatemp
                mcSeatempAverage=mcSeatempAverage/float(len(self.yearRange))    
                subtractedMc = mcSeatempAverage - totalSeatemp
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedTemp),numpy.absolute(subtractedMc)).astype(float)
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
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filex) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()

            
            ff2.write(json.dumps((seatempAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return seatempAverage, finalDensity
        elif montecarlo == False:
            return seatempAverage



###############################################################################
###############################################################################



    def seapressureReader(self, fName = 'totalLandMask', save = True, montecarlo=False,\
                        counts=None, highLowAndMonth = None, outputName = None):

        seaPressAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))

        for fileName in self.filex:
            loaded_file = netCDF4.Dataset(self.path+fileName)
            fastSeaPress = numpy.array(loaded_file.variables['SP_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            seaPressAverage = seaPressAverage + fastSeaPress/100.
            print fileName[3:8], '-----> DONE!!!'            
        seaPressAverage=seaPressAverage/float(len(self.yearRange))
            
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalSeaPress = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('SP.'+str(totalYear)+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastSeaPress = numpy.array(loaded_file.variables['SP_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalSeaPress = totalSeaPress + totalFastSeaPress/100.
            totalSeaPress = totalSeaPress/float(len(yArray))
            
            subtractedTemp = seaPressAverage - totalSeaPress
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcSeaPressAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append('SP.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcfastSeaPress = numpy.array(loaded_file.variables['SP_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcSeaPressAverage = mcSeaPressAverage + mcfastSeaPress/100.
                mcSeaPressAverage = mcSeaPressAverage/float(len(self.yearRange))    
                subtractedMc = mcSeaPressAverage - totalSeaPress
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedTemp),numpy.absolute(subtractedMc)).astype(float)
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
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filex) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            
            ff2.write(json.dumps((seaPressAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return seaPressAverage, finalDensity
        elif montecarlo == False:
            return seaPressAverage



###############################################################################
###############################################################################

    def fluxReader(self, fName = 'totalLandMask', fluxType = ['sws','lws','ssh','slh'],\
                    save = True, montecarlo=False, counts=None, accumulation=False,\
                    highLowAndMonth = None,outputName = None):
        
        filex2 = []
        fluxPair = {'sws':'SSR_GDS0_SFC_120','lws':'STR_GDS0_SFC_120','ssh':'SSHF_GDS0_SFC_120','slh':'SLHF_GDS0_SFC_120'}
        fluxAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        for ftype in fluxType:
            filex2.append(ftype)
        for year in self.yearRange:
            energyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
#            days = calendar.monthrange(int(year),int(self.month[0]))[1]
            secToDays = 86400.
            for fluxName in filex2:
                loaded_file = netCDF4.Dataset(self.path+'flux/'+fluxName+'/'+fluxName.upper()+'.'+str(year)+'.nc')
                fastFlux = numpy.array(loaded_file.variables[fluxPair[fluxName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                energyBudget = energyBudget + fastFlux/secToDays
            fluxAverage = fluxAverage + energyBudget
            print year, '-----> DONE!!!'
        fluxAverage=fluxAverage/float(len(self.yearRange))
        
        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFluxAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for year in yArray:
                totalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
#                days = calendar.monthrange(int(year),int(self.month[0]))[1]
                secToDays = 86400.
                for fluxName in filex2:
                    loaded_file = netCDF4.Dataset(self.path+'flux/'+fluxName+'/'+fluxName.upper()+'.'+str(year)+'.nc')
                    totalFastFlux = numpy.array(loaded_file.variables[fluxPair[fluxName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalEnergyBudget = totalEnergyBudget + totalFastFlux/secToDays
                totalFluxAverage = totalFluxAverage + totalEnergyBudget
            totalFluxAverage=totalFluxAverage/float(len(yArray))
            subtractedFlux = fluxAverage - totalFluxAverage
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                years = random.sample(yArray,len(self.yearRange))
                mcFluxAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for year in years:  #   Notice that year and years are different!!!
                    mcTotalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
#                    days = calendar.monthrange(int(year),int(self.month[0]))[1]
                    secToDays = 86400.
                    for fluxName in filex2:
                        loaded_file = netCDF4.Dataset(self.path+'flux/'+fluxName+'/'+fluxName.upper()+'.'+str(year)+'.nc')
                        mcTotalFastFlux = numpy.array(loaded_file.variables[fluxPair[fluxName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        mcTotalEnergyBudget = mcTotalEnergyBudget + mcTotalFastFlux/secToDays
                    mcFluxAverage = mcFluxAverage + mcTotalEnergyBudget
                mcFluxAverage = mcFluxAverage/float(len(self.yearRange))
                mcSubtractedFlux = mcFluxAverage - totalFluxAverage
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedFlux),numpy.absolute(mcSubtractedFlux)).astype(float)
            finalDensity = cellDensity/counts
            

        if accumulation == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFluxAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for year in yArray:
                totalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
#                days = calendar.monthrange(int(year),int(self.month[0]))[1]
                secToDays = 86400.
                for fluxName in filex2:
                    loaded_file = netCDF4.Dataset(self.path+'flux/'+fluxName+'/'+fluxName.upper()+'.'+str(year)+'.nc')
                    totalFastFlux = numpy.array(loaded_file.variables[fluxPair[fluxName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalEnergyBudget = totalEnergyBudget + totalFastFlux/secToDays
                totalFluxAverage = totalFluxAverage + totalEnergyBudget
            totalFluxAverage=totalFluxAverage/float(len(yArray))
            subtractedFlux = fluxAverage - totalFluxAverage
            
            totalAccuFlux = numpy.array([])
            latAccuFlux = numpy.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(latAccuFlux*(3.14/180.)).reshape(latAccuFlux.size,1)
            accuFluxNume = numpy.sum(numpy.sum(subtractedFlux*cosValue,0))
            accuFluxDeno = self.finalLon.size*numpy.sum(cosValue)
            accuFluxRatio = accuFluxNume/accuFluxDeno
            totalAccuFlux = numpy.append(totalAccuFlux,accuFluxRatio)


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
                nameToWrite4 = self.fileContent+'Accumulated'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth

            if len(self.yearRange) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.yearRange) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
                if accumulation == True:
                    ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    totalAccuFlux.tofile(ff4,sep=",")
                    ff4.close()
            
            ff2.write(json.dumps((fluxAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return fluxAverage, finalDensity
        elif montecarlo == False:
            return fluxAverage
###############################################################################
###############################################################################

    def precipReader(self, fName = 'totalLandMask', precipType = ['LargeScalePrec','ConvectivePrec'],\
                    save = True, montecarlo=False, counts=None, accumulation=False,\
                    highLowAndMonth = None,outputName = None):
        
        filex3 = []
        precipPair = {'LargeScalePrec':'LSP_GDS0_SFC_120','ConvectivePrec':'CP_GDS0_SFC_120','LargeScaleSnow':'LSF_GDS0_SFC_120','ConvectiveSnow':'CSF_GDS0_SFC_120'}
        precipAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        for ftype in precipType:
            filex3.append(ftype)
        for year in self.yearRange:
            energyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for precipName in filex3:
                loaded_file = netCDF4.Dataset(self.path+precipName+'/'+precipName+'.'+str(year)+'.nc')
                fastprecip = numpy.array(loaded_file.variables[precipPair[precipName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                energyBudget = energyBudget + fastprecip*1000.  #   Multiplied by 1000 to convert from meter to mm.
            precipAverage = precipAverage + energyBudget
            print year, '-----> DONE!!!'
        precipAverage=precipAverage/float(len(self.yearRange))
        
        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalprecipAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for year in yArray:
                totalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for precipName in filex3:
                    loaded_file = netCDF4.Dataset(self.path+precipName+'/'+precipName+'.'+str(year)+'.nc')
                    totalFastprecip = numpy.array(loaded_file.variables[precipPair[precipName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalEnergyBudget = totalEnergyBudget + totalFastprecip*1000.
                totalprecipAverage = totalprecipAverage + totalEnergyBudget
            totalprecipAverage=totalprecipAverage/float(len(yArray))
            
            subtractedprecip = precipAverage - totalprecipAverage
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                years = random.sample(yArray,len(self.yearRange))
                mcprecipAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for year in years:  #   Notice that year and years are different!!!
                    mcTotalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
                    for precipName in filex3:
                        loaded_file = netCDF4.Dataset(self.path+precipName+'/'+precipName+'.'+str(year)+'.nc')
                        mcTotalFastprecip = numpy.array(loaded_file.variables[precipPair[precipName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        mcTotalEnergyBudget = mcTotalEnergyBudget + mcTotalFastprecip*1000.
                    mcprecipAverage = mcprecipAverage + mcTotalEnergyBudget
                mcprecipAverage = mcprecipAverage/float(len(self.yearRange))
                
                mcSubtractedprecip = mcprecipAverage - totalprecipAverage
                
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedprecip),numpy.absolute(mcSubtractedprecip)).astype(float)
            finalDensity = cellDensity/counts
            

        if accumulation == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalprecipAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            for year in yArray:
                totalEnergyBudget = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for precipName in filex3:
                    loaded_file = netCDF4.Dataset(self.path+precipName+'/'+precipName+'.'+str(year)+'.nc')
                    totalFastprecip = numpy.array(loaded_file.variables[precipPair[precipName]][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalEnergyBudget = totalEnergyBudget + totalFastprecip*1000.
                totalprecipAverage = totalprecipAverage + totalEnergyBudget
            totalprecipAverage = totalprecipAverage/float(len(yArray))
            subtractedprecip = precipAverage - totalprecipAverage
            
            totalAccuprecip = numpy.array([])
            latAccuprecip = numpy.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
            cosValue = numpy.cos(latAccuprecip*(3.14/180.)).reshape(latAccuprecip.size,1)
            accuprecipNume = numpy.sum(numpy.sum(subtractedprecip*cosValue,0))
            accuprecipDeno = self.finalLon.size*numpy.sum(cosValue)
            accuprecipRatio = accuprecipNume/accuprecipDeno
            totalAccuprecip = numpy.append(totalAccuprecip,accuprecipRatio)


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
                nameToWrite4 = self.fileContent+'Accumulated'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth

            if len(self.yearRange) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.yearRange) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
                if accumulation == True:
                    ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    totalAccuprecip.tofile(ff4,sep=",")
                    ff4.close()
            
            ff2.write(json.dumps((precipAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return precipAverage, finalDensity
        elif montecarlo == False:
            return precipAverage

###############################################################################
###############################################################################

    def cloudReader(self, fName = 'totalLandMask', save = True, montecarlo=False, \
                    counts=100, highLowAndMonth = None, outputName = None):
        
        waterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        iceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
        cloudAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
 
        for fileNameWater, fileNameIce in zip(self.filexWater,self.filexIce):
            loaded_file_water = netCDF4.Dataset(self.pathWater+fileNameWater)
            loaded_file_ice = netCDF4.Dataset(self.pathIce+fileNameIce)
            if fileNameWater == 'LiqWatPath.2013.nc':
                if int(self.month[0]) in [1,2,3,4,5,6,7,8,9]:
                    fastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    fastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                else:
                    fastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    fastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            else:
                    fastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    fastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])

            waterAverage = fastWater
            iceAverage   = fastIce
            cloudAverage = cloudAverage + waterAverage + iceAverage

            print fileNameIce[12:16], '-----> DONE!!!'
        cloudAverage=cloudAverage/float(len(self.yearRange))    # this is the final average value!!!
        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalWaterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalIceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalCloudAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFilesWater = []
            totalFilesIce = []
            for totalYear in numpy.arange(1979,2015):
                totalFilesWater.append('LiqWatPath.'+str(totalYear)+'.nc')
                totalFilesIce.append('FrozWatPath.'+str(totalYear)+'.nc')
            for totalNameWater, totalNameIce in zip(totalFilesWater, totalFilesIce):
                loaded_file_water = netCDF4.Dataset(self.pathWater+totalNameWater)
                loaded_file_ice = netCDF4.Dataset(self.pathIce+totalNameIce)
                if totalNameWater == 'LiqWatPath.2013.nc':
                    if int(self.month[0]) in [1,2,3,4,5,6,7,8,9]:
                        totalFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        totalFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    else:
                        totalFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        totalFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                else:
                    totalFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    totalFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])

                totalWaterAverage = totalFastWater
                totalIceAverage   = totalFastIce
                totalCloudAverage = totalCloudAverage + totalWaterAverage + totalIceAverage
            totalCloudAverage=totalCloudAverage/float(len(yArray))
            subtractedCloud = cloudAverage - totalCloudAverage
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFilesWater = []
                mcFilesIce = []
                years = random.sample(yArray,len(self.filexWater))
                mcWaterAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcIceAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                mcCloudAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFilesWater.append('LiqWatPath.'+str(ycount)+'.nc')
                    mcFilesIce.append('FrozWatPath.'+str(ycount)+'.nc')
                for i,j in zip(mcFilesWater,mcFilesIce):
                    loaded_file_water = netCDF4.Dataset(self.pathWater+i)
                    loaded_file_ice = netCDF4.Dataset(self.pathIce+j)
                    if i == 'LiqWatPath.2013.nc':
                        if int(self.month[0]) in [1,2,3,4,5,6,7,8,9]:
                            mcFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                            mcFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123_1'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        else:
                            mcFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                            mcFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-10,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    else:
                        mcFastWater = numpy.array(loaded_file_water.variables['VITCLCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                        mcFastIce = numpy.array(loaded_file_ice.variables['VITCFCW_GDS0_EATM_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    
                    mcWaterAverage = mcFastWater
                    mcIceAverage   = mcFastIce
                    mcCloudAverage = mcCloudAverage + mcWaterAverage + mcIceAverage
                mcCloudAverage = mcCloudAverage/float(len(self.yearRange))
                subtractedMc = mcCloudAverage - totalCloudAverage
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedCloud),numpy.absolute(subtractedMc)).astype(float)
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

            if len(self.filexWater) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filexWater) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')    
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
                    
            ff2.write(json.dumps((cloudAverage).tolist()))
            ff2.close()
            
        if montecarlo == True:
            return cloudAverage, finalDensity
        elif montecarlo == False:
            return cloudAverage
            
            
###############################################################################
###############################################################################


    def moistureReader(self, geoLevel=15, save = True, accumulation=False,\
                    montecarlo=False, counts=None, highLowAndMonth = None, outputName = None):
        #   Pressure levels (starting from 0): 8=500mb, 15=850, 18=1000mb
        moistAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))

        for fileName in self.filex:
            print self.path+fileName
            loaded_file  = netCDF4.Dataset(self.path+fileName)
            fastmoist    = numpy.array(loaded_file.variables['Q_GDS0_ISBL_123'][int(self.month[0])-1,int(geoLevel),self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
            moistAverage = moistAverage + fastmoist
            print fileName[4:8], '-----> DONE!!!'            
        moistAverage = moistAverage/float(len(self.yearRange))
        
        
#==============================================================================
#==============================================================================
# #         if accumulation==True:
# #             totalFiles2 = []
# #             accuFinalIndex=numpy.array([])
# #             for totalYear in self.yearRange:
# #                 totalFiles2.append('SAT.'+str(totalYear)+'.nc')
# #             for totalName2 in totalFiles2:
# #                 loaded_file = netCDF4.Dataset(self.path+totalName2)
# #                 totalFastmoist2 = numpy.array(loaded_file.variables['2T_GDS0_SFC_123'][int(self.month[0])-1,self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
# #                 fastLat = numpy.array(loaded_file.variables['g0_lat_1'][self.finalLat[0]:self.finalLat[-1]+1])
# #                 cosValue = numpy.cos(fastLat*(3.14/180.))
# #                 cosValue = cosValue.reshape(self.finalLat.size,1)
# #                 accuIndex = numpy.array([])
# #                 for i in range(self.finalLat.size):
# #                     numerator = numpy.sum(numpy.sum(totalFastmoist2[i:i+1,:]*cosValue[i:i+1,:],0))
# #                     denominator = numpy.count_nonzero(totalFastmoist2[i:i+1,:])*numpy.sum(cosValue[i:i+1,0])
# #                     if denominator == 0.:
# #                         accuIndex = numpy.append(accuIndex,0.)
# #                     elif denominator != 0.:
# #                         accuIndex = numpy.append(accuIndex,numerator/denominator)
# #                 accuFinalIndex = numpy.append(accuFinalIndex,numpy.mean(numpy.ma.masked_equal(accuIndex,0.))-273.15)
# #                 print totalName2[4:8], '-----> moisterature accumulation DONE!!!'
# #==============================================================================
#==============================================================================

        
        if montecarlo == True:
            import random
            yArray = numpy.arange(1979,2015)
            cellDensity = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalmoist = numpy.zeros((self.finalLat.size,self.finalLon.size))
            totalFiles = []
            for totalYear in numpy.arange(1979,2015):
                totalFiles.append('SH.'+str(totalYear)+'.nc')
            for totalName in totalFiles:
                loaded_file = netCDF4.Dataset(self.path+totalName)
                totalFastmoist = numpy.array(loaded_file.variables['Q_GDS0_ISBL_123'][int(self.month[0])-1,int(geoLevel),self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                totalmoist = totalmoist + totalFastmoist                
            totalmoist = totalmoist/float(len(yArray))
            
            subtractedmoist = moistAverage - totalmoist
            
            for i in range(0,counts):
                print "Monte Carlo cycle: ", i
                mcFiles = []
                years = random.sample(yArray,len(self.filex))
                mcmoistAverage = numpy.zeros((self.finalLat.size,self.finalLon.size))
                for ycount in years:
                    mcFiles.append('SH.'+str(ycount)+'.nc')
                for j in mcFiles:
                    loaded_file = netCDF4.Dataset(self.path+j)
                    mcfastmoist = numpy.array(loaded_file.variables['Q_GDS0_ISBL_123'][int(self.month[0])-1,int(geoLevel),self.finalLat[0]:self.finalLat[-1]+1,self.finalLon[0]:self.finalLon[-1]+1])
                    mcmoistAverage = mcmoistAverage + mcfastmoist
                mcmoistAverage=mcmoistAverage/float(len(self.yearRange))
                
                subtractedMc = mcmoistAverage - totalmoist
                
                cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedmoist),numpy.absolute(subtractedMc)).astype(float)
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
                nameToWrite4 = self.fileContent+'Index'
                
            if highLowAndMonth is None:
                baseM = '---'
            elif highLowAndMonth is not None:
                baseM = highLowAndMonth
                
            if len(self.filex) == 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close()
            elif len(self.filex) > 1:
                ff2 = open(nameToWrite2+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                if montecarlo == True:
                    ff3 = open(nameToWrite3+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    finalDensity.tofile(ff3,sep=",")
                    ff3.close() 
                if accumulation==True:
                    ff4 = open(nameToWrite4+baseM+'_'+str(self.year1)+'_'+str(self.year2)+'_'+self.monthName+'_'+str(abs(self.lon5[0]))+lonFir+str(abs(self.lon5[1]))+lonSec+'_'+str(abs(self.lat5[0]))+latFir+str(abs(self.lat5[1]))+latSec,'w')
                    accuFinalIndex.tofile(ff4,sep=",")
                    ff4.close()
                    
            
            ff2.write(json.dumps((moistAverage).tolist()))
            ff2.close()
        if montecarlo == True:
            return moistAverage, finalDensity
        elif montecarlo==False and accumulation==True:
            return moistAverage,accuFinalIndex
        elif montecarlo==False and accumulation==False:
            return moistAverage


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
        totalAverage = totalAverage/float(len(yArray))
        
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
                
        EnergyAverage = EnergyAverage/float(len(self.yearRange))
        if indexSig == True:
            totalFinalIndexAve = numpy.average(totalFinalIndex)
            normalFinalIndex = numpy.average(normalFinalIndex)
            finalValue = normalFinalIndex - totalFinalIndexAve
            
        subtractedSpecific = EnergyAverage - totalAverage
        
#        yearCollector = numpy.array([])
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
            mcEnergyAverage = mcEnergyAverage/float(len(self.yearRange))
            if indexSig == True:
                monteValue = numpy.average(monteFinalIndex) - totalFinalIndexAve
                if abs(finalValue) > abs(monteValue):
                    monteCount += 1.
#                elif abs(finalValue) < abs(monteValue):
#                   yearCollector = numpy.append(yearCollector,years)
            subtractedMc = mcEnergyAverage - totalAverage
            cellDensity = cellDensity + numpy.greater(numpy.absolute(subtractedSpecific),numpy.absolute(subtractedMc)).astype(float)
        finalDensity = cellDensity/counts
        if indexSig == True:
            indexSignificance = monteCount/counts
        
#        uni, cou = numpy.unique(yearCollector, return_counts = True)
        
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
            return finalDensity , indexSignificance #, zip(uni,cou)
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
     plotType = 'contourf', montecarlo = False, save = False, name = None):

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
    if save == True:
        plt.savefig('/home/nba035/plot/'+name+'.eps',dpi = 80, format = 'eps')

    plt.show()


