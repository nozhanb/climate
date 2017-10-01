import numpy
import netCDF4
import matplotlib.pyplot as plt
#import numpy.ma as ma  #   I do not remember what it does! must be checked.
#numpy.set_printoptions(threshold=numpy.nan) #   if activated it prints out the numpy arrays in their entierty 
                                            #   regrdless of how big they are.

class index:

    def __init__(self,file_name,month,latitude,longitude, between = False):
        
#        self.file = file_name   #   must be a list of file names in string format. e.g. ['file1','file2','file3']
        self.month = month      #   month has to be given as a list of integers e.g. [1,5] which stands for Jan. and May.      
        self.lat = ((90-latitude[0])*2,(90-latitude[1])*2)    #   must be a tuple, float and in degree
        self.lon = ((abs(180)+longitude[0])*2,(abs(180)+longitude[1])*2) #   must be a tuple, float and in degree    
        self.lat5 = latitude
        self.lon5 = longitude
        if between == False:
            self.file = file_name
        elif between == True:
            self.file = []
            initial = int(file_name[0].split('SIC.')[1].split('.nc')[0])
            final = int(file_name[1].split('SIC.')[1].split('.nc')[0])
            for year in range(initial, final+1):
                self.file.append('SIC.'+str(year)+'.nc')
    def ice_reader(self):
        ''' This function reads in a .nc file and returns a dictionary of differnt indices for indicated months.'''
        self.output_dic = {}
        for fileName in self.file:
            for monthNo in self.month:
                loaded_file = netCDF4.Dataset(fileName)
                latiInverse = loaded_file.variables['g0_lat_1'][int(self.lat[1]):int(self.lat[0]) + 1]
                lati = numpy.flipud(latiInverse.reshape((latiInverse.size,1))).flatten('C')
                longi = loaded_file.variables['g0_lon_2'][int(self.lon[0]):int(self.lon[1]) + 1]

#                r = 6371    # this is the radius of the Earth.
                numerator = numpy.array([])
#                iceArea = numpy.array([])
                denominator = numpy.array([])
                cosval = numpy.array([])
                latCounter = numpy.arange(self.lat[0] + 1, self.lat[1], -1)
                lonCounter = numpy.arange(self.lon[0], self.lon[1] + 1, 1)
                counter1 = 0
                for i in lati:
                    lat2 = i + 0.5
                    counter2 = 0
                    for j in longi:
                        cosval = numpy.append(cosval,numpy.cos(lat2*(3.14/180)))
                        surat = numpy.array(loaded_file.variables['CI_GDS0_SFC_123'][monthNo,latCounter[counter1],lonCounter[counter2]]*(numpy.cos(lat2*(3.14/180))*0.5)*0.5)
                        makh = numpy.array((numpy.cos(lat2*(3.14/180))*0.5)*0.5)
                        if surat > 10000:
                            print "lat,long,surat: ------>", i,j,surat
                        numerator = numpy.append(numerator,numpy.array(surat)) 
                        denominator = numpy.append(denominator,makh)
           
                        if int(lat2) == 90:
                            numerator = numpy.append(numerator,0)
                            denominator = numpy.append(denominator,0)
                        counter2 += 1
                    counter1 += 1
                    
                sumnum = numpy.sum(numerator)
                sumdenom = numpy.sum(denominator)

#                print 'min, max', min(iceArea), max(iceArea)
                print 'Dic name --->', fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = float(sumnum/sumdenom)
#==============================================================================
#             print 'dss max and min:', max(dss),min(dss)
#             print 'cosval max and min:', max(cosval), min(cosval)
#             print 'lat max and min:', max(lati), min(lati)
#             print 'long max and min:', max(longi), min(longi)
#             print 'icy max and min:', max(icy), min(icy), icy.size
#             print max(areaTotal)
#==============================================================================
        return self.output_dic
        
def plotIndex(xdata, ydata):
    x_axis = xdata     #   just give the initial and final year in int and tuple form.
    y_axis = ydata     #   ydata is in dic format
    valueArray = numpy.array([])
    for key, value in y_axis.iteritems():
        valueArray = numpy.append(valueArray, value)
    
    yearArray = numpy.arange(x_axis[0], x_axis[1]+1) 
    fig1, ax1 = plt.subplots(1, figsize=(15,6))
 #   ax1.plot(yearArray, valueArray, color = 'k', linestyle ='--', linewidth = 3.0)
    ax1.plot(yearArray,valueArray,color = 'g', linestyle ='-', linewidth = 3.0)
    #   Aesthetic part
        #   xand y -axix adjustment
    ax1.set_xlim(x_axis[0], x_axis[1])
#    ax1.set_ylim(min(valueArray),max(valueArray))
    print '---->',type(fig1)
    plt.xlabel('Year', fontsize = 17)
    plt.ylabel('SSI Index', fontsize = 17)
    plt.title('Sea Surface Ice Index 1979-2016')

    ax1.grid('on')
   
    plt.show()















