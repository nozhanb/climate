import numpy
import netCDF4
import calendar
import collections  #   to make an ordered dictionary.
import matplotlib.pyplot as plt
#import numpy.ma as ma  #   I do not remember what it does! must be checked.
#numpy.set_printoptions(threshold=numpy.nan) #   if activated it prints out the numpy arrays in their entierty 
                                            #   regrdless of how big they are.

#   data = index.index(['SIC.1979.nc','SIC.2016.nc'],[1,2,3],(80,90),(0,30), between=True)


class index:

    def __init__(self,file_name,month,latitude,longitude, between = False):
        
        #  file_name must be a list of file names in string format. e.g. ['file1','file2','file3']
        #  month has to be given as a list of integers e.g. [1,5] which stands for Jan. and May.
        #  latitude and longitude must be tuples, float and in degree (0,90), from 0-degree to 90-degree.
        #  betwee: if True it will find the ice cover over all the years between two given years in the file_name.
        
        #   testsL Check for the last latitude at 90 and also at 0 and make sure the right value is generated.        
        
        self.month = month
        self.lat = ((90-latitude[0])*2,(90-latitude[1])*2)
        self.lon = ((abs(180)+longitude[0])*2,(abs(180)+longitude[1])*2)
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
        self.output_dic = collections.OrderedDict() #   Notice that it MUST be an Ordered Dictionary!!!
        for fileName in self.file:
            for monthNo in self.month:
                loaded_file = netCDF4.Dataset(fileName)
                latiInverse = loaded_file.variables['g0_lat_1'][int(self.lat[1]):int(self.lat[0]) + 1]
                lati = numpy.flipud(latiInverse.reshape((latiInverse.size,1))).flatten('C')
                longi = loaded_file.variables['g0_lon_2'][int(self.lon[0]):int(self.lon[1]) + 1]

                numerator = numpy.array([])
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

                print 'Dic name --->', fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'
                self.output_dic[fileName+'_'+str(monthNo)+'_'+'('+str(self.lat5[0])+','+str(self.lat5[1])+')'+'_'+'('+str(self.lon5[0])+','+str(self.lon5[1])+')'] = float(sumnum/sumdenom)
        return self.output_dic
        
        
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
#                first = ax1[row].XAxis()
 #               print '----->', first
                print type(ax1)
                print type(ax1[row])
            
            
        fig1.subplots_adjust(hspace=0)

    plt.xlabel('Year', fontsize = 17)
    plt.ylabel('SSI Index', fontsize = 17)
    plt.title('Sea Surface Ice Index 1979-2016')
    
    plt.show()















