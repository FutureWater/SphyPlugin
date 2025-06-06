from __future__ import print_function

# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
# Copyright (C) 2013-2019  FutureWater
# Email: sphy@futurewater.nl
#
# Authors (alphabetical order):
# P. Droogers, J. Eekhout, W. Immerzeel, S. Khanal, A. Lutz, G. Simons, W. Terink
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from builtins import object
#-Authorship information-###################################################################
__author__ = "FutureWater"
__copyright__ = "FutureWater"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "sphy@futurewater.nl"
__date__ ='1 February 2025'
############################################################################################
############################################################################################

import datetime, subprocess, os, glob, csv, math
from osgeo import osr

#-Class that defines the processing of the meteorological forcings
class processForcing(object):
    def __init__(self, resultsdir, t_srs, resolution, extent, startdate, enddate, \
            textlog, progbar, procsteps):
        #-Create forcing directory in results directory if it does not exist
        self.outdir = os.path.join(resultsdir, 'forcing/')
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        else: #-remove old forcing files
            f = glob.glob(self.outdir + '*' )
            for fi in f:
                os.remove(fi)
            
        #-Directory for temporary files
        self.tempdir = resultsdir
        #-Set the general settings
        self.t_srs = t_srs
        self.xMin = str(extent[0])
        self.yMin = str(extent[1])
        self.xMax = str(extent[2])
        self.yMax = str(extent[3])
        self.t_res = str(resolution)
        self.startDate = startdate
        self.endDate = enddate
        self.timeSteps = ((enddate-startdate).days + 1)
        #-Initialize the database settings
        self.dbSource = None
        self.dbTs = None
        self.dbSrs = None
        self.dbFormat = None
        self.precDBPath = None
        self.tavgDBPath = None
        self.tmaxDBPath = None
        self.tminDBPath = None
        #-Initialize the model and database DEMs
        self.modelDem = None
        self.dbDem = None
        #-Initialize the CSV file settings (if user stations will be used)
        self.precLocFile = None
        self.precDataFile = None
        self.tempLocFile = None
        self.tempDataFile = None
        #-Log
        self.textLog = textlog
        #-progressbar
        self.progBar = progbar
        self.procSteps = procsteps
        self.counter = 0.

        
    #-Create precipitation forcing based on the database    
    def createPrecDB(self):
        #-If the database is from WFDEI (Watch forcing)
        if self.dbSource == 'WFDEI':
            self.textLog.append('Processing precipitation from ' + self.dbSource + ' database...\n')
            for i in range(0, self.timeSteps):
                daynr = i+1 # required for pcraster extension
                curdate = self.startDate + datetime.timedelta(days=i)
                year = str(curdate.year)
                month = curdate.month
                if month<10:
                    month = '0'+str(month)
                else:
                    month = str(month)
                day = curdate.day
                if day<10:
                    day = '0'+str(day)
                else:
                    day = str(day)
                commands = []
                #-command for extracting band and converting to GTiff format
                com = 'gdal_translate -ot Float32 -of GTiff -b ' + day + ' ' + self.precDBPath + 'Prec_daily_WFDEI_GPCC_cl_' +\
                    year + month + '.nc ' + self.tempdir + 'temp.tif'
                commands.append(com)
                #-command to project, resample, and extract to the correct extent
                com = 'gdalwarp -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                     self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                     self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' + \
                     self.tempdir + 'temp2.tif'
                commands.append(com)
                #-command to convert to pcraster map
                com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                commands.append(com)
                pcrstr = self.pcrExtention(daynr)
                #-command to convert to pcraster map with extension and convert from mm/s to mm/d
                com = self.pcrasterModelFile('"' + self.outdir + 'prec' + pcrstr + '" = if("' + self.tempdir + 'clone.map", 3600 * 24 * "' + self.tempdir + 'temp.map")')
                commands.append('pcrcalc -f ' + com)
                self.subProcessing(commands)
                #-Progress bar
                self.counter += 1
                self.progBar.setValue(int(int(self.counter/self.procSteps*100)))
                #-Text log
                self.textLog.append('Prec ' + year + '-' + month + '-' + day)
                #-Remove unnecessary files
                self.removeFiles(self.tempdir, self.outdir)
                
            self.textLog.append('\nProcessing precipitation from ' + self.dbSource + ' database finished!')
        
        #-Else if the database is from FEWS_RFE2.0 (For South East Afrika Database) or from ERA-INTERIM (Used for Iberian Peninsula)
        elif self.dbSource == 'FEWS_RFE2.0_GSOD' or self.dbSource == 'ERA-INTERIM':
            self.textLog.append('Processing precipitation from ' + self.dbSource + ' database...\n')
            for i in range(0, self.timeSteps):
                daynr = i+1 # required for pcraster extension
                curdate = self.startDate + datetime.timedelta(days=i)
                year = str(curdate.year)
                month = curdate.month
                if month<10:
                    month = '0'+str(month)
                else:
                    month = str(month)
                day = curdate.day
                if day<10:
                    day = '0'+str(day)
                else:
                    day = str(day)
                commands = []
                #-command to project, resample, and extract to the correct extent
                com = 'gdalwarp -ot Float32 -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                     self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                     self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.precDBPath + year + month + \
                     day + '_prec.tif ' + self.tempdir + 'temp.tif'
                commands.append(com)
                #-command to convert to pcraster map
                com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp.tif ' + self.tempdir + 'temp.map'
                commands.append(com)
                pcrstr = self.pcrExtention(daynr)
                #-command to convert to pcraster map with extension
                com = self.pcrasterModelFile('"' + self.outdir + 'prec' + pcrstr + '" = if("' + self.tempdir + 'clone.map", "' + self.tempdir + 'temp.map")')
                commands.append('pcrcalc -f ' + com)
                self.subProcessing(commands)
                #-Progress bar
                self.counter += 1
                self.progBar.setValue(int(self.counter/self.procSteps*100))
                #-Text log
                self.textLog.append('Prec ' + year + '-' + month + '-' + day)
                #-Remove unnecessary files
                self.removeFiles(self.tempdir, self.outdir)
            self.textLog.append('\nProcessing precipitation from ' + self.dbSource + ' database finished!')
            
        elif self.dbSource == 'ERA5Land':
            self.textLog.append('Processing precipitation from ' + self.dbSource + ' database...\n')
            for i in range(0, self.timeSteps):
                daynr = i+1 # required for pcraster extension
                curdate = self.startDate + datetime.timedelta(days=i)
                year = str(curdate.year)
                month = curdate.month
                if month<10:
                    month = '0'+str(month)
                else:
                    month = str(month)
                day = curdate.day
                if day<10:
                    day = '0'+str(day)
                else:
                    day = str(day)
                commands = []
                #-command for extracting band and converting to GTiff format
                com = 'gdal_translate -ot Float32 -of GTiff -b ' + day + ' ' + self.precDBPath + 'ERA5_Land_' +\
                    year + '_' + month + '_PR.nc ' + self.tempdir + 'temp.tif'
                commands.append(com)
                #-command to project, resample, and extract to the correct extent
                com = 'gdalwarp -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                     self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                     self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' + \
                     self.tempdir + 'temp2.tif'
                commands.append(com)
                #-command to convert to pcraster map
                com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                commands.append(com)
                pcrstr = self.pcrExtention(daynr)
                #- no conversion is required here, ERA5 already in mm/day
                com = self.pcrasterModelFile('"' + self.outdir + 'prec' + pcrstr + '" = if("' + self.tempdir + 'clone.map","' + self.tempdir + 'temp.map")')
                commands.append('pcrcalc -f ' + com)
                self.subProcessing(commands)
                #-Progress bar
                self.counter += 1
                self.progBar.setValue(int(int(self.counter/self.procSteps*100)))
                #-Text log
                self.textLog.append('Prec ' + year + '-' + month + '-' + day)
                #-Remove unnecessary files
                self.removeFiles(self.tempdir, self.outdir)
                
            self.textLog.append('\nProcessing precipitation from ' + self.dbSource + ' database finished!')            
            
            
        else:
            self.textLog.append('Error: processing of precipitation from database not possible because database is not found')
        
    #-Create temperature forcing based on the database
    def createTempDB(self):
        #-If the database is from WFDEI (Watch forcing)
        if self.dbSource == 'WFDEI' or self.dbSource == 'ERA5Land': 
            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database...\n')
            commands = []
            if self.dbSource == 'WFDEI':
                #-Convert the WFDEI dem (NetCDF) to a GeoTiff
                com = 'gdal_translate -of GTiff ' + self.dbDem + ' ' + self.tempdir + 'temp.tif'
                commands.append(com)
                self.textLog.append('\nConverting WFDEI dem to GeoTiff')
                #-Resample the WFDEI dem to the clone projection, extent, and resolution
                com = 'gdalwarp -r bilinear -ot Float32 -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + \
                    ' -tr ' + self.t_res + ' ' + self.t_res + ' -te ' + self.xMin + ' ' + \
                    self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' ' + self.tempdir + \
                    'temp.tif ' + self.tempdir + 'temp2.tif'
                commands.append(com)
                self.textLog.append('\nResampling WFDEI dem to clone resolution, projection, and extent')
                #-Convert to PCRaster format
                com = 'gdal_translate -of PCRaster ' + self.tempdir + 'temp2.tif ' + \
                    self.tempdir + 'demWFDEI.map'
                commands.append(com)
                self.textLog.append('\nConverting WFDEI dem to pcraster format: demWFDEI.map')
                #-Difference between model dem and WFDEI dem
                com = self.pcrasterModelFile('"' + self.tempdir + 'demDiff.map' + '" = "' + self.modelDem + '" - "' + self.tempdir + 'demWFDEI.map"') 
                commands.append('pcrcalc -f ' + com)
                self.textLog.append('\nCalculating the difference between the model dem and WFDEI dem: demDiff.map\n')
                self.subProcessing(commands)
                #-Remove temporary files
                self.removeFiles(self.tempdir, self.outdir)
            #-Loop over the 3 temperature forcings and create daily maps
            forcings = {'Tair': self.tavgDBPath, 'Tmax': self.tmaxDBPath, 'Tmin': self.tminDBPath}
            for f in forcings:
                for i in range(0, self.timeSteps):
                    daynr = i+1 # required for pcraster extension
                    curdate = self.startDate + datetime.timedelta(days=i)
                    year = str(curdate.year)
                    month = curdate.month
                    if month<10:
                        month = '0'+str(month)
                    else:
                        month = str(month)
                    day = curdate.day
                    if day<10:
                        day = '0'+str(day)
                    else:
                        day = str(day)
                    self.textLog.append(f + ' ' + year + '-' + month + '-' + day)    
                    commands = []
                    #-command for extracting band and converting to GTiff format
                    if self.dbSource == 'WFDEI':
                        com = 'gdal_translate -ot Float32 -of GTiff -b ' + day + ' ' + forcings[f] + f + '_daily_WFDEI_cl_' +\
                            year + month + '.nc ' + self.tempdir + 'temp.tif'
                    
                    elif self.dbSource == 'ERA5Land':
                        com = 'gdal_translate -ot Float32 -of GTiff -b ' + day + ' ' + forcings[f] + f + '_Land_' +\
                            year + '_' + month + '.nc ' + self.tempdir + 'temp.tif'                    
                    
                    commands.append(com)
                    #-command to project, resample, and extract to the correct extent
                    com = 'gdalwarp -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                        self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                        self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' + \
                        self.tempdir + 'temp2.tif'
                    commands.append(com)
                    #-command to convert to pcraster map
                    com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                    commands.append(com)
                    pcrstr = self.pcrExtention(daynr)
                    
                    if self.dbSource == 'WFDEI':
                        #-command to calculate the correct temperature using a lapse rate and the demDiff.map
                        com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = ("' + self.tempdir + \
                            'temp.map" + (-0.0065 * "' + self.tempdir + 'demDiff.map")) - 273.15')
                            
                    elif self.dbSource == 'ERA5Land':
                        com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = ("' + self.tempdir + \
                            'temp.map" )')
                    
                                        
                    commands.append('pcrcalc -f ' + com)
                    self.subProcessing(commands)
                    #-Progress bar
                    self.counter += 1
                    self.progBar.setValue(int(self.counter/self.procSteps*100))
                    #-Remove temporary files
                    self.removeFiles(self.tempdir, self.outdir)
            
            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database finished!')
        
        #-Else if the database is GSOD interpolated stations (interpolated to reference elevation = 0 MASL).
        elif self.dbSource == 'FEWS_RFE2.0_GSOD':  
            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database...\n')
            #-Loop over the 3 temperature forcings and create daily maps
            forcings = {'tair': self.tavgDBPath, 'tmax': self.tmaxDBPath, 'tmin': self.tminDBPath}
            for f in forcings:
                for i in range(0, self.timeSteps):
                    daynr = i+1 # required for pcraster extension
                    curdate = self.startDate + datetime.timedelta(days=i)
                    year = str(curdate.year)
                    month = curdate.month
                    if month<10:
                        month = '0'+str(month)
                    else:
                        month = str(month)
                    day = curdate.day
                    if day<10:
                        day = '0'+str(day)
                    else:
                        day = str(day)
                    self.textLog.append(f + ' ' + year + '-' + month + '-' + day)    
                    commands = []
                    #-command to project, resample, and extract to the correct extent
                    com = 'gdalwarp -ot Float32 -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                        self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                        self.t_res + ' ' + self.t_res + ' -r bilinear ' + forcings[f] + f + '_' + \
                        year + month + day + '.tif ' + self.tempdir + 'temp.tif'
                    commands.append(com)
                    #-command to convert to pcraster map
                    com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp.tif ' + self.tempdir + 'temp.map'
                    commands.append(com)
                    pcrstr = self.pcrExtention(daynr)
                    #-command to calculate the correct temperature using a lapse rate and the dem.map
                    com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = ("' + self.tempdir + \
                        'temp.map" + (-0.0065 * "' + self.modelDem +'"))')
                    commands.append('pcrcalc -f ' + com)
                    self.subProcessing(commands)
                    #-Progress bar
                    self.counter += 1
                    self.progBar.setValue(int(self.counter/self.procSteps*100))
                    #-Remove temporary files
                    self.removeFiles(self.tempdir, self.outdir)
            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database finished!')
        
        #-If the database is from ERA-INTERIM (Used for Iberian Peninsula)
        elif self.dbSource == 'ERA-INTERIM':
            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database...\n')
            #-Loop over the 3 temperature forcings and create daily maps
            forcings = {'tavg': self.tavgDBPath, 'tmax': self.tmaxDBPath, 'tmin': self.tminDBPath}
            for f in forcings:
                for i in range(0, self.timeSteps):
                    daynr = i+1 # required for pcraster extension
                    curdate = self.startDate + datetime.timedelta(days=i)
                    year = str(curdate.year)
                    month = curdate.month
                    if month<10:
                        month = '0'+str(month)
                    else:
                        month = str(month)
                    day = curdate.day
                    if day<10:
                        day = '0'+str(day)
                    else:
                        day = str(day)
                    self.textLog.append(f + ' ' + year + '-' + month + '-' + day)    
                    commands = []
                    #-Calculate reference elevation temperature
                    com = 'gdal_calc.py -A ' + forcings[f] + year + month + day + '_' + f + '.tif -B ' + \
                        self.dbDem + ' --outfile=' + self.tempdir + 'temp.tif --calc="A+(B*0.0065)"'
                    commands.append(com)
                    #-command to project, resample, and extract to the correct extent
                    com = 'gdalwarp -s_srs ' + self.dbSrs + ' -t_srs ' + self.t_srs + ' -te ' + \
                        self.xMin + ' ' + self.yMin + ' ' + self.xMax + ' ' + self.yMax + ' -tr ' + \
                        self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' + \
                        self.tempdir + 'temp2.tif'
                    commands.append(com)
                    #-Convert reference temperature map to pcraster map
                    com = 'gdal_translate -of PCRaster -ot Float32 ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                    commands.append(com)
                    #-command to calculate the correct temperature using a lapse rate
                    pcrstr = self.pcrExtention(daynr)
                    com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = "' + self.tempdir + \
                        'temp.map" - (0.0065 * "' + self.modelDem + '")')
                    commands.append('pcrcalc -f ' + com)
                    self.subProcessing(commands)
                    #-Progress bar
                    self.counter += 1
                    self.progBar.setValue(int(self.counter/self.procSteps*100))
                    #-Remove temporary files
                    self.removeFiles(self.tempdir, self.outdir)

            self.textLog.append('\nProcessing temperature from ' + self.dbSource + ' database finished!')    
        else:
            self.textLog.append('\nError: processing of temperature from database not possible because database is not found')
        
    #-Create precipitation forcing based on user-defined stations        
    def createPrecCSV(self):
        if os.path.isfile(self.precLocFile) and os.path.isfile(self.precDataFile):
            self.textLog.append('Processing precipitation from user-defined CSV files...\n')
            stations = self.readStationsLoc(self.precLocFile)
            #-outsize in pixels
            outsizeX = int(math.floor((int(self.xMax) - int(self.xMin)) / float(self.t_res)))
            outsizeY = int(math.floor((int(self.yMax) - int(self.yMin)) / float(self.t_res)))
            with open(self.precDataFile, 'r') as csvfile:
                rain_data = csv.reader(csvfile, delimiter=',')
                next(rain_data)
                for row in rain_data:
                    date = datetime.datetime.strptime(row[0], '%d-%m-%Y')
                    date = datetime.date(date.year, date.month, date.day)
                    if date > self.startDate or date > self.endDate: #-check if period can be found in user defined data
                        self.textLog.append('\nError: Defined period to process can not be found in data csv file!')
                        return 0
                    elif date == self.startDate: #-stop if start date is found
                        daynr = 1
                        break
                self.textLog.append('Prec ' + str(date.year) + '-' + str(date.month) + '-' + str(date.day))
                #-Create output csv file
                f = open(os.path.join(self.tempdir, 'temp_prec.csv'), 'w')
                f.write('X,Y,temp_prec\n')
                #-Write first day of data to file
                for i in range(1, len(row)):
                    ##-TEST FOR missing values
                    if int(float(row[i])) != -9999:
                        f.write('%f,%f,%f\n' %(float(stations[i-1][1]), float(stations[i-1][2]), float(row[i])))
                f.close()
                #-Create a VRT file
                self.vrtCreate('temp_prec')
                commands = []
                #-Convert to a grid
                com = 'gdal_grid -a invdist:power=2.0:smoothing=0.0 -ot Float32 -txe ' + self.xMin + \
                    ' ' + self.xMax + ' -tye ' + self.yMin + ' ' + self.yMax + ' -outsize ' + str(outsizeX) + \
                    ' ' + str(outsizeY) + ' -l temp_prec ' + self.tempdir + 'temp_prec.vrt ' + self.tempdir + 'temp.tif'
                commands.append(com)
                #-Warp to extent and resolution
                com = 'gdalwarp -te ' + self.xMin + ' ' + self.yMin + ' '  + self.xMax + ' ' + self.yMax +\
                 ' -tr ' + self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' +\
                 self.tempdir + 'temp2.tif'
                commands.append(com)
                #-Convert to pcraster map
                com = 'gdal_translate -of PCRaster ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                commands.append(com)
                #-Clip from clone
                pcrstr = self.pcrExtention(daynr)
                com = self.pcrasterModelFile('"' + self.outdir + 'prec' + pcrstr + '" = if("' + self.tempdir + 'clone.map", "' + self.tempdir + 'temp.map")')
                commands.append('pcrcalc -f ' + com)
                self.subProcessing(commands)
                #-Progress bar
                self.counter += 1
                self.progBar.setValue(int(self.counter/self.procSteps*100))
                #-Remove temporary files
                self.removeFiles(self.tempdir, self.outdir)
                #-Continue looping until the end date
                for row in rain_data:
                    daynr+=1
                    date = datetime.datetime.strptime(row[0],'%d-%m-%Y')
                    date = datetime.date(date.year, date.month, date.day)
                    if date>self.endDate:
                        self.textLog.append('\nError: Defined period to process can not be found in data csv file!')
                        return 0
                    self.textLog.append('Prec ' + str(date.year) + '-' + str(date.month) + '-' + str(date.day))
                    #-Create output csv file
                    f = open(os.path.join(self.tempdir, 'temp_prec.csv'), 'w')
                    f.write('X,Y,temp_prec\n')
                    for i in range(1, len(row)):
                        ##-TEST FOR missing values
                        if int(float(row[i])) != -9999:
                            f.write('%f,%f,%f\n' %(float(stations[i-1][1]), float(stations[i-1][2]), float(row[i])))
                    f.close()
                    #-Create a VRT File
                    self.vrtCreate('temp_prec')
                    commands = []
                    #-Convert to a grid
                    com = 'gdal_grid -a invdist:power=2.0:smoothing=0.0 -ot Float32 -txe ' + self.xMin + \
                        ' ' + self.xMax + ' -tye ' + self.yMin + ' ' + self.yMax + ' -outsize ' + str(outsizeX) + \
                        ' ' + str(outsizeY) + ' -l temp_prec ' + self.tempdir + 'temp_prec.vrt ' + self.tempdir + 'temp.tif'
                    commands.append(com)
                    #-Warp to extent and resolution
                    com = 'gdalwarp -te ' + self.xMin + ' ' + self.yMin + ' '  + self.xMax + ' ' + self.yMax +\
                     ' -tr ' + self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' +\
                     self.tempdir + 'temp2.tif'
                    commands.append(com)
                    #-Convert to pcraster map
                    com = 'gdal_translate -of PCRaster ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                    commands.append(com)
                    #-Clip from clone
                    pcrstr = self.pcrExtention(daynr)
                    com = self.pcrasterModelFile('"' + self.outdir + 'prec' + pcrstr + '" = if("' + self.tempdir + 'clone.map", "' + self.tempdir + 'temp.map")')
                    commands.append('pcrcalc -f ' + com)
                    self.subProcessing(commands)
                    #-Progress bar
                    self.counter += 1
                    self.progBar.setValue(int(self.counter/self.procSteps*100))
                    #-Remove temporary files
                    self.removeFiles(self.tempdir, self.outdir)
                    #-Stop if end date is found
                    if date == self.endDate:
                        self.textLog.append('\nProcessing precipitation from user-defined CSV files finished!')
                        break
                #-if the CSV file has less records than the user defined end date
                if date<self.endDate:
                    self.textLog.append('\nProcessing precipitation from user-defined CSV files finished, but not all dates are processed because the data CSV file contains a shorter period than defined by the user!')
        else:
            self.textLog.append('\nError: processing of precipitation not possible because CSV files not found')

    #-Create temperature forcing based on user-defined stations        
    def createTempCSV(self):
        if os.path.isfile(self.tempLocFile) and os.path.isfile(self.tempDataFile):
            self.textLog.append('\nProcessing temperature from user-defined CSV files...\n')
            stations = self.readStationsLoc(self.tempLocFile)
            #-outsize in pixels
            outsizeX = int(math.floor((int(self.xMax) - int(self.xMin)) / float(self.t_res)))
            outsizeY = int(math.floor((int(self.yMax) - int(self.yMin)) / float(self.t_res)))
            with open(self.tempDataFile, 'r') as csvfile:
                temp_data = csv.reader(csvfile, delimiter=',')
                next(temp_data)
                for row in temp_data:
                    date = datetime.datetime.strptime(row[0], '%d-%m-%Y')
                    date = datetime.date(date.year, date.month, date.day)
                    if date > self.startDate or date > self.endDate: #-check if period can be found in user defined data
                        self.textLog.append('\nError: Defined period to process can not be found in data csv file!')
                        return 0
                    elif date == self.startDate: #-stop if start date is found
                        daynr = 1
                        break
                #-forcing CSV should be in order Tair, Tmax, Tmin
                forcings = {'Tmin': 'temp_tmin', 'Tmax': 'temp_tmax', 'Tair': 'temp_tair'}
                s = 0  # counter required for knowing which column to start in the data file
                for f in forcings:
                    self.textLog.append(f + ' ' + str(date.year) + '-' + str(date.month) + '-' + str(date.day))
                    #-Create output csv file for temperature
                    fvar = forcings[f]
                    fi = open(os.path.join(self.tempdir,  fvar + '.csv'), 'w')
                    fi.write('X,Y,' + fvar + '\n')
                    #-Write first day of data to file
                    for i in range(1, len(stations) + 1):
                        ##-TEST FOR missing values
                        if int(float(row[i+s])) != -9999:
                            Tref = float(stations[i-1][3]) * 0.0065 + float(row[i+s])  #-convert temperature to reference level temperature
                            fi.write('%f,%f,%f\n' %(float(stations[i-1][1]), float(stations[i-1][2]), Tref))
                    fi.close()
                    #-Create a VRT file
                    self.vrtCreate(fvar)
                    commands = []
                    #-Convert to a grid
                    com = 'gdal_grid -a invdist:power=2.0:smoothing=0.0 -ot Float32 -txe ' + self.xMin + \
                        ' ' + self.xMax + ' -tye ' + self.yMin + ' ' + self.yMax + ' -outsize ' + str(outsizeX) + \
                        ' ' + str(outsizeY) + ' -l ' + fvar + ' ' + self.tempdir + fvar + '.vrt ' + self.tempdir + 'temp.tif'
                    commands.append(com)
                    #-Warp to extent and resolution
                    com = 'gdalwarp -te ' + self.xMin + ' ' + self.yMin + ' '  + self.xMax + ' ' + self.yMax +\
                     ' -tr ' + self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' +\
                     self.tempdir + 'temp2.tif'
                    commands.append(com)
                    #-Convert to pcraster map
                    com = 'gdal_translate -of PCRaster ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                    commands.append(com)
                    #-Calculate the temperature using the DEM and the reference temperature
                    pcrstr = self.pcrExtention(daynr)
                    com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = ("' + self.tempdir + 'dem.map" * -0.0065) + "' + self.tempdir + 'temp.map"')
                    commands.append('pcrcalc -f ' + com)
                    self.subProcessing(commands)
                    #-Progress bar
                    self.counter += 1
                    self.progBar.setValue(int(self.counter/self.procSteps*100))
                    #-Remove temporary files
                    self.removeFiles(self.tempdir, self.outdir)
                    #-update column number
                    s = s + len(stations)
                #-Continue looping until the end date
                for row in temp_data:
                    daynr+=1
                    date = datetime.datetime.strptime(row[0],'%d-%m-%Y')
                    date = datetime.date(date.year, date.month, date.day)
                    if date>self.endDate:
                        self.textLog.append('\nError: Defined period to process can not be found in data csv file!')
                        return 0
                    s = 0  # counter required for knowing which column to start in the data file
                    for f in forcings:
                        self.textLog.append(f + ' ' + str(date.year) + '-' + str(date.month) + '-' + str(date.day))
                        #-Create output csv file for temperature
                        fvar = forcings[f]
                        fi = open(os.path.join(self.tempdir,  fvar + '.csv'), 'w')
                        fi.write('X,Y,' + fvar + '\n')
                        for i in range(1, len(stations) + 1):
                            ##-TEST FOR missing values
                            if int(float(row[i+s])) != -9999:
                                Tref = float(stations[i-1][3]) * 0.0065 + float(row[i+s])  #-convert temperature to reference level temperature
                                fi.write('%f,%f,%f\n' %(float(stations[i-1][1]), float(stations[i-1][2]), Tref))
                        fi.close()
                        #-Create a VRT file
                        self.vrtCreate(fvar)
                        commands = []
                        #-Convert to a grid
                        com = 'gdal_grid -a invdist:power=2.0:smoothing=0.0 -ot Float32 -txe ' + self.xMin + \
                            ' ' + self.xMax + ' -tye ' + self.yMin + ' ' + self.yMax + ' -outsize ' + str(outsizeX) + \
                            ' ' + str(outsizeY) + ' -l ' + fvar + ' ' + self.tempdir + fvar + '.vrt ' + self.tempdir + 'temp.tif'
                        commands.append(com)
                        #-Warp to extent and resolution
                        com = 'gdalwarp -te ' + self.xMin + ' ' + self.yMin + ' '  + self.xMax + ' ' + self.yMax +\
                         ' -tr ' + self.t_res + ' ' + self.t_res + ' -r bilinear ' + self.tempdir + 'temp.tif ' +\
                         self.tempdir + 'temp2.tif'
                        commands.append(com)
                        #-Convert to pcraster map
                        com = 'gdal_translate -of PCRaster ' + self.tempdir + 'temp2.tif ' + self.tempdir + 'temp.map'
                        commands.append(com)
                        #-Calculate the temperature using the DEM and the reference temperature
                        pcrstr = self.pcrExtention(daynr)
                        com = self.pcrasterModelFile('"' + self.outdir + f + pcrstr + '" = ("' + self.tempdir + 'dem.map" * -0.0065) + "' + self.tempdir + 'temp.map"')
                        commands.append('pcrcalc -f ' + com)
                        self.subProcessing(commands)
                        #-Progress bar
                        self.counter += 1
                        self.progBar.setValue(int(self.counter/self.procSteps*100))
                        #-Remove temporary files
                        self.removeFiles(self.tempdir, self.outdir)
                        #-update column number
                        s = s + len(stations)
                    #-Stop if end date is found
                    if date == self.endDate:
                        self.textLog.append('\nProcessing temperature from user-defined CSV files finished!')
                        break
                #-if the CSV file has less records than the user defined end date
                if date<self.endDate:
                    self.textLog.append('\nProcessing temperature from user-defined CSV files finished, but not all dates are processed because the data CSV file contains a shorter period than defined by the user!')
        else:
            self.textLog.append('\nError: processing of temperature not possible because CSV files not found')   
        
    #-Function to determine pcraster extentsion number
    def pcrExtention(self, day):
        #-Make pcraster string
        if day < 10:
            pcrstr = "0000.00"+str(day)
        elif day < 100:
            pcrstr = "0000.0"+str(day)
        elif day < 1000:
            pcrstr = "0000."+str(day)
        elif day < 10000:
            nr = str(day)
            thous = nr[0]
            hund = nr[1:4]
            pcrstr = "000"+thous+"."+hund
        else:
            nr = str(day)
            thous = nr[0:2]
            hund = nr[2:5]
            pcrstr = "00"+thous+"."+hund
        return pcrstr

    #-Remove temporary files    
    def removeFiles(self, d1, d2):
        f = glob.glob(d1 + 'temp*')
        for fi in f:
            os.remove(fi)
        f = glob.glob(d2 + '*xml')
        for fi in f:
            os.remove(fi)
            
    #-Function that creates a pcraster model file (*.mod)
    def pcrasterModelFile(self, command):
        #-Create a batch file
        batchfile = os.path.dirname(__file__) + '/pcraster/run.mod'
        f = open(batchfile, "w")
        f.write(command)
        f.close()
        return batchfile
    
    #-Subprocessing the commands to execute
    def subProcessing(self, commands):    
        for command in commands:
            if command.find('pcrcalc')>=0:
                environ = None
            else:
                environ = None
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                #universal_newlines=False,
                env=environ,
                text= True
                )
            proc = process.communicate()[0]
            if proc.find('Warning')>=0:
                # fix_print_with_import
                print('...')
            elif proc.find('WARNING')>=0:
                # fix_print_with_import
                print('...')
            else:
                # fix_print_with_import
                print(proc)
                
    #-Function to transform coordinates of stations from lat/lon to user-defined CRS
    def coordinateTransform(self, s_srs, t_srs, lon, lat):
        #-Source coordinate system
        source = osr.SpatialReference()
        source.ImportFromEPSG(s_srs)
        #-Target coordinate system
        target = osr.SpatialReference()
        target.ImportFromEPSG(t_srs)
        #-Transfrom
        transform = osr.CoordinateTransformation(source, target)
        (X, Y, Height) = transform.TransformPoint(lon, lat)
        return X, Y
    
    #-Function to read the ID, location, and elevation for all stations
    def readStationsLoc(self, f):
        stations = []
        t_srs = int((self.t_srs).split('EPSG:')[1])
        with open(f, 'r') as csvfile:
            locations = csv.reader(csvfile, delimiter=',')
            next(locations)
            for row in locations:
                locXY = self.coordinateTransform(4326, t_srs, float(row[3]), float(row[2]))
                X = locXY[0]
                Y = locXY[1]
                #statdata = [int(row[0]), X, Y, float(row[4])]
                statdata = [str(row[0]), X, Y, float(row[4])]
                stations.append(statdata)
        return stations
                
    #-Function to create VRT file
    def vrtCreate(self, fvar):
        epsg = int((self.t_srs).split('EPSG:')[1])
        f = open(self.tempdir + fvar + '.vrt', 'w')
        f.write('<OGRVRTDataSource>\n')
        f.write('<OGRVRTLayer name="' + fvar + '">\n')
        f.write('<SrcDataSource>' + self.tempdir + fvar + '.csv</SrcDataSource>\n')
        f.write('<GeometryType>wkbPoint</GeometryType>\n')
        f.write('<LayerSRS>epsg:'+str(epsg)+'</LayerSRS>\n')
        f.write('<GeometryField encoding="PointFromColumns" x="X" y="Y" z="' + fvar + '"/>\n')
        f.write('</OGRVRTLayer>\n')
        f.write('</OGRVRTDataSource>\n')
        f.close()