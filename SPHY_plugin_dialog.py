# -*- coding: utf-8 -*-
# The SPHY model Pre-Processor interface plugin for QGIS:
# A QGIS plugin that allows the user to create SPHY model input data based on a database. 
#
# Copyright (C) 2015  Wilco Terink
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
#
# Email: w.terink@futurewater.nl OR terinkw@gmail.com

#-Authorship information-###################################################################
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
__author__ = "Wilco Terink"
__copyright__ = "Wilco Terink"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "w.terink@futurewater.nl, terinkw@gmail.com"
__date__ ='1 January 2017'
############################################################################################

"""
/***************************************************************************
 SphyPreProcessDialog
                                 A QGIS plugin
 A tool to convert raw data into SPHY model input data
                             -------------------
        begin                : 2015-06-23
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Wilco Terink
        email                : w.terink@futurewater.nl
 ***************************************************************************/

"""

import os, subprocess, configparser, sqlite3, datetime, math, glob, time, processing

from qgis.PyQt import QtGui, QtCore, QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal
from qgis.core import *
from qgis.utils import iface, plugins
from qgis.gui import QgsMessageBar, QgsMapToolEmitPoint, QgsRubberBand 
from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingParameterRasterLayer
from qgis.core import QgsProcessingParameterVectorLayer
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProcessingParameterCrs
from qgis.core import QgsProcessingParameterFile
from qgis.core import QgsProcessingParameterFeatureSink
from qgis.core import QgsCoordinateReferenceSystem
from qgis.core import QgsExpression
import processing


from SphyPlugin.gui.generated.SPHY_plugin_dialog_base import Ui_SphyPluginDialog

#-Import spatial processing class with gdal commands
from SphyPlugin.aux_scripts.spatial_processing import SpatialProcessing
#-Import forcing processing 
from SphyPlugin.aux_scripts.forcing import processForcing
from SphyPlugin.aux_scripts.glaciers import Glaciers_model

#-Class that allows to drag a rectangle on the map canvas
class RectangleMapTool(QgsMapToolEmitPoint):

    # This signal is emitted when the rectangle drawing is completed, passing the drawn rectangle as an argument.
    deactivated = pyqtSignal()
    finished = QtCore.pyqtSignal(object)

    def __init__(self, canvas):
        self.canvas = canvas
        QgsMapToolEmitPoint.__init__(self, self.canvas)
        self.rubberBand = QgsRubberBand(self.canvas, QgsWkbTypes.PolygonGeometry) # A class for drawing transient features
        clr = QtGui.QColor('red')
        clr.setAlpha(50)
        self.rubberBand.setFillColor(clr)
        self.rubberBand.setWidth(1)
        self.reset()
    
    def reset(self):
        self.startPoint = self.endPoint = None
        self.isEmittingPoint = False
        self.rubberBand.reset(QgsWkbTypes.PolygonGeometry)
    
    # the event when the mouse button is pressed
    def canvasPressEvent(self, e):
        self.startPoint = self.toMapCoordinates(e.pos())
        self.endPoint = self.startPoint
        self.isEmittingPoint = True
        self.showRect(self.startPoint, self.endPoint)
    
    # handles the event when the mouse button is released
    def canvasReleaseEvent(self, e):
        self.isEmittingPoint = False
        r = self.rectangle()
        if r is not None:
            self.rubberBand.reset()
            self.finished.emit(r)
    
    # It updates the end point and redraws the rectangle.
    def canvasMoveEvent(self, e):
        if not self.isEmittingPoint:
            return
        self.endPoint = self.toMapCoordinates(e.pos())
        self.showRect(self.startPoint, self.endPoint)
    
    #  Draws the rectangle on the canvas by adding points to the rubberBand.
    def showRect(self, startPoint, endPoint):
        self.rubberBand.reset(QgsWkbTypes.PolygonGeometry)
        if startPoint.x() == endPoint.x() or startPoint.y() == endPoint.y():
            return
        point1 = QgsPointXY(startPoint.x(), startPoint.y())
        point2 = QgsPointXY(startPoint.x(), endPoint.y())
        point3 = QgsPointXY(endPoint.x(), endPoint.y())
        point4 = QgsPointXY(endPoint.x(), startPoint.y())

        self.rubberBand.addPoint(point1, False)
        self.rubberBand.addPoint(point2, False)
        self.rubberBand.addPoint(point3, False)
        self.rubberBand.addPoint(point4, True)    # true to update canvas
        self.rubberBand.show()
    
    # Returns the QgsRectangle object representing the drawn rectangle, or None if the rectangle is not properly defined.
    def rectangle(self):
        if self.startPoint is None or self.endPoint is None:
            return None
        elif self.startPoint.x() == self.endPoint.x() or self.startPoint.y() == self.endPoint.y():
            return None
        return QgsRectangle(self.startPoint, self.endPoint)
    
    # Deactivates the tool and emits a deactivated signal.
    def deactivate(self):
        super(RectangleMapTool, self).deactivate()
        self.deactivated.emit()

#-Preprocessor
class SphyPluginDialog(QtWidgets.QDialog, Ui_SphyPluginDialog):
    def __init__(self, parent=None):
        """Constructor."""
        super(SphyPluginDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        self.setupUi(self)
        #-Define the path where the plugins are installed
        self.pluginPath = os.path.dirname(__file__) + '/'
        
        #-Trigger the PYthon console in order to let the print statements work in the thread subprocessing part
        iface.actionShowPythonDialog().trigger()
        iface.actionShowPythonDialog().trigger()
        
        """
        If QGIS is loaded, check if there is a recent SPHY preprocessor config file in the registry
        if not, then create a reference to the SPHY preprocessor config template and initialize the plugin
        with that template file.
        """
        # Check if an existing config file is present from the most recent project
        self.currentConfig = configparser.ConfigParser(allow_no_value = True)
        self.settings = QtCore.QSettings()
        self.currentConfigFileName = self.settings.value("sphyPreProcessPlugin/currentConfig")
        try:
            self.currentConfig.read(self.currentConfigFileName)
            self.projectDir = os.path.dirname(self.currentConfigFileName) + '/'
            self.currentProject = True
            self.Tab.setEnabled(1)
            #self.Tab.setEnabled(self, 3, False)
            self.initGuiConfigMap()
            self.saveAsButton.setDisabled(0)
        except:
            self.currentProject = False
            self.Tab.setEnabled(0)
            self.projectDir = './'
            self.databasePath = './'
            self.resultsPath = './'
            self.saveAsButton.setDisabled(1)
            
        #self.saveButton.setDisabled(1)
        
        #-Detect and set coordinate systems
        self.mapCrs = iface.mapCanvas().mapSettings().destinationCrs()
        try:
            crs = self.lookupUTM(self.currentConfig.getint('GENERAL', 'utmZoneNr'), self.currentConfig.get('GENERAL', 'utmZoneStr'))
            self.userCRS = QgsCoordinateReferenceSystem("EPSG:" + str(crs))
        except:
            self.userCRS = self.mapCrs
        
        """ Basin delineation Tab """
        self.delineateButton.setEnabled(1)
                
        #-clear the process log text widget
        self.processLog1TextEdit.clear()
        self.processLog2TextEdit.clear()
        self.processLog3TextEdit.clear()
        self.processLog4TextEdit.clear()
        
    #-Initialize the GUI
    def initGuiConfigMap(self):
        #####-Dictionary for General settings Tab and Basin delineation Tab
        self.configDict = {'databaseLineEdit':('GENERAL', 'Database_dir'), 'resultsLineEdit':('GENERAL', 'Results_dir'),\
                           'utmSpinBox': ('GENERAL', 'utmZoneNr'),'startDateEdit': ('GENERAL', ('startyear', 'startmonth', 'startday')), 'endDateEdit': ('GENERAL', \
                           ('endyear', 'endmonth', 'endday')), 'outletsLineEdit' : ('DELINEATION', 'outlets_shp'), 'clipMaskCheckBox' : ('DELINEATION', 'clip'),\
                           'createSubBasinCheckBox' : ('DELINEATION', 'subbasins'), 'stationsLineEdit' : ('STATIONS', 'stations_shp')}
        self.setGui()
        #####-Dictionary for UTM coordinate system
        self.configRadioDict = {'utmNRadioButton': ('GENERAL', 'utmZoneStr'), 'utmSRadioButton': ('GENERAL', 'utmZoneStr')}
        self.setRadioGui()
        #####-Dictionary for Area selection Tab        
        self.configAreaDict = {'selectedAreaMapLineEdit': ('AREA', 'clone_shp'), 'spatialResolutionSpinBox': ('AREA', 'resolution'), 'numberCellsLineEdit': ('AREA', 'cells'),\
                               'areaSizeLineEdit': ('AREA', 'area'),'xminLineEdit': ('AREA', 'xmin'),'xmaxLineEdit': ('AREA', 'xmax'),'ymaxLineEdit': ('AREA', 'ymax'),\
                               'yminLineEdit': ('AREA', 'ymin'), 'columnsLineEdit': ('AREA', 'cols'), 'rowsLineEdit': ('AREA', 'rows'), 'cloneLineEdit': ('AREA', 'clone_grid')}
        self.setAreaDict()
        #####-Dictionary for Modules Tab
        self.configModulesDict = {'glacierModCheckBox': ('MODULES', 'glacier'), 'snowModCheckBox': ('MODULES', 'snow'), 'groundwaterModCheckBox': ('MODULES', 'groundwater'),\
                                  'routingModCheckBox': ('MODULES', 'routing')}
        #-general maps are always created in the "Create initial maps" Tab 
        self.generalMaps = {'DEM': 'dem.map', 'Slope': 'slope.map', 'Root_field': 'root_field.map', 'Root_sat': 'root_sat.map',\
                            'Root_dry': 'root_dry.map', 'Root_wilt': 'root_wilt.map', 'Root_Ksat': 'root_ksat.map', 'Sub_field': 'sub_field.map', 'Sub_sat': 'sub_sat.map',\
                            'Sub_Ksat': 'sub_ksat.map', 'LandUse': 'landuse.map', 'Latitudes': 'latitude.map'}
        #-glacier and routing maps are only created if these modules are turned on. Snow and groundwater modules don't require the creation of maps, but are implemented for possible
        # future developments. The Gui doesn't do anything with these two modules yet.
        self.glacierMaps = {'Glaciers Table': 'glaciers.csv'} #╔{'GlacFrac': 'glacfrac.map', 'GlacFracCI': 'glac_cleanice.map', 'GlacFracDB': 'glac_debris.map'}
        self.routingMaps = {'LDD': 'ldd.map', 'Outlets': 'outlets.map', 'Rivers': 'river.map', 'AccuFlux': 'accuflux.map', 'Sub-basins': 'subbasins.map'}
        self.setModulesDict()
        #-Dictionary for the Meteorological forcing Tab
        self.forcingDict = {'FlagCheckBox': 'FLAG', 'DBRadioButton': 'DB', 'LocFileLineEdit': 'LocFile', 'DataFileLineEdit': 'DataFile'}        
        self.setForcingDict()
    
    #-Set the GUI with the correct values (obtained from config file)
    def setGui(self):
        for key in self.configDict:
            widget = eval('self.' + key)
            i = self.configDict[key]
            module = i[0]
            pars = i[1]
            if module == 'GENERAL' and (pars[0] == 'startyear' or pars[0] == 'endyear'): 
                #self.setGui(widget, module, pars[0], pars[1], pars[2])
                widget.setDate(QtCore.QDate(self.currentConfig.getint(module, pars[0]),self.currentConfig.getint(module, pars[1]),self.currentConfig.getint(module, pars[2])))
                #self.updateDate()
            else:
                if isinstance(widget, QtWidgets.QLineEdit):
                    widget.setText(self.currentConfig.get(module, pars))
                elif isinstance(widget, QtWidgets.QCheckBox):
                    if self.currentConfig.getint(module, pars) == 1:
                        widget.setChecked(1)
                elif isinstance(widget, QtWidgets.QDoubleSpinBox):
                    widget.setValue(self.currentConfig.getfloat(module, pars))
                elif isinstance(widget, QtWidgets.QSpinBox):
                    widget.setValue(self.currentConfig.getint(module, pars))
                # define the variables self.databasePath and self.resultsPath and pcraster bin path. Set the database metadata config file as well.
                if widget == self.databaseLineEdit:
                    self.databasePath = self.currentConfig.get(module, pars)
                    if os.path.isfile(os.path.join(self.databasePath, 'metadata.cfg')):
                        self.databaseConfig = configparser.ConfigParser(allow_no_value = True)
                        self.databaseConfig.read(os.path.join(self.databasePath, 'metadata.cfg'))
                    else:
                        self.databaseConfig = False
                elif widget == self.resultsLineEdit:
                    self.resultsPath = self.currentConfig.get(module, pars)
                 #-Check if a outlet(s) file has been defined by the user
                elif widget == self.outletsLineEdit:
                    if os.path.exists(self.currentConfig.get(module, pars)):
                        self.outletsShp = self.currentConfig.get(module, pars)
                        widget.setText(self.outletsShp)
                    else:
                        self.outletsShp = False
                #-Check if a station(s) file has been defined by the user
                elif widget == self.stationsLineEdit:
                    if os.path.exists(self.currentConfig.get(module, pars)):
                        self.stationsShp = self.currentConfig.get(module, pars)
                        widget.setText(self.stationsShp)
                    else:
                        self.stationsShp = False
        #-initialize the enddate and startdate for processing forcing data
        self.enddate = datetime.date(QtCore.QDate.year(self.endDateEdit.date()),QtCore.QDate.month(self.endDateEdit.date()),QtCore.QDate.day(self.endDateEdit.date()))
        self.startdate = datetime.date(QtCore.QDate.year(self.startDateEdit.date()),QtCore.QDate.month(self.startDateEdit.date()),QtCore.QDate.day(self.startDateEdit.date()))
                
    #-Set the GUI radio buttons (obtained from config file)
    def setRadioGui(self):
        for key in self.configRadioDict:
            i = self.configRadioDict[key]
            module = i[0]
            pars = i[1]
            utmStr = self.currentConfig.get(module, pars)
            if key == 'utmNRadioButton' and utmStr == 'N':
                self.utmNRadioButton.setChecked(1)
            elif key == 'utmSRadioButton' and utmStr == 'S':
                self.utmSRadioButton.setChecked(1)
                
    #-Set the GUI for the Area selection (obtained from config file)
    def setAreaDict(self):
        for key in self.configAreaDict:
            widget = eval('self.' + key)
            i = self.configAreaDict[key]
            module = i[0]
            pars = i[1]
            if isinstance(widget, QtWidgets.QLineEdit):
                widget.setText(self.currentConfig.get(module, pars))
            elif isinstance(widget, QtWidgets.QSpinBox):
                widget.setValue(self.currentConfig.getint(module, pars))
            #-check if the user already has a clone shapefile of the area
            if widget == self.selectedAreaMapLineEdit:
                self.selectedAreaShp = self.currentConfig.get(module, pars)
                if os.path.exists(self.selectedAreaShp):
                    self.areaPropertiesGroupBox.setEnabled(1)
                else:
                    self.selectedAreaShp = False
                    self.areaPropertiesGroupBox.setEnabled(0)
            #-set the desired spatial resolution and extents
            elif widget == self.spatialResolutionSpinBox:
                self.spatialRes = self.currentConfig.getint(module, pars)
            elif widget == self.xminLineEdit:
                self.xMin = self.currentConfig.getint(module, pars)
            elif widget == self.xmaxLineEdit:
                self.xMax = self.currentConfig.getint(module, pars)
            elif widget == self.yminLineEdit:
                self.yMin = self.currentConfig.getint(module, pars)
            elif widget == self.ymaxLineEdit:
                self.yMax = self.currentConfig.getint(module, pars)
            elif widget == self.rowsLineEdit:
                self.rows = self.currentConfig.getint(module, pars)
            elif widget == self.columnsLineEdit:
                self.cols = self.currentConfig.getint(module, pars)
            #-set the clone
            elif widget == self.cloneLineEdit:
                self.clone = self.currentConfig.get(module, pars)
    
    #-Set the GUI for the Modules
    def setModulesDict(self):
        #-Clear the widget list
        self.modulesListWidget.clear()
        #-First set the modules to False (only routing and glacier are implemented because they require creation of maps)
        self.routing = False
        self.Tab.setTabEnabled(3, False)
        self.glacier = False
        #-Add the general maps to the list widget
        for k in self.generalMaps:
            self.modulesListWidget.addItem(k)
        #-check with modules to use
        for key in self.configModulesDict:
            widget = eval('self.' + key)
            i = self.configModulesDict[key]
            module = i[0]
            pars = i[1]
            widget.setChecked(self.currentConfig.getint(module, pars))
            if widget == self.routingModCheckBox and self.routingModCheckBox.checkState():
                self.routing = True
                self.Tab.setTabEnabled(3, True)
                for k in self.routingMaps:
                    self.modulesListWidget.addItem(k)
            elif widget == self.glacierModCheckBox and self.glacierModCheckBox.checkState():
                self.glacier = True
                for k in self.glacierMaps:
                    self.modulesListWidget.addItem(k)
        #-Sort the list
        self.modulesListWidget.sortItems()
        
    #-Set the GUI for the Meteorological forcing
    def setForcingDict(self):
        forcings = ['prec','temp']
        for f in forcings:
            for key in self.forcingDict:
                widget = eval('self.' + f +  key)
                pars = f + self.forcingDict[key]
                if isinstance(widget, QtWidgets.QCheckBox) or isinstance(widget, QtWidgets.QRadioButton):
                    checked = self.currentConfig.getint('FORCING', pars)
                    widget.setChecked(checked)
                    if widget == eval('self.' + f + 'FlagCheckBox'):
                        eval('self.' + f + 'GroupBox.setEnabled(checked)')
                        if checked == 1:
                            setattr(self, f + 'FLAG', True)
                        else:
                            setattr(self, f + 'FLAG', False)
                    elif widget == eval('self.' + f + 'DBRadioButton'):
                        if checked == 1:
                            setattr(self, f + 'DB', True)
                            setattr(self, f + 'CSV', False)
                            eval('self.' + f + 'FilesGroupBox.setEnabled(0)')
                        else:
                            setattr(self, f + 'DB', False)
                            setattr(self, f + 'CSV', True)
                            eval('self.' + f + 'FilesGroupBox.setEnabled(1)')
                            if os.path.isfile(self.currentConfig.get('FORCING', f + 'LocFile')) and os.path.isfile(self.currentConfig.get('FORCING', f + 'DataFile')):
                                setattr(self, f + 'LocFile', self.currentConfig.get('FORCING', f + 'LocFile'))
                                setattr(self, f + 'DataFile', self.currentConfig.get('FORCING', f + 'DataFile'))
                            else:
                                setattr(self, f + 'LocFile', False)
                                setattr(self, f + 'DataFile', False)
                            eval('self.' + f + 'CSVRadioButton.setChecked(1)')
                else:
                    widget.setText(self.currentConfig.get('FORCING', pars))
                
    #-Function that updates the paths in the GUI, and updates the config file
    def updatePath(self):
        sender =  self.sender().objectName()
        if sender == 'databaseFolderButton':
            tempname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the database folder for your area of interest', self.projectDir, QtWidgets.QFileDialog.ShowDirsOnly)
            if tempname and os.path.isfile(os.path.join(tempname, 'metadata.cfg')):
                self.databaseLineEdit.setText(tempname.replace('\\', '/') + '/')
                self.updateConfig('GENERAL', 'Database_dir', tempname.replace('\\', '/') + '/')
            else:
                iface.messageBar().pushMessage('Error:', 'No database found in the specified folder.\nSelect a different folder.' , QgsMessageBar.CRITICAL, 10)
        elif sender == 'resultsFolderButton':
            tempname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the folder where processed model input should be written', self.projectDir, QtWidgets.QFileDialog.ShowDirsOnly)
            if tempname:
                self.resultsLineEdit.setText(tempname.replace('\\', '/') + '/')
                self.updateConfig('GENERAL', 'Results_dir', tempname.replace('\\', '/') + '/')
        elif sender == 'pcrasterBinFolderButton':
            tempname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select the PCRaster bin folder', self.projectDir, QtWidgets.QFileDialog.ShowDirsOnly)
            if os.path.isfile(os.path.join(tempname, 'pcrcalc.exe')) == False:
                QtWidgets.QMessageBox.warning(self, 'PCRaster bin path error', 'PCRaster is not found in the specified folder. \nSelect a different folder.')
            else:
                self.updateConfig('GENERAL', 'Pcraster_dir', tempname.replace('\\', '/') + '/')
        self.saveProject()

    #-Update the config file        
    def updateConfig(self, module, par, value):
        self.currentConfig.set(module, par, str(value))
        self.updateSaveButtons(1)
        
    #-Update config with area properties
    def updateAreaConfig(self):
        self.updateConfig('AREA', 'clone_shp', self.selectedAreaShp)
        self.updateConfig('AREA', 'area', '%s' %int(self.area))
        self.updateConfig('AREA', 'cells', '%s' %int(self.cells))
        self.updateConfig('AREA', 'xmin', '%s' %int(self.xMin))
        self.updateConfig('AREA', 'xmax', '%s' %int(self.xMax))
        self.updateConfig('AREA', 'ymin', '%s' %int(self.yMin))
        self.updateConfig('AREA', 'ymax', '%s' %int(self.yMax))
        self.updateConfig('AREA', 'cols', '%s' %int(self.cols))
        self.updateConfig('AREA', 'rows', '%s' %int(self.rows))
        self.updateConfig('AREA', 'resolution', self.spatialRes)
                    
       
    #-Function that shows a UTM jpg image to help the user identify which UTM zone to use    
    def showUTM(self):
        subprocess.call(self.pluginPath + 'images/utm_zones.jpg', shell=True)
    
    #-Add or remove background layer (from NaturalEarth)    
    def showBackground(self, state):
        if state==QtCore.Qt.Unchecked:
            #-try to remove the background maps if they already exist in the canvas
            try:
                canvas = iface.mapCanvas()
                allLayers = canvas.layers()
                for i in allLayers:
                    if i.name() == 'OpenTopoMap' or i.name() == 'countries':
                        QgsProject.instance().removeMapLayer(i.id())
            except:
                passs
        #-Show the background layers
        else:
            #raster = iface.addRasterLayer(self.pluginPath + 'NaturalEarthData/HYP_50M_SR_W/HYP_50M_SR_W.tif', 'shaded relief')

            # URL template for OpenTopoMap
            url = "https://tile.opentopomap.org/{z}/{x}/{y}.png"

            # Create the raster layer
            layer_name = 'OpenTopoMap'
            raster = QgsRasterLayer(f"type=xyz&url={url}", layer_name, "wms")

            # Add the raster layer to the QGIS interface
            if raster.isValid():
                raster.setCrs(QgsCoordinateReferenceSystem(3857, QgsCoordinateReferenceSystem.EpsgCrsId))
                QgsProject.instance().addMapLayer(raster)
                raster.renderer().setOpacity(0.5)
                print(f"{layer_name} layer added successfully")
            else:
                print(f"Failed to add {layer_name} layer")

            countries = iface.addVectorLayer(self.pluginPath + 'images/ne_10m_admin_0_countries.shp', 'countries', 'ogr')
            countries.setCrs(QgsCoordinateReferenceSystem(4326, QgsCoordinateReferenceSystem.EpsgCrsId))
            fill_style = QgsSimpleFillSymbolLayer()
            fill_style.setBrushStyle(QtCore.Qt.NoBrush)
            countries.renderer().symbol().changeSymbolLayer(0, fill_style)
            iface.mapCanvas().refresh()
            iface.layerTreeView().refreshLayerSymbology(countries.id())

    #-Function to select a rectangle on the map for the area of interest.
    def selectArea(self):
        self.hide()
        self.deleteSelectedArea()
        iface.messageBar().pushMessage('Action:', 'Drag a rectangle for your area of interest.', 0)
        canvas = iface.mapCanvas()        
        self.selectRectangle = RectangleMapTool(canvas)
        canvas.setMapTool(self.selectRectangle)
        self.selectRectangle.finished.connect(self.areaSelectionFinished)
    
    #-Function that is called as soon as a rectangle is drawn from the selectArea function
    def areaSelectionFinished(self, rectangle):
        canvas = iface.mapCanvas()
        canvas.unsetMapTool(self.selectRectangle)
        iface.messageBar().clearWidgets()
        self.show()
        #-project rectangle to user defined utm to get the correct coordinates
        self.mapCrs = iface.mapCanvas().mapSettings().destinationCrs()
        rectangle = self.coordinateTransform(self.mapCrs, self.userCRS, rectangle)
        #-if rectangle is valid, then create clone polygon shapefile
        if rectangle:
            #-calculate correct extent properties
            rectangle = self.calculateExtent(rectangle)
            #-Transform the coordinates of the correct rectangle extent to the mapCrs for drawing
            rectangleMapCrs = self.coordinateTransform(self.userCRS, self.mapCrs, rectangle)
            #-Create polygone shapefile of the selected area
            self.createPolygonArea(rectangleMapCrs)
            #-Update the config for the area settings 
            self.updateAreaConfig()
        else:
            iface.messageBar().pushMessage('Error:', 'UTM zone is not valid for selected area. Select a different area.', QgsMessageBar.CRITICAL, 10)
            #-deactivate the lineedit and groupbox for area properties
            self.updateConfig('AREA', 'clone_shp', '')
        self.saveProject()
    
    #-Create shapefile of selected area
    def createPolygonArea(self, rectangle):
        #-Create a memory vector layer
        vl = QgsVectorLayer('Polygon', 'Selected area', 'memory')
        vl.startEditing()
        fet = QgsFeature()
        fet.setGeometry(QgsGeometry.fromRect(rectangle))
        vl.addFeature(fet)
        iface.messageBar().popWidget() #   -> uncomment to remove epsg crs warning
        #-Write the polygone to a vector shapefile
        QgsVectorFileWriter.writeAsVectorFormat(vl, self.selectedAreaShp, 'Selected area', self.mapCrs, 'ESRI Shapefile')
        #-Add the created vector shapefile to the canvas
        self.addSelectedArea()

    #-Delete the selected area shapefile if it already exists (from canvas and from disk)       
    def deleteSelectedArea(self):
        if self.selectedAreaShp:
            for l in iface.mapCanvas().layers():
                #-Remove the layer if it already exists
                if l.name() == 'Selected area':
                    QgsProject.instance().removeMapLayer(l.id())
                    #QgsVectorFileWriter.deleteShapeFile(l.source())
        else:
            #-Define the shapefile
            self.selectedAreaShp = self.resultsPath + 'area.shp'
            
    #-Add selected polygon area to canvas
    def addSelectedArea(self):
        fill_style = QgsSimpleFillSymbolLayer()
        clr = QtGui.QColor('blue')
        clr.setAlpha(50)
        fill_style.setColor(clr)
        wb = QgsVectorLayer(self.selectedAreaShp, 'Selected area', 'ogr')
        #wb.setCrs(self.mapCrs)
        wb.renderer().symbol().changeSymbolLayer(0, fill_style)
        QgsProject.instance().addMapLayer(wb, False)
        root = QgsProject.instance().layerTreeRoot()
        root.insertLayer(0, wb)
        iface.mapCanvas().refresh()
        iface.layerTreeView().refreshLayerSymbology(wb.id())
        
    #-Function to calculate selected area properties
    def calculateExtent(self, rectangle):
        self.xMin = math.floor(rectangle.xMinimum())
        self.xMax = math.ceil(rectangle.xMaximum())
        self.yMin = math.floor(rectangle.yMinimum())
        self.yMax = math.ceil(rectangle.yMaximum())
        self.cols = math.ceil((self.xMax - self.xMin) / self.spatialRes)
        self.rows = math.ceil((self.yMax - self.yMin) / self.spatialRes)
        self.xMax = self.xMin + self.cols * self.spatialRes
        self.yMax = self.yMin + self.rows * self.spatialRes
        self.cells = self.cols * self.rows
        #-give warning if cells are too much
        if self.cells > 500000 and self.cells <= 1000000:
            iface.messageBar().pushMessage('Warning:', 'Your model has > 500,000 cells. This means that model run-time will likely be > 1 hour for a 10-year simulation period. Choose a larger spatial resolution to reduce model run-time.', 10)
        if self.cells > 1000000 and self.cells <=2000000:
            iface.messageBar().pushMessage('Warning:', 'Your model has > 1,000,000 cells. This means that model run-time will likely be > 2 hours for a 10-year simulation period. Choose a larger spatial resolution to reduce model run-time.', 10)
        elif self.cells > 2000000:
            iface.messageBar().pushMessage('Warning:', 'Your model has > 2,000,000 cells. This means that model run-time will be very long!! Choose a larger spatial resolution to reduce model run-time.', 10)
        self.area = (self.cells * self.spatialRes**2) / 1000000  # to km2
        rectangle.setXMinimum(self.xMin)
        rectangle.setXMaximum(self.xMax)
        rectangle.setYMinimum(self.yMin)
        rectangle.setYMaximum(self.yMax)
        return rectangle
    
    #-Transform point from one coordinate system to the other    
    def coordinateTransform(self, s_srs, t_srs, rectangle):
        xform = QgsCoordinateTransform(s_srs, t_srs,QgsProject.instance())
        try:
            rectangle = xform.transform(rectangle)
        except:
            rectangle = None
        #-This function can result in strange coordinates is the rectangle is dragged for a different utm area than is specified
        return rectangle
    
    #-Re-create the selected area based on the user defined resolution
    def recreateArea(self):
        #-Get the spatial resolution from the GUI
        self.spatialRes = self.spatialResolutionSpinBox.value()
        #-Define the rectangle
        rectangle = QgsRectangle(self.xMin, self.yMin, self.xMax, self.yMax)
        #-calculate correct extent properties
        rectangle = self.calculateExtent(rectangle)
        #-Transform the coordinates of the correct rectangle extent to the mapCrs for drawing
        rectangleMapCrs = self.coordinateTransform(self.userCRS, self.mapCrs, rectangle)
        #-delete the current area
        self.deleteSelectedArea()
        #-Create polygone shapefile of the selected area
        self.createPolygonArea(rectangleMapCrs)
        #-Update the config file
        self.updateAreaConfig()
        self.saveProject()

    #-Function that creates a clone from the given extent
    def createClone(self):
        if os.path.isfile(self.resultsPath + 'clone.map'):
            #-check if file exists in map canvas and if so, then remove
            for l in iface.mapCanvas().layers():
                if l.source() == self.clone:
                    QgsProject.instance().removeMapLayer(l.id()) 
                    #QgsVectorFileWriter.deleteShapeFile(l.source())
            os.remove(self.resultsPath + 'clone.map')
        #-Command for creating new clone
        command = ('mapattr -s -R ' + str(int(self.rows)) + ' -C ' + str(int(self.cols)) + ' -P yb2t -B -x '\
                         + str(int(self.xMin)) + ' -y ' + str(int(self.yMax)) + ' -l ' + str(self.spatialRes)\
                         + ' ' + self.resultsPath + 'clone.map')
        #command = self.pcrasterBatchFile(command)
        subprocess.Popen(command, shell=True).wait()
        #-Check if clone was succesfully created
        if os.path.isfile(self.resultsPath + 'clone.map'):
            iface.messageBar().pushSuccess('Info:', 'Clone map was successfully created.')
            self.updateConfig('AREA', 'clone_grid', self.resultsPath + 'clone.map')
        else:
            iface.messageBar().pushCritical('Error:', 'Clone map was NOT created!')
            self.updateConfig('AREA', 'clone_grid', '')
        self.saveProject()
        
    #-Function to update the modules
    def updateModules(self, state):
        sender = self.sender().objectName()
        module = self.configModulesDict[sender][0]
        par = self.configModulesDict[sender][1]
        if state == QtCore.Qt.Unchecked:
            self.updateConfig(module, par, 0)
        else:
            self.updateConfig(module, par, 1)
        self.saveProject()
        
    #-Function that creates the initial maps based on the selected modules
    def createInitMaps(self):
        #-clear the process log text widget
        self.processLog1TextEdit.clear()   
        # maps to process to define progress in progressbar
        maps = 13
        mm = 0.
        # recreate clone -> is necessary if initial maps are created again after delineation has been performed
        # because delineation results in a smaller clone -> clipped to basin
        self.createClone()
        mm+=1
        self.initialMapsProgressBar.setValue(int(int(mm/maps*100)))
        if self.routing:
            maps += 5
        if self.glacier:
            maps += 1 #3
        #-set the map properties of resulting maps
        t_srs =  self.userCRS.authid()
        res = self.spatialRes
        extent = '-te ' + str(self.xMin) + ' ' + str(self.yMin) + ' ' + str(self.xMax) + ' ' + str(self.yMax)
        #-delete old raster layers from canvas and disk if exists
        for k in self.generalMaps:
            self.deleteLayer(os.path.join(self.resultsPath, self.generalMaps[k]), 'raster')
        for k in self.routingMaps:
            self.deleteLayer(os.path.join(self.resultsPath, self.routingMaps[k]), 'raster')
        for k in self.glacierMaps:
            self.deleteLayer(os.path.join(self.resultsPath, self.glacierMaps[k]), 'raster')

        #-remove temporary tiffs from results dir
        fi = glob.glob(os.path.join(self.resultsPath, 'temp*'))
        for f in fi:
            os.remove(f)

        # %% 1. CREATING DEM ------------------------------------------------------------
        #-create the commands to execute
        print('1. CREATING DEM')
        commands = []
        infile = os.path.join(self.databasePath, self.databaseConfig.get('DEM', 'file'))
        outfile = os.path.join(self.resultsPath, 'temp.tif')
        s_srs = 'EPSG:' + self.databaseConfig.get('DEM', 'EPSG')
        #-Create a class with the gdal methods
        m = SpatialProcessing(infile, outfile, s_srs, t_srs, res, extra=extent)
        #-Project, clip and resample
        commands.append(m.reproject())
        #-Convert to PCRaster map
        m.extra = '-of PCRaster'
        m.input = m.output
        m.output = os.path.join(self.resultsPath, self.generalMaps['DEM'])
        commands.append(m.rasterTranslate())
        #-Execute the command(s) 
        self.runCommands(commands)
        self.processLog1TextEdit.append('DEM is created')
        self.addCanvasLayer(os.path.join(self.resultsPath, 'dem.map'), 'DEM', 'raster')
        #-set progress bar value
        mm+=1
        self.initialMapsProgressBar.setValue(int(mm/maps*100))
        # %% 2. CREATING SLOPE ------------------------------------------------------------
        print('2. CREATING SLOPE')
        command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.generalMaps['Slope']) + '"'\
                                        + ' = slope(' + '"' + os.path.join(self.resultsPath, self.generalMaps['DEM']) + '"' + ')')
        #-Execute the command(s) 
        self.runCommands(['pcrcalc -f ' + command])
        self.processLog1TextEdit.append('SLOPE is created')
        self.addCanvasLayer(os.path.join(self.resultsPath, 'slope.map'), 'Slope', 'raster')

        fi = glob.glob(os.path.join(self.resultsPath, 'temp.*'))
        for f in fi:
            os.remove(f)
        #-set progress bar value
        mm+=1
        self.initialMapsProgressBar.setValue(int(mm/maps*100))
        # %% 3. CREATING LATITUDE ------------------------------------------------------------
        #-create the commands to execute
        print('3. CREATING LATITUDE')
        commands = []
        infile = os.path.join(self.databasePath, self.databaseConfig.get('LATITUDE', 'file'))
        outfile = os.path.join(self.resultsPath, 'temp.tif')
        s_srs = 'EPSG:' + self.databaseConfig.get('LATITUDE', 'EPSG')
        #-Create a class with the gdal methods
        m = SpatialProcessing(infile, outfile, s_srs, t_srs, res, extra=extent)
        #-Project, clip and resample
        commands.append(m.reproject())
        #-Convert to PCRaster map
        m.extra = '-of PCRaster'
        m.input = m.output
        m.output = os.path.join(self.resultsPath, self.generalMaps['Latitudes'])
        commands.append(m.rasterTranslate())
        #-Execute the command(s) 
        self.runCommands(commands)
        self.processLog1TextEdit.append('LATITUDE is created')
        self.addCanvasLayer(os.path.join(self.resultsPath, 'latitude.map'), 'Latitude', 'raster')
        fi = glob.glob(os.path.join(self.resultsPath, 'temp.*'))
        for f in fi:
            os.remove(f)
        #-set progress bar value
        mm+=1
        self.initialMapsProgressBar.setValue(int(mm/maps*100))
        # %% 4. CREATING LANDUSE ------------------------------------------------------------
        print('4. CREATING LANDUSE')
        #-create the commands to execute
        commands = []
        infile = os.path.join(self.databasePath, self.databaseConfig.get('LANDUSE', 'file'))
        outfile = os.path.join(self.resultsPath, 'temp.tif')
        s_srs = 'EPSG:' + self.databaseConfig.get('LANDUSE', 'EPSG')
        #-Create a class with the gdal methods
        m = SpatialProcessing(infile, outfile, s_srs, t_srs, res, resampling='mode', rtype='Int32', extra=extent)
        #-Project, clip and resample
        commands.append(m.reproject())
        #-Convert to PCRaster map
        m.extra = '-of PCRaster'
        m.input = m.output
        m.output = os.path.join(self.resultsPath, self.generalMaps['LandUse'])
        commands.append(m.rasterTranslate())
        #-Execute the command(s) 
        self.runCommands(commands)
        self.processLog1TextEdit.append('LANDUSE is created')
        self.addCanvasLayer(os.path.join(self.resultsPath, 'landuse.map'), 'Landuse', 'raster')

        #-remove temporary tiffs from results dir
        fi = glob.glob(os.path.join(self.resultsPath, 'temp.*'))
        for f in fi:
            os.remove(f)
        #-set progress bar value
        mm+=1
        self.initialMapsProgressBar.setValue(int(mm/maps*100))
        # %% 5. CREATING SOIL MAPS ------------------------------------------------------------
        print('5. CREATING SOIL MAPS')
        soilMapTiffs = {'Root_field': 'root_field_file', 'Root_sat': 'root_sat_file', 'Root_dry': 'root_dry_file'\
                        ,'Root_wilt': 'root_wilt_file', 'Root_Ksat': 'root_ksat_file', 'Sub_field': 'sub_field_file'\
                        ,'Sub_sat': 'sub_sat_file', 'Sub_Ksat': 'sub_ksat_file'}
        for smap in soilMapTiffs:
            #-create the commands to execute
            commands = []
            infile = os.path.join(self.databasePath, self.databaseConfig.get('SOIL', soilMapTiffs[smap]))
            outfile = os.path.join(self.resultsPath, 'temp.tif')
            s_srs = 'EPSG:' + self.databaseConfig.get('SOIL', 'EPSG')
            #-Create a class with the gdal methods
            m = SpatialProcessing(infile, outfile, s_srs, t_srs, res, extra=extent)
            #-Project, clip and resample
            commands.append(m.reproject())
            #-Convert to PCRaster map
            m.extra = '-of PCRaster'
            m.input = m.output
            m.output = os.path.join(self.resultsPath, self.generalMaps[smap])
            commands.append(m.rasterTranslate())
            #-Execute the command(s) 
            self.runCommands(commands)
            self.processLog1TextEdit.append(smap.upper() + ' is created')
            self.addCanvasLayer(os.path.join(self.resultsPath, smap + '.map'), smap, 'raster')

                #-remove temporary tiffs from results dir
            fi = glob.glob(os.path.join(self.resultsPath, 'temp.*'))
            for f in fi:
                os.remove(f)
            #-set progress bar value
            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
            

        # %% 6. CREATING ROUTING MAPS ------------------------------------------------------------
        print('6. CREATING ROUTING MAPS')
        if self.currentConfig.getint('MODULES', 'routing') == 1:

            ### LDD map #######
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"'\
                            + ' = lddcreate(' + '"' + os.path.join(self.resultsPath, self.generalMaps['DEM']) + '"' + ', 1e31, 1e31, 1e31, 1e31)')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('LDD is created')
            ### LDD REPAIRED map ######
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"'\
                            + ' = lddrepair(' + '"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"' + ')')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('LDD repaired is created')
            self.addCanvasLayer(os.path.join(self.resultsPath, self.routingMaps['LDD']), 'LDD', 'raster')

            #-set progress bar value
            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
            ###### Accuflux map #############
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['AccuFlux']) + '"'\
                            + ' = accuflux(' + '"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"' + ',1)')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('Accuflux is created')
            self.addCanvasLayer(os.path.join(self.resultsPath, self.routingMaps['AccuFlux']), 'AccuFlux', 'raster')

            #-set progress bar value
            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
            ####### Rivers map ########## assumed >50 cells
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['Rivers']) + '"'\
                            + ' = ' + '"' + os.path.join(self.resultsPath, self.routingMaps['AccuFlux']) + '"' + '> 50')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('Rivers is created')
            self.addCanvasLayer(os.path.join(self.resultsPath, self.routingMaps['Rivers']), 'Rivers', 'raster')

            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
            ######## Outlets #########  is initially the pits in the ldd (later the user can specify outlets)
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['Outlets']) + '"'\
                            + ' = pit(' + '"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"' + ')')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('Outlets is created')
            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
            ######## Sub-basins #########  is initially based on outlets (pits in the ldd (later the user can specify outlets and re-create sub-basins))
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['Sub-basins']) + '"'\
                            + ' = subcatchment(' + '"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"' + ',' + '"' + os.path.join(self.resultsPath, self.routingMaps['Outlets']) + '"' + ')')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            self.processLog1TextEdit.append('SubBasins is created')
            mm+=1
            self.initialMapsProgressBar.setValue(int(mm/maps*100))
                
        # %% 7. CREATING GLACIER MAPS ------------------------------------------------------------
        print('CREATING GLACIER MAPS')
        if self.currentConfig.getint('MODULES', 'glacier') == 1:
            print('Running glaciers model')
            processing.run("model:glaciers_model", 
                           {'clone_map': os.path.join(self.resultsPath, 'clone.map'),
                            'rgi_shapefile':os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'rgi_file')),
                            'debris_tiff': os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'debris_file')),
                            'dem':os.path.join(self.databasePath, self.databaseConfig.get('DEM', 'file')),
                            'ferrinoti_tiff':os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'ferrinoti_file')),
                            'model_resolution':self.spatialRes,'model_crs':t_srs,
                            'finer_resolution':self.spatialRes/10,'output_folder':self.resultsPath,
                            'glaciers':'TEMPORARY_OUTPUT','rgi_clipped_reproject_glac_id':'TEMPORARY_OUTPUT',
                            'intersection_glaciers_uid':'TEMPORARY_OUTPUT','ice_depth':'TEMPORARY_OUTPUT',
                            'debris':'TEMPORARY_OUTPUT','frac_glac':'TEMPORARY_OUTPUT','mod_id':'TEMPORARY_OUTPUT',
                            'modid_int_glacid':'TEMPORARY_OUTPUT','u_id':'TEMPORARY_OUTPUT','modid_int_glacid_inclmodh':'TEMPORARY_OUTPUT',
                            'intersection_glaciers_uid_hglac':'TEMPORARY_OUTPUT','debris_geom':'TEMPORARY_OUTPUT'})
            print('Glaciers Module done')

            # It does not work for now
            # create an instance of the glaciers model class
            # model = Glaciers_model()
            # parameters = {'clone_map': os.path.join(self.resultsPath, 'clone.map'),
            #                 'rgi_shapefile':os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'rgi_file')),
            #                 'debris_tiff': os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'debris_file')),
            #                 'dem':os.path.join(self.databasePath, self.databaseConfig.get('DEM', 'file')),
            #                 'ferrinoti_tiff':os.path.join(self.databasePath, self.databaseConfig.get('GLACIER', 'ferrinoti_file')),
            #                 'model_resolution':self.spatialRes,'model_crs':t_srs,
            #                 'finer_resolution':self.spatialRes/10,'output_folder':self.resultsPath,
            #                 'glaciers':'TEMPORARY_OUTPUT','rgi_clipped_reproject_glac_id':'TEMPORARY_OUTPUT',
            #                 'intersection_glaciers_uid':'TEMPORARY_OUTPUT','ice_depth':'TEMPORARY_OUTPUT',
            #                 'debris':'TEMPORARY_OUTPUT','frac_glac':'TEMPORARY_OUTPUT','mod_id':'TEMPORARY_OUTPUT',
            #                 'modid_int_glacid':'TEMPORARY_OUTPUT','u_id':'TEMPORARY_OUTPUT','modid_int_glacid_inclmodh':'TEMPORARY_OUTPUT',
            #                 'intersection_glaciers_uid_hglac':'TEMPORARY_OUTPUT','debris_geom':'TEMPORARY_OUTPUT'}
            
            # # Prepare the context and feedback
            # context = QgsProcessingContext()
            # feedback = QgsProcessingFeedback()

            # # Run the model
            # try:
            #     results = model.processAlgorithm(parameters, context, feedback)
            #     print('Model run successfully!')
            #     print('Results:', results)
            # except Exception as e:
            #     print('Error:', str(e))
                   
        self.initialMapsProgressBar.setValue(100)
        time.sleep(1)
        self.processLog1TextEdit.append('Processing is finished')
        self.initialMapsProgressBar.setValue(0)
        #-Try to delete the temporary layers from canvas and disk
        self.deleteLayer(os.path.join(self.resultsPath, 'temp.shp'), 'shape')
        self.deleteLayer(os.path.join(self.resultsPath, 'temp.tif'), 'raster')
        self.deleteLayer(os.path.join(self.resultsPath, 'temp2.tif'), 'raster')
        self.deleteLayer(os.path.join(self.resultsPath, 'temp3.tif'), 'raster')
        #-Activate the delineation button in the "Basin delineation" Tab
        self.delineateButton.setEnabled(1)

        #-remove AUX FILES
        fi = glob.glob(os.path.join(self.resultsPath, '*.aux.xml'))
        for f in fi:
            os.remove(f)
        
    # %% END OF INITAL MAPS CREATION ------------------------------------------------------------

    #-Function to update the basin delineation settings        
    def updateDelineation(self, state):
        sender = self.sender().objectName()
        if sender == 'selectOutletsButton':
            outlets = QtWidgets.QFileDialog.getOpenFileName(self, "Select the Outlet(s) shapefile", self.resultsPath, "outlets.shp")[0]
            if outlets:
                self.deleteLayer(outlets, 'shape', remLayerDisk=False)
                self.outletsShp = outlets
                self.updateConfig('DELINEATION', 'outlets_shp', outlets)
                self.addCanvasLayer(outlets, 'Outlet(s)', 'shape')
        else:
            module = self.configDict[sender][0]
            par = self.configDict[sender][1]
            if state == QtCore.Qt.Unchecked:
                self.updateConfig(module, par, 0)
            else:
                self.updateConfig(module, par, 1)
        self.saveProject()
        self.Tab.setCurrentIndex(3)
        
    #-Function to Delineate basin and create sub-basins and clipped maps
    def delineate(self):
        self.processLog2TextEdit.clear()
        #-Check if outlet(s) shapefile and ldd are present
        if self.outletsShp and os.path.isfile(os.path.join(self.resultsPath, self.routingMaps['LDD'])):
            #-set counters for progress bar
            mm = 0.
            maps = 2 # basin and outlets
            if self.currentConfig.getint('DELINEATION', 'subbasins') == 1:
                maps += 1
            if self.currentConfig.getint('DELINEATION', 'clip') == 1:
                maps += len(self.generalMaps) + len(self.routingMaps) + 2
            if self.glacier:
                maps += len(self.glacierMaps)
            #-delete old raster layers from canvas
            for k in self.routingMaps:
                    self.deleteLayer(os.path.join(self.resultsPath, self.routingMaps[k]), 'raster', remLayerDisk=False)
            #-delete basin map if exists
            self.deleteLayer(os.path.join(self.resultsPath, 'basin.map'), 'raster')
            #####-First convert outlet(s) shapefile to GeoTiff
            extent = str(self.xMin) + ',' + str(self.xMax) + ',' + str(self.yMin) + ',' + str(self.yMax)
            extent_mod = str(self.xMin) + ',' + str(self.yMin) + ',' + str(self.xMax) + ',' + str(self.yMax)
            outfile = os.path.join(self.resultsPath, 'temp.tif')
            self.processLog2TextEdit.append('Converting Outlet(s) to raster...')

            # Rasterize outlets.shp
            command = 'gdal_rasterize -l outlets ' + ' -a id -tr ' + str(self.spatialRes) + ' ' + str(self.spatialRes) + ' -te ' + extent_mod.replace(',', ' ') + ' -ot Float32 -of GTiff ' + '-a_nodata 9999.0 ' + self.outletsShp + ' ' + outfile
            #-Execute the command(s) 
            self.runCommands([command])

            #####-Translate GeoTiff to PCRaster map
            infile = outfile
            outfile = os.path.join(self.resultsPath, 'temp.map')
            #-Create a class with the gdal methods
            m = SpatialProcessing(infile, outfile, None, None, None, extra='-of PCRaster -ot Float32')
            self.runCommands([m.rasterTranslate()])

            #-Convert to nominal PCRaster map
            infile = outfile 
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['Outlets']) + '"'\
                + ' = nominal("' + infile + '")')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])

            mm += 1
            self.delineateProgressBar.setValue(int(mm/maps*100))
            self.processLog2TextEdit.append('Delineating basin...')
            #-Delineate the basin based on the ldd and the defined outlets
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, 'basin.map') + '"'\
                + ' = catchment("' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '","' + os.path.join(self.resultsPath, self.routingMaps['Outlets']) + '")')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])
            #-update progressbar
            mm += 1
            self.delineateProgressBar.setValue(int(mm/maps*100))
            #-Create sub-basins
            if self.currentConfig.getint('DELINEATION', 'subbasins') == 1:
                self.processLog2TextEdit.append('Creating sub-basins...')
                command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['Sub-basins']) + '"'\
                    + ' = subcatchment("' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '","' + os.path.join(self.resultsPath, self.routingMaps['Outlets']) + '")')
                #-Execute the command(s) 
                self.runCommands(['pcrcalc -f ' + command])
                #-update progressbar
                mm += 1
                self.delineateProgressBar.setValue(int(mm/maps*100))
                    
            #----------------------------- Clip all the maps to the basin outline  ------------------------------------------------------------
            if self.currentConfig.getint('DELINEATION', 'clip') == 1:
                self.processLog2TextEdit.append('Clipping maps to basin outline...')
                
                #-re-create the clone, based on the delineated basin map
                command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, 'clone.map') + '"'\
                    + ' = boolean("' + os.path.join(self.resultsPath, 'basin.map') +  '")')
                #-Execute the command(s) 
                self.runCommands(['pcrcalc -f ' + command])
                mm += 1
                self.delineateProgressBar.setValue(int(mm/maps*100))
                #-delete old general raster layers from canvas and clip general maps to basin outline and add to canvas
                for k in self.generalMaps:
                    self.deleteLayer(os.path.join(self.resultsPath, self.generalMaps[k]), 'raster', remLayerDisk=False)
                    command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.generalMaps[k]) + '"'\
                        + ' = if("' + os.path.join(self.resultsPath, 'clone.map') + '","' + os.path.join(self.resultsPath, self.generalMaps[k]) +  '")')
                    #-Execute the command(s) 
                    self.runCommands(['pcrcalc -f ' + command])
                    #-update progressbar
                    mm += 1
                    self.delineateProgressBar.setValue(int(mm/maps*100))
                         
                #-delete old routing raster layers from canvas and clip routing maps to basin outline and add to canvas
                for k in self.routingMaps:
                    command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps[k]) + '"'\
                        + ' = if("' + os.path.join(self.resultsPath, 'clone.map') + '","' + os.path.join(self.resultsPath, self.routingMaps[k]) +  '")')
                    #-Execute the command(s) 
                    self.runCommands(['pcrcalc -f ' + command])
                    #-update progressbar
                    mm += 1
                    self.delineateProgressBar.setValue(int(mm/maps*100))    
                
                #-repair ldd, because clipping may result in unsound ldd
                self.deleteLayer(os.path.join(self.resultsPath, self.routingMaps['LDD']), 'raster', remLayerDisk=False)
                command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"'\
                            + ' = lddrepair(' + '"' + os.path.join(self.resultsPath, self.routingMaps['LDD']) + '"' + ')')
                #-Execute the command(s) 
                self.runCommands(['pcrcalc -f ' + command])
                
                #-Clip basin map to outline
                command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, 'basin.map') + '"'\
                    + ' = if("' + os.path.join(self.resultsPath, 'clone.map') + '","' + os.path.join(self.resultsPath, 'basin.map') +  '")')
                #-Execute the command(s) 
                self.runCommands(['pcrcalc -f ' + command])
                #-update progressbar
                mm += 1
                self.delineateProgressBar.setValue(int(mm/maps*100))
                #-Disable the delineation button if clipped option is used, because otherwise errors occur in output maps if clicked again.
                #-Therefore the modules section (create initial maps) needs to be run first again in order to re-activate this button
                self.delineateButton.setEnabled(0)
            
            #-if no clipping, then add the routing maps again to the mapcanvas
            else:
                for k in self.routingMaps:
                    self.addCanvasLayer(os.path.join(self.resultsPath, self.routingMaps[k]), k, 'raster')
                self.addCanvasLayer(os.path.join(self.resultsPath, 'basin.map'), 'Basin', 'raster')
                if self.currentConfig.getint('DELINEATION', 'subbasins') == 0:
                    self.deleteLayer(os.path.join(self.resultsPath, 'subbasins.map'), 'raster')
                
            self.delineateProgressBar.setValue(100)
            time.sleep(1)
            self.processLog2TextEdit.append('Basin delineation finished.')
            self.delineateProgressBar.setValue(0)
   
        else:
            self.processLog2TextEdit.append('Error: missing outlets.shp and/or ldd.map in output folder.')
        #-remove temporary tiffs from results dir
        fi = glob.glob(os.path.join(self.resultsPath, 'temp*'))
        for f in fi:
            os.remove(f)
        #-remove AUX FILES
        fi = glob.glob(os.path.join(self.resultsPath, '*.aux.xml'))
        for f in fi:
            os.remove(f)

    #-Function to update the station settings
    def updateStations(self):
        stations = QtWidgets.QFileDialog.getOpenFileName(self, "Select the station(s) shapefile", self.resultsPath, "stations.shp")[0]
        if stations:
            self.deleteLayer(stations, 'shape', remLayerDisk=False)
            self.stationsShp = stations
            self.updateConfig('STATIONS', 'stations_shp', stations)
            self.addCanvasLayer(stations, 'Station(s)', 'shape')
        self.saveProject()
        self.Tab.setCurrentIndex(4)

    #-Function to convert stations.shp to pcraster nominal map with stations        
    def createStations(self):
        self.processLog3TextEdit.clear()
        #-Check if station(s) shapefile is present
        if self.stationsShp:
            #-delete stations map if exists
            self.deleteLayer(os.path.join(self.resultsPath, 'stations.map'), 'raster')
            #####-First convert station(s) shapefile to GeoTiff
            extent = str(self.xMin) + ',' + str(self.xMax) + ',' + str(self.yMin) + ',' + str(self.yMax)
            extent_mod = str(self.xMin) + ',' + str(self.yMin) + ',' + str(self.xMax) + ',' + str(self.yMax)

            outfile = os.path.join(self.resultsPath, 'temp.tif')
            self.processLog3TextEdit.append('Converting Station(s) to raster...')
            # Rasterize stations.shp
            command = 'gdal_rasterize -l stations ' + ' -a id -tr ' + str(self.spatialRes) + ' ' + str(self.spatialRes) + ' -te ' + extent_mod.replace(',', ' ') + ' -ot Float32 -of GTiff ' + '-a_nodata 9999.0 ' + self.stationsShp + ' ' + outfile
            #-Execute the command(s) 
            self.runCommands([command])

            #####-Translate GeoTiff to PCRaster map
            infile = outfile
            outfile = os.path.join(self.resultsPath, 'temp.map')
            #-Create a class with the gdal methods
            m = SpatialProcessing(infile, outfile, None, None, None, extra='-of PCRaster -ot Float32')
            self.runCommands([m.rasterTranslate()])

            #-Convert to nominal PCRaster map
            infile = outfile 
            command = self.pcrasterModelFile('"' + os.path.join(self.resultsPath, 'stations.map') + '"'\
                + ' = nominal("' + infile + '")')
            #-Execute the command(s) 
            self.runCommands(['pcrcalc -f ' + command])                
                
            self.processLog3TextEdit.append('Station creation finished.')
        else:
            self.processLog3TextEdit.append('Error: missing stations.shp in output folder.')
        #-remove temporary tiffs from results dir
        fi = glob.glob(os.path.join(self.resultsPath, 'temp*'))
        for f in fi:
            os.remove(f)

        #-remove AUX FILES
        fi = glob.glob(os.path.join(self.resultsPath, '*.aux.xml'))
        for f in fi:
            os.remove(f)
    
    #-Function to update the forcing settings
    def updateForcing(self, state):
        sender = self.sender()
        senderName = sender.objectName()
        var = senderName[0:4]
        key = senderName.split(var)[1]
        #-If the select csv buttons are clicked
        if isinstance(sender, QtWidgets.QToolButton): 
            toolDict = {'LocToolButton': 'LocFile', 'DataToolButton': 'DataFile'}
            pars = var + toolDict[key]
            if key[0:3] == 'Loc':
                f = QtWidgets.QFileDialog.getOpenFileName(self, "Select the location csv file", self.resultsPath, "*.csv")[0]
            else:
                f = QtWidgets.QFileDialog.getOpenFileName(self, "Select the data csv file", self.resultsPath, "*.csv")[0]
            if os.path.isfile(f):
                self.updateConfig('FORCING', pars, f)
        #-Otherwise update the checkbox or radiobutton values in the config file
        else:
            pars = var + self.forcingDict[key]
            if isinstance(sender, QtWidgets.QCheckBox):
                if state == QtCore.Qt.Checked:
                    self.updateConfig('FORCING', var + 'FLAG', 1)
                else:
                    self.updateConfig('FORCING', var + 'FLAG', 0)
            else:
                if state:
                    self.updateConfig('FORCING', var + 'DB', 1)
                else:
                    self.updateConfig('FORCING', var + 'DB', 0)
        self.saveProject()
        
    def createForcing(self):
        self.processLog4TextEdit.clear()
        #-settings for progressbar
        timeSteps = ((self.enddate-self.startdate).days + 1)
        procSteps = 0.
        if self.precFLAG:
            procSteps+=timeSteps
        if self.tempFLAG:
            procSteps+=(timeSteps * 3)
        #-Create instance of processForcing
        f = processForcing(self.resultsPath, self.userCRS.authid(), self.spatialRes, [self.xMin, self.yMin, self.xMax, self.yMax],\
            self.startdate, self.enddate, self.processLog4TextEdit, self.forcingProgressBar, procSteps)
        if self.precDB or self.tempDB:
            #-Database properties
            f.dbSource = self.databaseConfig.get('METEO', 'source')
            f.dbTs = self.databaseConfig.get('METEO', 'file_timestep')
            f.dbSrs = 'EPSG:' + self.databaseConfig.get('METEO', 'EPSG')
            f.dbFormat = self.databaseConfig.get('METEO', 'format')
        if self.precFLAG:
            if self.precDB:
                f.precDBPath = os.path.join(self.databasePath, self.databaseConfig.get('METEO', 'prec_folder'))
                f.createPrecDB() 
            else:
                f.precLocFile = self.precLocFile
                f.precDataFile = self.precDataFile 
                f.createPrecCSV()
        if self.tempFLAG:
            if self.tempDB:
                f.tavgDBPath = os.path.join(self.databasePath, self.databaseConfig.get('METEO', 'tavg_folder'))
                f.tmaxDBPath = os.path.join(self.databasePath, self.databaseConfig.get('METEO', 'tmax_folder'))
                f.tminDBPath = os.path.join(self.databasePath, self.databaseConfig.get('METEO', 'tmin_folder'))
                f.dbDem = os.path.join(self.databasePath, self.databaseConfig.get('METEO', 'dem'))
                f.modelDem = os.path.join(self.resultsPath, self.generalMaps['DEM'])
                f.createTempDB() 
            else:
                f.tempLocFile = self.tempLocFile
                f.tempDataFile = self.tempDataFile 
                f.createTempCSV()
        if not self.precFLAG and not self.tempFLAG:
            self.processLog4TextEdit.append('Nothing to process.')
        self.forcingProgressBar.setValue(100)
        time.sleep(.5)
        self.forcingProgressBar.setValue(0)
        

    def runCommands(self,commands):
        for command in commands:
            result = subprocess.run(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check= True
            )

        if result.returncode != 0:
            print([f"Command '{command}' failed with return code {result.returncode}"])
        else:
            print([f"Command '{command}' completed successfully"])
        
    #-Function that creates a pcraster model file (*.mod)
    def pcrasterModelFile(self, command):
        #-Create a batch file
        batchfile = os.path.dirname(__file__) + '/aux_scripts/pcraster/run.mod'
        f = open(batchfile, "w")
        f.write(command)
        f.close()
        return batchfile
        
    #-Change user specified CRS if spinbox value is changed or radio button is toggled
    def changeCRS(self, enabled):
        sender = self.sender().objectName()
        widget = eval('self.' + sender)
        if isinstance(widget, QtWidgets.QSpinBox):
            value = widget.value()
            self.updateConfig('GENERAL', 'utmZoneNr', value)
        else: #-then it is a radio button
            if sender == "utmNRadioButton" and enabled:
                self.updateConfig('GENERAL', 'utmZoneStr', 'N')
            else:
                self.updateConfig('GENERAL', 'utmZoneStr', 'S')
        crs = self.lookupUTM(self.currentConfig.getint('GENERAL', 'utmZoneNr'), self.currentConfig.get('GENERAL', 'utmZoneStr'))
        self.userCRS = QgsCoordinateReferenceSystem("EPSG:" + str(crs))
        self.saveProject()  #-if this causes hanging, then remove this line

    #-function to lookup the WGS84 / UTM zone from the sqlite database in QGIS
    def lookupUTM(self, utmNr, utmStr):
        con = sqlite3.connect(QgsApplication.srsDatabaseFilePath())
        cur = con.cursor()
        cur.execute("select * from tbl_srs WHERE description = 'WGS 84 / UTM zone %s%s'" %(str(utmNr), utmStr))
        epsg = cur.fetchone()
        cur.close()
        con.close()
        #-return the EPSG code
        return epsg[5]

    # validate start and enddate and set in config        
    def updateDate(self):
#         if self.exitDate: # don't execute this function if GUI is initialized during new project or open project creation.
#             return
        # validate if simulation settings are ok
        startdate = self.startDateEdit.date()
        enddate = self.endDateEdit.date()
        if startdate >= enddate:
            QtWidgets.QMessageBox.warning(self, "Date error", "End date should be larger than start date")
            if self.sender().objectName() == "startDateEdit":
                enddate = QtCore.QDate(startdate).addDays(1)
                self.endDateEdit.setDate(enddate)
            else:
                startdate = QtCore.QDate(enddate).addDays(-1)
                self.startDateEdit.setDate(startdate)
        self.updateConfig("GENERAL", "startyear", QtCore.QDate.year(startdate))
        self.updateConfig("GENERAL", "startmonth", QtCore.QDate.month(startdate))
        self.updateConfig("GENERAL", "startday", QtCore.QDate.day(startdate))
        self.updateConfig("GENERAL", "endyear", QtCore.QDate.year(enddate))
        self.updateConfig("GENERAL", "endmonth", QtCore.QDate.month(enddate))
        self.updateConfig("GENERAL", "endday", QtCore.QDate.day(enddate))
        self.saveProject()
        
    ############### DELETE AND ADD LAYERS FROM CANVAS AND DISK #############
    
    #-Add a shapefile or raster layer to the canvas (specified by name and type = vector or raster)
    def addCanvasLayer(self, filename, name, ftype):
        if ftype == 'raster':
            layer = QgsRasterLayer(filename, name)
            layer.setCrs(QgsCoordinateReferenceSystem(self.userCRS))
            QgsProject.instance().addMapLayer(layer, False)
            root = QgsProject.instance().layerTreeRoot()
            root.insertLayer(0, layer)
            self.rasterSymbology(filename, layer)
            iface.messageBar().popWidget() #   -> uncomment to remove epsg crs warning
        elif ftype == 'shape':
            layer = QgsVectorLayer(filename, name, 'ogr')
            #layer.setCrs(QgsCoordinateReferenceSystem(self.userCRS))
            QgsProject.instance().addMapLayer(layer, False)
            root = QgsProject.instance().layerTreeRoot()
            root.insertLayer(0, layer)
        
    #-Delete a shapefile or raster layer from the canvas and disk (specified by name and type = vector or raster)
    def deleteLayer(self, filename, ftype, remLayerDisk=True):
        for l in list(QgsProject.instance().mapLayers().values()):
        #for l in iface.legendInterface().layers(): old QGIS2             
            # fix_print_with_import
            print(l.name())
            if l.source() == filename:
                QgsProject.instance().removeMapLayer(l.id())
                # fix_print_with_import
                # fix_print_with_import
                print(filename + ' removed')
                break
        #-Check if file should also be removed from disk
        if os.path.isfile(filename) and remLayerDisk:
            if ftype == 'raster':
                try:
                    os.remove(filename)
                except:
                    pass
            else:
                QgsVectorFileWriter.deleteShapeFile(filename)
        
    #-Function that applies a certain color scheme to a raster layer        
    def rasterSymbology(self, filename, layer):
        #-Define the lower, mid, and high value colors
        lcolor = [239, 90, 36] # sort of red
        mcolor = [255, 255, 153] #  sort of yellow
        hcolor = [0, 102, 204] # sort of blue
        #-Get layer properties
        provider = layer.dataProvider()
        extent = layer.extent()
        stats = provider.bandStatistics(1, QgsRasterBandStats.All,extent, 0)
        #-Calculate some raster statistics
        mean_value = stats.mean
        std2_value  = stats.stdDev * 2 # 2 x std
        minvalue = max(mean_value - std2_value, 0)
        maxvalue = max(mean_value + std2_value, 0)
        diff = maxvalue - minvalue
        classes = 5
        it_step = diff/(classes-1)
        #-Number of color classes to show
        classes = 5
        #-The middle class number
        mclass = math.ceil(classes/2.0)
        #-Iterative colorsteps between the first and middle color class, and middle and high color class, respectively
        cstep1 = []
        cstep2 = []
        for c in range(0, 3):
            cstep1.append((mcolor[c]-lcolor[c])/(mclass-1))
            cstep2.append((hcolor[c]-mcolor[c])/(mclass-1))
        lst = [] #-Empty list for adding the colors and values
        for c in range(1, classes+1):
            if c==1: #- first color class = minimum value
                color  = lcolor
                rastervalue = minvalue
            elif c<mclass:
                color = [lcolor[0] + (c-1) * cstep1[0], lcolor[1] + (c-1) * cstep1[1], lcolor[2] + (c-1) * cstep1[2]]
                rastervalue = minvalue + it_step * (c-1)
            elif c==mclass: #-middle color class
                color = mcolor
                rastervalue = minvalue + it_step * (c-1)
            elif c<classes:
                color = [mcolor[0] + (c-mclass) * cstep2[0], mcolor[1] + (c-mclass) * cstep2[1], mcolor[2] + (c-mclass) * cstep2[2]]
                rastervalue = minvalue + it_step * (c-1)
            elif c==classes:
                color = hcolor
                rastervalue = maxvalue
            lst.append(QgsColorRampShader.ColorRampItem(rastervalue, QtGui.QColor(int(color[0]), int(color[1]), int(color[2])),str(rastervalue)))

        myRasterShader = QgsRasterShader()
        myColorRamp = QgsColorRampShader(minimumValue = minvalue, maximumValue = maxvalue)
        
        myColorRamp.setColorRampItemList(lst)
        myColorRamp.setColorRampType(QgsColorRampShader.Interpolated)
        myRasterShader.setRasterShaderFunction(myColorRamp)
        
        myPseudoRenderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 
                                                            layer.type(),  
                                                            myRasterShader)
        
        layer.setRenderer(myPseudoRenderer)
        layer.triggerRepaint()
        #iface.legendInterface().refreshLayerSymbology(layer) old QGIS2
        iface.layerTreeView().refreshLayerSymbology(layer.id())



    ################ NEW PROJECT, OPEN PROJECT, SAVE AND SAVE AS ############

    # function launched when new project button is clicked
    def createNewProject(self):
        # check for current project and ask to save
        if self.currentProject:
            mes = QtWidgets.QMessageBox()
            mes.setWindowTitle("Save current project")
            mes.setText("Do you want to save the current project?")
            mes.setStandardButtons(QtWidgets.QMessageBox.Save | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)
            ret = mes.exec_()
            if ret == QtWidgets.QMessageBox.Save:
                self.saveProject()
                newproject = True
            elif ret == QtWidgets.QMessageBox.No: # create new project without saving current one
                newproject = True
            else:
                newproject = False
        else:
            newproject = True
        # check if a new project needs/can be created based on the criteria tested above    
        if newproject:
            self.currentConfig.read(os.path.join(os.path.dirname(__file__), "config", "preprocess_config_template.cfg"))
            # clear project canvas
            #qgsProject = QgsProject.instance()
            #qgsProject.clear()
            # save as new project
            self.saveAsProject("new")
            self.Tab.setCurrentIndex(0)
            
            
    # function launched when existing project is opened
    def openProject(self):
        # check for the current project and ask to save
        if self.currentProject:
            mes = QtWidgets.QMessageBox()
            mes.setWindowTitle("Save current project")
            mes.setText("Do you want to save the current project?")
            mes.setStandardButtons(QtWidgets.QMessageBox.Save | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)
            ret = mes.exec_()
            if ret == QtWidgets.QMessageBox.Save:
                self.saveProject()
                openproject = True
            elif ret == QtWidgets.QMessageBox.No: # open project without saving current one
                openproject = True
            else:
                openproject = False
        else:
            openproject = True
        # check if a project can be opened based on the criteria tested above
        if openproject:
            tempname = QtWidgets.QFileDialog.getOpenFileName(self, "Open project *.cfg", self.projectDir,"*.cfg")[0]
            if tempname:
                # set the new config file
                self.currentConfigFileName = tempname
                self.currentConfig.read(self.currentConfigFileName)
                
#                 # open the corresponding qgs file
#                 qgsProjectFileName = ((self.currentConfigFileName).split(".cfg")[0]) + ".qgs"
#                 qgsProject = QgsProject.instance()
#                 qgsProject.clear()
#                 qgsProject.setFileName(qgsProjectFileName)
#                 qgsProject.read()
                # save the project
                self.saveProject()
            
    # Save as project
    def saveAsProject(self, ptype=False):
        if ptype:
            tempname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save '+ptype+' project as', self.projectDir, '*.cfg')[0]
        else:
            tempname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save current project as', self.projectDir, '*.cfg')[0]
        if tempname:
            self.currentConfigFileName = tempname
            self.saveProject()
            
    # Save the project
    def saveProject(self):
        with open(self.currentConfigFileName, 'w') as f:
            self.currentConfig.write(f)
        
#         if self.currentProject is False:
#             temp = self.currentConfigFileName
#             self.sphyLocationPath = temp.split(":")[0] + ":"

        
        self.settings.setValue("sphyPreProcessPlugin/currentConfig", self.currentConfigFileName)
        self.projectDir = os.path.dirname(self.currentConfigFileName[0])
        
#         # write the qgs project file
#         qgsProjectFileName = ((self.currentConfigFileName).split(".cfg")[0]) + ".qgs"
#         qgsProject = QgsProject.instance()
#         qgsProject.setFileName(qgsProjectFileName)
#         qgsProject.write()
#         self.settings.setValue("sphyplugin/qgsProjectFileName", qgsProjectFileName.replace("\\","/"))
        
        # update tab, project settings, and gui
        self.currentProject = True
        self.Tab.setEnabled(1)
        #self.exitDate = True
        self.initGuiConfigMap()
        #self.exitDate = False
        self.updateSaveButtons(0)
        
    # enable save buttons after gui has been modified
    def updateSaveButtons(self, flag):
        if self.currentProject:
            self.saveAsButton.setEnabled(1)
#             self.saveButton.setEnabled(flag) 
            
