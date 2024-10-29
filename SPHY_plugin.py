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
from __future__ import absolute_import
from builtins import object
__author__ = "Wilco Terink"
__copyright__ = "Wilco Terink"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "w.terink@futurewater.nl, terinkw@gmail.com"
__date__ ='1 January 2017'
############################################################################################

"""
/***************************************************************************
 SphyPlugin
                                 A QGIS plugin
 A tool to convert raw data into SPHY model input data
                              -------------------
        begin                : 2015-06-23
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Wilco Terink
        email                : w.terink@futurewater.nl
 ***************************************************************************/

"""
from qgis.PyQt.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtGui import QIcon
# Initialize Qt resources from file resources.py
from SphyPlugin.gui.generated import resources_rc
# Import the code for the dialog
from .SPHY_plugin_dialog import SphyPluginDialog
import os.path

try:
    from pydevd import *
except:
    None

class SphyPlugin(object):
    def __init__(self, iface):
        self.iface = iface
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'SphyPreProcess_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = SphyPluginDialog()

    def initGui(self):
        self.action = QAction(QIcon(":/plugins/SphyPlugin/images/icon.png"), u"SPHY PLUGIN", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu(u"&SPHY PLUGIN", self.action)


    def unload(self):
        self.iface.removePluginMenu(u"&SPHY PLUGIN", self.action)
        self.iface.removeToolBarIcon(self.action)
 
    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
