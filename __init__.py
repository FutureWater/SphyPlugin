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

#-Authorship information-###################################################################
__author__ = "FutureWater"
__copyright__ = "FutureWater"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "sphy@futurewater.nl"
__date__ ='1 February 2025'
############################################################################################
"""
/***************************************************************************
 SphyPreProcess
                                 A QGIS plugin
 A tool to convert raw data into SPHY model input data
                             -------------------
        begin                : 2015-06-23
        copyright            : (C) 2015 by Wilco Terink
        email                : w.terink@futurewater.nl
        git sha              : $Format:%H$
 ***************************************************************************/

 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load SphyPreProcess class from file SphyPreProcess.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .SPHY_plugin import SphyPlugin
    return SphyPlugin(iface)
