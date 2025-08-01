# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
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
from builtins import str
from builtins import object

############################################################################################

import os

#-Class with gdal processing commands
class SpatialProcessing(object):
    def __init__(self, infile, outfile, s_srs, t_srs, resolution, resampling='bilinear', rtype='Float32', extra=None, projwin=None):
        self.input = infile
        self.output = outfile
        self.s_srs = s_srs
        self.t_srs = t_srs
        self.t_res = str(resolution)
        self.resampling = resampling
        self.rtype = rtype
        self.extra = extra
        self.projwin = projwin
        
    #-Reproject, resample, and clip to extent
    def reproject(self):
        command = 'gdalwarp -s_srs ' + self.s_srs + ' -t_srs ' + self.t_srs + ' -r ' + self.resampling\
                        + ' -tr ' + self.t_res + ' ' + self.t_res +' -ot ' + self.rtype + ' ' + self.extra\
                        + ' ' + self.input + ' ' + self.output 
        return command
    
    #-Convert raster format
    def rasterTranslate(self):
        command = 'gdal_translate ' + self.extra + ' ' + self.input + ' ' + self.output
        return command
    
    #-Rasterize Vector
    def rasterize(self):
        d = os.path.dirname(self.input)
        fi = os.path.relpath(self.input, d)
        layer = fi.split('.shp')[0]
        command = 'gdal_rasterize ' + self.extra + ' -l ' + layer + ' -tr ' + self.t_res + ' ' + \
            self.t_res + ' ' + self.input + ' ' + self.output
        return command
    
