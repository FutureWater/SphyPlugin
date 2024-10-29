"""
Model exported as python.
Name : glaciers_model
Group : 
With QGIS : 33409
"""

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


class Glaciers_model(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer('clone_map', 'Clone Map', defaultValue=None))
        self.addParameter(QgsProcessingParameterVectorLayer('rgi_shapefile', 'RGI Shapefile', types=[QgsProcessing.TypeVectorPolygon], defaultValue=None))
        self.addParameter(QgsProcessingParameterRasterLayer('debris_tiff', 'Debris Tiff', defaultValue=None))
        self.addParameter(QgsProcessingParameterRasterLayer('dem', 'DEM', defaultValue=None))
        self.addParameter(QgsProcessingParameterRasterLayer('ferrinoti_tiff', 'Ferrinoti Tiff', defaultValue=None))
        self.addParameter(QgsProcessingParameterNumber('model_resolution', 'Model Resolution', type=QgsProcessingParameterNumber.Integer, defaultValue=1000))
        self.addParameter(QgsProcessingParameterCrs('model_crs', 'Model CRS', defaultValue='ESRI:102025'))
        self.addParameter(QgsProcessingParameterNumber('finer_resolution', 'Finer Resolution', type=QgsProcessingParameterNumber.Double, defaultValue=100))
        self.addParameter(QgsProcessingParameterFile('output_folder', 'Output folder', behavior=QgsProcessingParameterFile.Folder, fileFilter='All files (*.*)', defaultValue='E:\\amelia\\SPHY_demo\\glaciers_tutorial\\qgis_model\\outputs\\'))
        self.addParameter(QgsProcessingParameterFeatureSink('Glaciers', 'glaciers', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Rgi_clipped_reproject_glac_id', 'RGI_Clipped_Reproject_GLAC_ID', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Intersection_glaciers_uid', 'intersection_glaciers_uid', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Ice_depth', 'ICE_DEPTH', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Debris', 'DEBRIS', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Frac_glac', 'FRAC_GLAC', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Mod_id', 'MOD_ID', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Modid_int_glacid', 'MODID_int_GLACID', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('U_id', 'U_ID', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Modid_int_glacid_inclmodh', 'MODID_int_GLACID_inclMODH', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Intersection_glaciers_uid_hglac', 'intersection_glaciers_uid_HGLAC', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Debris_geom', 'DEBRIS_geom', type=QgsProcessing.TypeVectorAnyGeometry, createByDefault=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(30, model_feedback)
        results = {}
        outputs = {}

        # Create U_ID grid
        alg_params = {
            'CRS': parameters['model_crs'],
            'EXTENT': '-95331.0000,119169.0000,3213564.0000,3443814.0000', #QgsExpression("layer_property(@clone_map,'extent') ").evaluate(),
            'HOVERLAY': 0,
            'HSPACING': parameters['finer_resolution'],
            'TYPE': 2,  # Rectangle (Polygon)
            'VOVERLAY': 0,
            'VSPACING': parameters['finer_resolution'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['CreateU_idGrid'] = processing.run('native:creategrid', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Reproject DEM
        alg_params = {
            'DATA_TYPE': 0,  # Use Input Layer Data Type
            'EXTRA': '',
            'INPUT': parameters['dem'],
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'RESAMPLING': 0,  # Nearest Neighbour
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:4326'),
            'TARGET_CRS': parameters['model_crs'],
            'TARGET_EXTENT': QgsExpression("layer_property(@clone_map,'extent') ").evaluate(),
            'TARGET_EXTENT_CRS': parameters['model_crs'],
            'TARGET_RESOLUTION': None,
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ReprojectDem'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Reproject RGI
        alg_params = {
            'CONVERT_CURVED_GEOMETRIES': False,
            'INPUT': parameters['rgi_shapefile'],
            'OPERATION': '',
            'TARGET_CRS': parameters['model_crs'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ReprojectRgi'] = processing.run('native:reprojectlayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Create grid
        alg_params = {
            'CRS': parameters['model_crs'],
            'EXTENT': QgsExpression("layer_property(@clone_map,'extent') ").evaluate(),
            'HOVERLAY': 0,
            'HSPACING': parameters['model_resolution'],
            'TYPE': 2,  # Rectangle (Polygon)
            'VOVERLAY': 0,
            'VSPACING': parameters['model_resolution'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['CreateGrid'] = processing.run('native:creategrid', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Clip RGI to clone
        alg_params = {
            'EXTENT': QgsExpression("layer_property(@clone_map,'extent') ").evaluate(),
            'INPUT': outputs['ReprojectRgi']['OUTPUT'],
            'OPTIONS': '',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ClipRgiToClone'] = processing.run('gdal:clipvectorbyextent', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}

        # Reproject Debris
        alg_params = {
            'DATA_TYPE': 0,  # Use Input Layer Data Type
            'EXTRA': '',
            'INPUT': parameters['debris_tiff'],
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'RESAMPLING': 0,  # Nearest Neighbour
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:4326'),
            'TARGET_CRS': parameters['model_crs'],
            'TARGET_EXTENT': QgsExpression("layer_property(@clone_map,'extent')\n\n\n\n\n\n\n\n").evaluate(),
            'TARGET_EXTENT_CRS': parameters['model_crs'],
            'TARGET_RESOLUTION': None,
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ReprojectDebris'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        # Rename U_ID
        alg_params = {
            'FIELD': 'id',
            'INPUT': outputs['CreateU_idGrid']['OUTPUT'],
            'NEW_NAME': 'U_ID',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameU_id'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(7)
        if feedback.isCanceled():
            return {}

        # Retain U_ID
        alg_params = {
            'FIELDS': ['U_ID'],
            'INPUT': outputs['RenameU_id']['OUTPUT'],
            'OUTPUT': parameters['U_id']
        }
        outputs['RetainU_id'] = processing.run('native:retainfields', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['U_id'] = outputs['RetainU_id']['OUTPUT']

        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}

        # Compute GLACID
        alg_params = {
            'FIELD_LENGTH': 0,
            'FIELD_NAME': 'GLAC_ID',
            'FIELD_PRECISION': 0,
            'FIELD_TYPE': 0,  # Decimal (double)
            'FORMULA': '@row_number + 1',
            'INPUT': outputs['ClipRgiToClone']['OUTPUT'],
            'OUTPUT': parameters['Rgi_clipped_reproject_glac_id']
        }
        outputs['ComputeGlacid'] = processing.run('native:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Rgi_clipped_reproject_glac_id'] = outputs['ComputeGlacid']['OUTPUT']

        feedback.setCurrentStep(9)
        if feedback.isCanceled():
            return {}

        # Rename field
        alg_params = {
            'FIELD': 'id',
            'INPUT': outputs['CreateGrid']['OUTPUT'],
            'NEW_NAME': 'MOD_ID',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameField'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(10)
        if feedback.isCanceled():
            return {}

        # Retain MOD_ID
        alg_params = {
            'FIELDS': ['MOD_ID'],
            'INPUT': outputs['RenameField']['OUTPUT'],
            'OUTPUT': parameters['Mod_id']
        }
        outputs['RetainMod_id'] = processing.run('native:retainfields', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Mod_id'] = outputs['RetainMod_id']['OUTPUT']

        feedback.setCurrentStep(11)
        if feedback.isCanceled():
            return {}

        # Intersection
        alg_params = {
            'GRID_SIZE': None,
            'INPUT': outputs['ComputeGlacid']['OUTPUT'],
            'INPUT_FIELDS': [''],
            'OVERLAY': outputs['RetainMod_id']['OUTPUT'],
            'OVERLAY_FIELDS': [''],
            'OVERLAY_FIELDS_PREFIX': '',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['Intersection'] = processing.run('native:intersection', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(12)
        if feedback.isCanceled():
            return {}

        # Retain MOD_ID - GLACID
        alg_params = {
            'FIELDS': ['MOD_ID','GLAC_ID'],
            'INPUT': outputs['Intersection']['OUTPUT'],
            'OUTPUT': parameters['Modid_int_glacid']
        }
        outputs['RetainMod_idGlacid'] = processing.run('native:retainfields', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Modid_int_glacid'] = outputs['RetainMod_idGlacid']['OUTPUT']

        feedback.setCurrentStep(13)
        if feedback.isCanceled():
            return {}

        # Mean DEM
        alg_params = {
            'COLUMN_PREFIX': 'MOD_H_',
            'INPUT': outputs['RetainMod_idGlacid']['OUTPUT'],
            'INPUT_RASTER': outputs['ReprojectDem']['OUTPUT'],
            'RASTER_BAND': 1,
            'STATISTICS': [2],  # Mean
            'OUTPUT': parameters['Modid_int_glacid_inclmodh']
        }
        outputs['MeanDem'] = processing.run('native:zonalstatisticsfb', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Modid_int_glacid_inclmodh'] = outputs['MeanDem']['OUTPUT']

        feedback.setCurrentStep(14)
        if feedback.isCanceled():
            return {}

        # Intersection Finer Grid
        alg_params = {
            'GRID_SIZE': None,
            'INPUT': outputs['MeanDem']['OUTPUT'],
            'INPUT_FIELDS': [''],
            'OVERLAY': outputs['RetainU_id']['OUTPUT'],
            'OVERLAY_FIELDS': [''],
            'OVERLAY_FIELDS_PREFIX': '',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['IntersectionFinerGrid'] = processing.run('native:intersection', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(15)
        if feedback.isCanceled():
            return {}

        # Drop field UID
        alg_params = {
            'COLUMN': ['U_ID'],
            'INPUT': outputs['IntersectionFinerGrid']['OUTPUT'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['DropFieldUid'] = processing.run('native:deletecolumn', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(16)
        if feedback.isCanceled():
            return {}

        # Recompute U_ID
        alg_params = {
            'FIELD_LENGTH': 0,
            'FIELD_NAME': 'U_ID',
            'FIELD_PRECISION': 0,
            'FIELD_TYPE': 0,  # Decimal (double)
            'FORMULA': '@row_number + 1',
            'INPUT': outputs['DropFieldUid']['OUTPUT'],
            'OUTPUT': parameters['Intersection_glaciers_uid']
        }
        outputs['RecomputeU_id'] = processing.run('native:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Intersection_glaciers_uid'] = outputs['RecomputeU_id']['OUTPUT']

        feedback.setCurrentStep(17)
        if feedback.isCanceled():
            return {}

        # Zonal statistics GLAC_H
        alg_params = {
            'COLUMN_PREFIX': 'GLAC_H_',
            'INPUT': outputs['RecomputeU_id']['OUTPUT'],
            'INPUT_RASTER': outputs['ReprojectDem']['OUTPUT'],
            'RASTER_BAND': 1,
            'STATISTICS': [2],  # Mean
            'OUTPUT': parameters['Intersection_glaciers_uid_hglac']
        }
        outputs['ZonalStatisticsGlac_h'] = processing.run('native:zonalstatisticsfb', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Intersection_glaciers_uid_hglac'] = outputs['ZonalStatisticsGlac_h']['OUTPUT']

        feedback.setCurrentStep(18)
        if feedback.isCanceled():
            return {}

        # Zonal statistics ICE_DEPTH
        alg_params = {
            'COLUMN_PREFIX': 'ICE_DEPTH_',
            'INPUT': outputs['ZonalStatisticsGlac_h']['OUTPUT'],
            'INPUT_RASTER': parameters['ferrinoti_tiff'],
            'RASTER_BAND': 1,
            'STATISTICS': [2],  # Mean
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ZonalStatisticsIce_depth'] = processing.run('native:zonalstatisticsfb', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(19)
        if feedback.isCanceled():
            return {}

        # FILL ICE_DEPTH
        alg_params = {
            'FIELD_LENGTH': 0,
            'FIELD_NAME': 'ICE_DEPTH_mean',
            'FIELD_PRECISION': 0,
            'FIELD_TYPE': 0,  # Decimal (double)
            'FORMULA': 'if("ICE_DEPTH_mean" is NULL, 25, "ICE_DEPTH_mean")\r\n',
            'INPUT': outputs['ZonalStatisticsIce_depth']['OUTPUT'],
            'OUTPUT': parameters['Ice_depth']
        }
        outputs['FillIce_depth'] = processing.run('native:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Ice_depth'] = outputs['FillIce_depth']['OUTPUT']

        feedback.setCurrentStep(20)
        if feedback.isCanceled():
            return {}

        # Zonal statistics DEBRIS
        alg_params = {
            'COLUMN_PREFIX': 'DEBRIS_',
            'INPUT': outputs['FillIce_depth']['OUTPUT'],
            'INPUT_RASTER': outputs['ReprojectDebris']['OUTPUT'],
            'RASTER_BAND': 1,
            'STATISTICS': [9],  # Majority
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ZonalStatisticsDebris'] = processing.run('native:zonalstatisticsfb', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(21)
        if feedback.isCanceled():
            return {}

        # Fill DEBRIS
        alg_params = {
            'FIELD_LENGTH': 0,
            'FIELD_NAME': 'DEBRIS_majority',
            'FIELD_PRECISION': 0,
            'FIELD_TYPE': 0,  # Decimal (double)
            'FORMULA': 'if("DEBRIS_majority" is NULL, 0, "DEBRIS_majority")',
            'INPUT': outputs['ZonalStatisticsDebris']['OUTPUT'],
            'OUTPUT': parameters['Debris']
        }
        outputs['FillDebris'] = processing.run('native:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Debris'] = outputs['FillDebris']['OUTPUT']

        feedback.setCurrentStep(22)
        if feedback.isCanceled():
            return {}

        # Add geometry attributes
        alg_params = {
            'CALC_METHOD': 0,  # Layer CRS
            'INPUT': outputs['FillDebris']['OUTPUT'],
            'OUTPUT': parameters['Debris_geom']
        }
        outputs['AddGeometryAttributes'] = processing.run('qgis:exportaddgeometrycolumns', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Debris_geom'] = outputs['AddGeometryAttributes']['OUTPUT']

        feedback.setCurrentStep(23)
        if feedback.isCanceled():
            return {}

        # Compute GLAC FRAC
        alg_params = {
            'FIELD_LENGTH': 0,
            'FIELD_NAME': 'FRAC_GLAC',
            'FIELD_PRECISION': 0,
            'FIELD_TYPE': 0,  # Decimal (double)
            'FORMULA': '$area/((@model_resolution / @finer_resolution) * @model_resolution * @finer_resolution) \r\n',
            'INPUT': outputs['AddGeometryAttributes']['OUTPUT'],
            'OUTPUT': parameters['Frac_glac']
        }
        outputs['ComputeGlacFrac'] = processing.run('native:fieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Frac_glac'] = outputs['ComputeGlacFrac']['OUTPUT']

        feedback.setCurrentStep(24)
        if feedback.isCanceled():
            return {}

        # Rename final table
        alg_params = {
            'FIELD': 'MOD_H_mean',
            'INPUT': outputs['ComputeGlacFrac']['OUTPUT'],
            'NEW_NAME': 'MOD_H',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameFinalTable'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(25)
        if feedback.isCanceled():
            return {}

        # Rename field 2
        alg_params = {
            'FIELD': 'GLAC_H_mean',
            'INPUT': outputs['RenameFinalTable']['OUTPUT'],
            'NEW_NAME': 'GLAC_H',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameField2'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(26)
        if feedback.isCanceled():
            return {}

        # Rename field 3
        alg_params = {
            'FIELD': 'ICE_DEPTH_mean',
            'INPUT': outputs['RenameField2']['OUTPUT'],
            'NEW_NAME': 'ICE_DEPTH',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameField3'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(27)
        if feedback.isCanceled():
            return {}

        # Rename field 4
        alg_params = {
            'FIELD': 'DEBRIS_majority',
            'INPUT': outputs['RenameField3']['OUTPUT'],
            'NEW_NAME': 'DEBRIS',
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['RenameField4'] = processing.run('native:renametablefield', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(28)
        if feedback.isCanceled():
            return {}

        # Drop field(s)
        alg_params = {
            'COLUMN': ['area','perimeter'],
            'INPUT': outputs['RenameField4']['OUTPUT'],
            'OUTPUT': parameters['Glaciers']
        }
        outputs['DropFields'] = processing.run('native:deletecolumn', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Glaciers'] = outputs['DropFields']['OUTPUT']

        feedback.setCurrentStep(29)
        if feedback.isCanceled():
            return {}

        # Export to spreadsheet
        alg_params = {
            'FORMATTED_VALUES': False,
            'LAYERS': outputs['DropFields']['OUTPUT'],
            'OUTPUT': QgsExpression(" concat(@output_folder,'/glaciers.csv')").evaluate(),
            'OVERWRITE': True,
            'USE_ALIAS': False,
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['ExportToSpreadsheet'] = processing.run('native:exporttospreadsheet', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        return results

    def name(self):
        return 'glaciers_model'

    def displayName(self):
        return 'glaciers_model'

    def group(self):
        return ''

    def groupId(self):
        return ''

    def createInstance(self):
        return Glaciers_model()
