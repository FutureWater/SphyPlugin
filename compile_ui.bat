@ECHO OFF

::Poner la ruta de la instalacion de OsGeo
set OSGEO4W_ROOT=C:\OSGeo4W

set PATH=%OSGEO4W_ROOT%\bin;%PATH%
set PATH=%PATH%;%OSGEO4W_ROOT%\apps\qgis-ltr\bin

@echo off
call "%OSGEO4W_ROOT%\bin\o4w_env.bat"
::call "%OSGEO4W_ROOT%\bin\qt5_env.bat"
::call "%OSGEO4W_ROOT%\bin\py3_env.bat"
@echo off
path %OSGEO4W_ROOT%\apps\qgis-ltr\bin;%OSGEO4W_ROOT%\apps\grass\grass-83\lib;%OSGEO4W_ROOT%\apps\grass\grass-83\bin;%PATH%

::cd /c %~dp0

@ECHO ON
::Ui Compilation
call pyuic5 --import-from SphyPlugin.gui.generated ui\SPHY_plugin_dialog_base.ui -o gui\generated\SPHY_plugin_dialog_base.py    

::Resources
call pyrcc5 ui\resources.qrc -o gui\generated\resources_rc.py

@ECHO OFF
GOTO END

:ERROR
   echo "Failed!"
   set ERRORLEVEL=%ERRORLEVEL%
   pause

:END
@ECHO ON

@PAUSE
 