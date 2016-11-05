REM run tcq to compute terrain correction

REM inputdata
set dx_input=1m
set dy_input=1m
set inputrange=49.6/50.7/-37.9/-37.4
set ingrd=G:\Data\1-GlobalGeophysicalData\ETOPO1\ETOPO1_Bed_g_gdal.grd\ETOPO1_Bed_g_gdal
set ext=.nc
set inputgrd=%inputpath%input.grd
gmt grdsample %ingrd%%ext% -R%inputrange% -nb -I%dx_input%/%dy_input% -G%inputgrd%
gmt grd2xyz %inputgrd% >%inputpath%input.txt
REM cut the faa from the same rigion and the same nodes
set faa=%inputpath%faa.grd
set faagrddata=G:\Data\1-GlobalGeophysicalData\BGA-WGM2012\WGM2012_Freeair_ponc_2min
gmt grdsample %faagrddata%.grd -R%inputrange% -nb -I%dx_input%/%dy_input% -G%faa%


del fort.60 gmt.history 
pause