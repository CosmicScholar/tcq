REM 对网格数据重新采样，给定范围个间距，为了使两个网格数据范围和
REM 网格节点一一对应
set innerrange=48.6/51.69/-38.89/-36.4
set dx_inner=1m
set dy_inner=1m
set ingrd=G:\Data\1-GlobalGeophysicalData\ETOPO1\ETOPO1_Bed_g_gdal.grd\ETOPO1_Bed_g_gdal
set ext=.nc

set outgrd=inner.grd
gmt grdsample %ingrd%%ext% -R%innerrange% -nb -I%dx_inner%/%dy_inner% -G%outgrd%
gmt grdinfo -L1 -L2 %outgrd% >%outgrd%_info.txt

REM transform inner grd to grd3 format
gmt grd2xyz %outgrd% -Z |z2grd3 -R%innerrange% -I%dx_inner%/%dy_inner% -Ginner.grd3


REM plot
del %outgrd%.ps
gmtset PSIMAGE_FORMAT hex
gmtset HEADER_FONT_SIZE 16p
gmtset LABEL_FONT_SIZE=12p
gmtset ANNOT_FONT_SIZE_PRIMARY = 12p
gmt grd2cpt %outgrd% -Crainbow -Z >tmp.cpt
gmt grdimage %outgrd% -JM12c -R%innerrange% -Ctmp.cpt -I -Ba -BWseN+t"Inner Zone" -K >%outgrd%.ps
gmt psbasemap -J -R -Lfx10c/0.8c/49/-37.5/40+l+fgrey@20+jr -O -K >>%outgrd%.ps
gmt psscale -Ctmp.cpt -Dx6c/-0.4c+w12c/0.5c+jTC+h -B500 -By+lm -O -P >> %outgrd%.ps
ps2raster %outgrd%.ps -AU -P -Tg
gsview32 %outgrd%.ps
del gmt.history tmp.cpt gmt.conf *.eps
pause