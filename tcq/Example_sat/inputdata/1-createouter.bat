REM 对网格数据重新采样，给定范围个间距，为了使两个网格数据范围和
REM 网格节点一一对应
set outerrange=45/55/-42/-32
set dxdy_outer=0.016666666666667
set ingrd=G:\Data\1-GlobalGeophysicalData\ETOPO1\ETOPO1_Bed_g_gdal.grd\ETOPO1_Bed_g_gdal
set ext=.nc

set outgrd=outer.grd
gmt grdsample %ingrd%%ext% -R%outerrange% -nb -I%dxdy_outer% -G%outgrd%
gmt grdinfo -L1 -L2 %outgrd% >%outgrd%_info.txt

REM transform outer grd to grd3 format
gmt grd2xyz %outgrd% -Z |z2grd3 -R%outerrange% -I%dxdy_outer%/%dxdy_outer% -Gouter.grd3

REM plot
del %outgrd%.ps
gmtset PSIMAGE_FORMAT hex
gmtset HEADER_FONT_SIZE 16p
gmtset LABEL_FONT_SIZE=12p
gmtset ANNOT_FONT_SIZE_PRIMARY = 12p
gmt grd2cpt %outgrd% -Crainbow -Z >tmp.cpt
gmt grdimage %outgrd% -JM12c -R%outerrange% -Ctmp.cpt -I -Ba -BWseN+t"Outer Zone" -K >%outgrd%.ps
gmt grdimage inner.grd -J -R -Cbathy -I -O -K >>%outgrd%.ps
gmt psbasemap -J -R -Lfx10c/0.8c/49/-37.5/200+l+fgrey@20+jr -O -K >>%outgrd%.ps
gmt psscale -Ctmp.cpt -Dx6c/-0.4c+w12c/0.5c+jTC+h -B1000 -By+lm -O -P >> %outgrd%.ps
ps2raster %outgrd%.ps -AU -P -Tg

del gmt.history tmp.cpt gmt.conf *.eps
gsview32 %outgrd%.ps
pause