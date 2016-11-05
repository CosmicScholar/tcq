# mgtc
A fortran program for **M**arine **G**ravity **T**errain **C**orrection

# Synopsis
**mgtc** [**-C**_input.txt_] [**-I**_inner.grd3_] [**-O**_outer.grd3_] [**-D**_radius_inner_] [**-R**_radius_outer_] [**-T** [**0**|**1**] ] [**-G**_output.txt_]

# Description
**mgtc** calculate gravitational effect of terrain( or bathymetry) and add it to free-air gravity anomaly, then we can obtain Bouguer anomaly. The following three files must be available: **input.txt**, **inner.grd3** and **outer.grd3**.

# Required Arguments
_**input.txt**_ 
 
　　The file of coordinates [**longitude** **latitude**  **elevation**] (on land is positive and in ocean is negetive) of calculated points.

_**inner.grd3**_

　　A binary file of more refined bathymetry for computing the inner zone effect, which can be convered from surfer grd file using [grd2xyz](http://gmt.soest.hawaii.edu/doc/latest/grd2xyz.html) and [z2grd3](). For example: ``
gmt grd2xyz %outgrd% -Z |z2grd3 -R%outerrange% -I%dxdy_outer%/%dxdy_outer% -Gouter.grd3
``

_**outer.grd3**_

　　A binary file of topography (such as [ETOPO1](http://www.ngdc.noaa.gov/)) for computing the outer zone effect, which can be convered from surfer grd file using [grd2xyz](http://gmt.soest.hawaii.edu/doc/latest/grd2xyz.html) and [z2grd3](). Transform method as above.

_**-D**_

　　Radius (in km) of inner zone. [[Hwang et al., 2003]](http://xueshu.baidu.com/s?wd=paperuri%3A%28d8ea03100103cf353252ce3c661f0b31%29&filter=sc_long_sign&tn=SE_xueshusource_2kduw22v&sc_vurl=http%3A%2F%2Fwww.sciencedirect.com%2Fscience%2Farticle%2Fpii%2FS0098300403001560&ie=utf-8&sc_us=16355709148949826115) suggent that the best inner radii for terrain correction computation is 20 km.

_**-R**_

　　Radius (in km) of outer zone. [[Hwang et al., 2003]](http://xueshu.baidu.com/s?wd=paperuri%3A%28d8ea03100103cf353252ce3c661f0b31%29&filter=sc_long_sign&tn=SE_xueshusource_2kduw22v&sc_vurl=http%3A%2F%2Fwww.sciencedirect.com%2Fscience%2Farticle%2Fpii%2FS0098300403001560&ie=utf-8&sc_us=16355709148949826115) suggent that the best outer radii for terrain correction computation is 200 km.

_**-G**_

　　Output file of *longitude* *latitude* *bathymetry* ***innermost* *inner* *outer* *total* *slab*** effect.

# Optional Arguments

_**-T**_

　　**0** for terrain correction at sea level. **1** for terrain correction at terrain surface. The latter is suitable for land.


# Example

## 　1. Topography data for inner and outer zone respectively

### 　　(1) Create inner.grd3

```batch
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
```

### 　　(2) Create outer.grd3

```batch
REM 对网格数据重新采样，给定范围个间距，为了使两个网格数据范围和
REM 网格节点一一对应
set outerrange=44/56/-42/-32
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
```

## 　2. Create input.txt file for calculated point

```batch
REM run tcq to compute terrain correction

REM inputdata
set dx_input=0.2m
set dy_input=0.2m
set inputrange=49.08/50.92/-38/-37.3

set intopo=G:\Research\5-Gravity_In_Ga\data\dragon_grd\geoCoord\bathy_dragon
set ext=.grd
set inputgrd=%inputpath%input.grd

gmt grdsample %intopo%%ext% -R%inputrange% -nb -I%dx_input%/%dy_input% -G%inputgrd%
gmt grd2xyz %inputgrd% >%inputpath%input.txt
gmt grdinfo -L1 -L2 %inputgrd% >%inputgrd%_info.txt

REM cut the faa from the same rigion and the same nodes
set faa=%inputpath%faa.grd
set infaa=G:\Research\5-Gravity_In_Ga\data\dragon_grd\geoCoord\faa_dragon
gmt grdsample %infaa%%ext% -R%inputrange% -nb -I%dx_input%/%dy_input% -G%faa%
gmt grdinfo -L1 -L2 %faa% >%faa%_info.txt

del fort.60 gmt.history 
pause
```

## 　3. Run mgtc and output the results
The parameters of ``dx_input`` ``dy_input`` and ``inputrange`` must the same as above (**Create input.txt file**)
```batch
REM run tcq to compute terrain correction
set inputpath=inputdata\
set outputpath=outputdata\
set faa=%inputpath%faa.grd
set dx_input=0.2m
set dy_input=0.2m
set inputrange=49.08/50.92/-38/-37.3
REM del %outputpath%*.*
REM run
set Rdistance=200
set output=%outputpath%topo_grav_R%Rdistance%.txt

tcq -C%inputpath%input.txt -I%inputpath%inner.grd3 -O%inputpath%outer.grd3 -D20 -R%Rdistance% -T0 -G%output%

gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,2 -G%outputpath%Bathy.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,3 -G%outputpath%TC_Inmost_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,4 -G%outputpath%TC_Inner_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,5 -G%outputpath%TC_Outer_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,6 -G%outputpath%TC_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,7 -G%outputpath%slab.grd
REM compute bouguer anomaly
gmt grdmath %faa% %outputpath%TC_R%Rdistance%.grd ADD = %outputpath%Boug_tcq_R%Rdistance%.grd

del fort.60 gmt.history 
pause
```