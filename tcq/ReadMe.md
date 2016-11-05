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
