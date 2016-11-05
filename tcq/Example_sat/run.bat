REM run tcq to compute terrain correction
set inputpath=inputdata\
set outputpath=outputdata\
set faa=%inputpath%faa.grd
set inputrange=49.6/50.7/-37.9/-37.4
set dx_input=1m
set dy_input=1m
del %outputpath%*.*
REM run
set Rdistance=200
set output=%outputpath%topo_grav_R%Rdistance%.txt

tcq -C%inputpath%input.txt -I%inputpath%inner.grd3 -O%inputpath%outer.grd3 -D20 -R%Rdistance% -T0 -G%output%

REM convert result
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,3 -G%outputpath%TC_Inmost_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,4 -G%outputpath%TC_Inner_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,5 -G%outputpath%TC_Outer_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,6 -G%outputpath%TC_R%Rdistance%.grd
gmt xyz2grd %output% -R%inputrange% -I%dx_input%/%dy_input% -i0,1,7 -G%outputpath%slab.grd
REM compute bouguer anomaly
gmt grdmath %faa% %outputpath%TC_R%Rdistance%.grd ADD = %outputpath%Boug_tcq_R%Rdistance%.grd

del fort.60 gmt.history 
pause