r = 695700
R = 695701
#color definitions
set border lw 1.5
set style line 1 lc rgb '#000000' lt 1 lw 0.5
set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5
unset key; unset border; set tics scale 0
set format ''
set angles degrees
set xyplane at -1
set view 56,81
set lmargin screen 0
set bmargin screen 0
set rmargin screen 1
set tmargin screen 1
set parametric
set urange[0:360]
set vrange[-90:90]
set isosamples 25
set hidden3d
#since we are using Cartesian coordinates, we don't want this
#set mapping spherical
splot   r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 2, 'OP.dat' u 1:2:3 w l lw 5 lc 5, 'OP1.dat' u 1:2:3 w l lw 5 lc 10
pause mouse
