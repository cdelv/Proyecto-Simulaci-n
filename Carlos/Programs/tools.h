//tools.h
#pragma once
#include "master.h"

template <typename T>
void Plot_orbit(std::vector<OrbitPropagator> &OP,const T &cb, std::string title);

//--------------------------implementar funciones----------------
template <typename T>
void Plot_orbit(std::vector<OrbitPropagator> &OP,const T &cb, std::string title)
{
 double r=cb.radius;
 int t=r/1000;
 int j=0;
 std::ofstream out;
 out.open("plot.gp");

 if(cb.name=="EARTH")
   {
     out<<"r = "<<cb.radius<<std::endl;
     out<<"R = "<<cb.radius+1<<std::endl;
     //out<<"set terminal pngcairo"<<std::endl;
     // out<<"set output 'fig.png'"<<std::endl;
     //out<<"set xr [-2*r:2*r]"<<std::endl;
     //out<<"set yr [-2*r:2*r]"<<std::endl;
     //out<<"set zr [-2*r:2*r]"<<std::endl;
     out<<"#color definitions"<<std::endl;
     out<<"set border lw 1.5"<<std::endl;
     out<<"set style line 1 lc rgb '#000000' lt 1 lw 0.5"<<std::endl;
     out<<"set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5"<<std::endl;
     out<<"unset key; unset border; set tics scale 0"<<std::endl;
     out<<"set format ''"<<std::endl;
     out<<"set angles degrees"<<std::endl;
     out<<"set xyplane at -1"<<std::endl;
     out<<"set view 56,81"<<std::endl;
     out<<"set lmargin screen 0"<<std::endl;
     out<<"set bmargin screen 0"<<std::endl;
     out<<"set rmargin screen 1"<<std::endl;
     out<<"set tmargin screen 1"<<std::endl;
     out<<"set parametric"<<std::endl;
     out<<"set urange[0:360]"<<std::endl;
     out<<"set vrange[-90:90]"<<std::endl;
     out<<"set isosamples 25"<<std::endl;
     out<<"set hidden3d"<<std::endl;
     out<<"#since we are using Cartesian coordinates, we don't want this"<<std::endl;
     out<<"#set mapping spherical"<<std::endl;
out<<"splot \
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 2, \
  'world_110m.txt' u (R*cos($1)*cos($2)):(R*sin($1)*cos($2)):(R*sin($2)) w l lw 3 lc rgb 'black'";
 for(auto i: OP){
   out <<", '"<< i.file<<".dat' u 1:2:3 w l lw 5 lc "<<5+j*5;
   j+=1;
 }
  out<<std::endl;
 out<<"pause mouse"<<std::endl;
   }
 else
   {
     out<<"r = "<<cb.radius<<std::endl;
     out<<"R = "<<cb.radius+1<<std::endl;
     //out<<"set terminal pngcairo"<<std::endl;
     // out<<"set output 'fig.png'"<<std::endl;
     //out<<"set xr [-2*r:2*r]"<<std::endl;
     //out<<"set yr [-2*r:2*r]"<<std::endl;
     //out<<"set zr [-2*r:2*r]"<<std::endl;
     out<<"#color definitions"<<std::endl;
     out<<"set border lw 1.5"<<std::endl;
     out<<"set style line 1 lc rgb '#000000' lt 1 lw 0.5"<<std::endl;
     out<<"set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5"<<std::endl;
     out<<"unset key; unset border; set tics scale 0"<<std::endl;
     out<<"set format ''"<<std::endl;
     out<<"set angles degrees"<<std::endl;
     out<<"set xyplane at -1"<<std::endl;
     out<<"set view 56,81"<<std::endl;
     out<<"set lmargin screen 0"<<std::endl;
     out<<"set bmargin screen 0"<<std::endl;
     out<<"set rmargin screen 1"<<std::endl;
     out<<"set tmargin screen 1"<<std::endl;
     out<<"set parametric"<<std::endl;
     out<<"set urange[0:360]"<<std::endl;
     out<<"set vrange[-90:90]"<<std::endl;
     out<<"set isosamples 25"<<std::endl;
     out<<"set hidden3d"<<std::endl;
     out<<"#since we are using Cartesian coordinates, we don't want this"<<std::endl;
     out<<"#set mapping spherical"<<std::endl;
out<<"splot \
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 2";
 for(auto i: OP){
   out <<", '"<< i.file<<".dat' u 1:2:3 w l lw 5 lc "<<5+j*5;
   j+=1;
 }
  out<<std::endl;
 out<<"pause mouse"<<std::endl;

   }
 
 out.close();
 system("gnuplot plot.gp");
}
