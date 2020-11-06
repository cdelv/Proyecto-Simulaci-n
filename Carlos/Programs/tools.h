//tools.h
#pragma once
#include "master.h"

const double d2r=M_PI/180;


template <typename T>
void Plot_orbit(std::vector<OrbitPropagator> &OP,const T &cb, std::string title, bool save);
template <typename T>
std::vector<double> coes2rv(std::vector<double> &coes, bool deg, const T& cb);
double ecc_anomaly(double ta, double e, std::string method);

//--------------------------implementar funciones----------------
template <typename T>
void Plot_orbit(std::vector<OrbitPropagator> &OP,const T &cb, std::string title, bool save)
{
 double r=cb.radius;
 int t=r/1000;
 int j=0;
 std::ofstream out;
 out.open("plot.gp");

 if(save){
   out<<"set terminal pngcairo"<<std::endl;
   out<<"set output 'fig.png'"<<std::endl;
 }
 
 if(cb.name=="EARTH")
   {
     out<<"r = "<<cb.radius<<std::endl;
     out<<"R = "<<cb.radius+1<<std::endl;
     out<<"set title '"<<title<<"'"<<std::endl;
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
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 2,	\
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
     out<<"set title '"<<title<<"'"<<std::endl;
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


template <typename T>
std::vector<double> coes2rv(std::vector<double> &coes, bool deg, const T &cb)
{
  double a=0, e=0, i=0, ta=0, aop=0, raan=0, r_norm=0, E=0, rxp=0,ryp=0,rzp=0,vxp=0,vyp=0,vzp=0;
  std::vector <double> rv(6,0);
  std::vector <std::vector <double>> row(3,std::vector<double>(3,0));
  a=coes[0]; e=coes[1]; i=coes[2]; ta=coes[3]; aop=coes[4]; raan=coes[5];
    
  if(deg) //lo pasa a radianes
    {
      i*=d2r;
      ta*=d2r;
      aop*=d2r;
      raan*=d2r;
    }
  std::cout <<"a "<<a<<" e "<<e<<" i "<<i<<" ta "<<ta<<" aop "<<aop<<" raan "<<raan<<std::endl;
  
  E=ecc_anomaly(ta,e,"tae"); //calcula la anomalía con newton rapshon
  
  r_norm=a*(1-e*e)/(1+e*std::cos(ta));

  //calcular r y v en el sistema perifocal
  rxp=r_norm*std::cos(ta);  vxp=r_norm*std::sin(E)*std::sqrt(cb.mu*a)/r_norm;
  ryp=r_norm*std::sin(ta);  vyp=std::cos(E)*std::sqrt(1-e*e)*std::sqrt(cb.mu*a)/r_norm;
  rzp=0;                    vzp=0;
  std::cout <<"pxr "<<rxp<<std::endl;
  std::cout <<"pyr "<<ryp<<std::endl;
  std::cout <<"pxv "<<vxp<<std::endl;
  std::cout <<"pyv "<<vyp<<std::endl;
 
  std::vector<double> rperif{rxp,ryp,rzp};
  std::vector<double> vperif{vxp,vyp,vzp};
 
  //matriz de rotación 
  row[0]={-std::sin(raan)*std::cos(i)*std::sin(aop)+std::cos(raan)*std::cos(aop),std::cos(raan)*std::cos(i)*std::sin(aop)+std::sin(raan)*std::cos(aop),
	  std::sin(i)*std::sin(aop)};
  row[1]={-std::sin(raan)*std::cos(i)*std::cos(aop)-std::cos(raan)*std::sin(aop),std::cos(raan)*std::cos(i)*std::cos(aop)-std::sin(raan)*std::sin(aop),
	  std::sin(i)*std::cos(aop)};
  row[2]={std::sin(raan)*std::sin(i),-std::cos(raan)*std::sin(i),std::cos(i)};

    for(int jj=0; jj<3; jj++)
      for(int ii=0; ii<3; ii++)
	std::cout <<"vperif"<<jj<<","<<ii<<"="<<row[jj][ii]<<std::endl;

  //multiclicación de la matriz de rotacion TRANSPUESTA por rperif y vperif
  //para calcular r y v en el sistema inercial del cuerpo central
    
    for(int ii=0; ii<3; ii++)
      for(int jj=0; jj<3; jj++)
	rv[ii]+=rperif[jj]*row[jj][ii];

    
    for(int ii=0; ii<3; ii++)
    for(int jj=0; jj<3; jj++)
    rv[ii+3]+=vperif[jj]*row[jj][ii];
  
  
    /*for(int ii=0; ii<3; ii++)
    for(int jj=0; jj<3; jj++)
      rv[ii]+=row[ii][jj]*rperif[jj];
  for(int ii=0; ii<3; ii++)
    for(int jj=0; jj<3; jj++)
    rv[ii+3]+=row[ii][jj]*vperif[jj];*/

  for(int jj=0; jj<6; jj++)
    std::cout <<"rv["<<jj<<"]="<<rv[jj]<<std::endl;
  return rv;
}
double ecc_anomaly(double ta, double e, std::string method)
{
  double  Me=ta, E0=0, E1=0, ratio=0, tol=1e-8;
  
  //newton's method for iteratively finding E
  if (method=="newton")
    {      
      if (Me<M_PI/2.0){
	E0=Me+e/2.0;}
      else{
	E0=Me-e;}
      
      for(int n=0; n<200; n++){
	
	ratio=(E0-e*std::sin(E0)-Me)/(1-e*std::cos(E0));
	if (abs(ratio)<tol)
	  {
	    if (n==0){
	      return E0;}
	    else{
	      return E1;}
	  }
	else{
	  E1=E0-ratio;
	  E0=E1;
	    }
      }
    }
  
  if (method=="tae")
    return 2*std::atan(std::sqrt((1-e)/(1+e))*std::tan(ta/2.0));
  
  else
    std::cout <<"método invalido para calcular la anomalía" << std::endl;
  
  //did not converge
  std::cout << "did not converge" << std::endl;
  return 0;
}
