//main.cpp

#include "master.h"

int main(void)
{
  //----------------parametros del propagador--------------------
  
  int N=1;                    //numero de orbitas
  earth cb;                  //cuerpo central
  perturbations perts;      //diccionario de perturbaciones
  double tspan=3600*24*1;  //tmax en segundos
  double dt=1;            // paso de tiempo en segundos
  bool coes=true;        //si condiciones iniciales son coes
  bool deg=true;        //si están en grados
 
  //------------------definir las perturbaciones----------------
  
  perts.J2=true;
  
  //------------------definir el propagador---------------------
  
  std::vector <OrbitPropagator> OP(N);

  //-----------------definir condiciones iniciales-------------
  
  // si coes=false state0{rx,ry,rx,vx,vy,vz}
  //si coes=true   stateo{a,e,i,ta,aop,raan}

  double r0=cb.radius+414;
  std::vector <double> state0{r0,0.0006189,51.6393,0.0,234.1955,105.6372};
  std::vector <double> state1(6,0);
  std::vector <double> state2(6,0);

 //state1= tlecoes("iss.txt",cb); //sacar condiciones iniciales de los TLE
  state2=coes2rv(state0, deg, cb);
  
  state1=rv2coes(state2,cb,deg,false);
  
 //------------------parametros de Plot_orbit--------------------
  
  bool save=false; //guarda la gráfica
  
 //---------------------- Propagar orbitas ----------------------

  //si hay varios propagadores el nombre "OP" debe ser distinto
  //es el archivo .dat donde se guardan los datos
  
  OP[0].inicie(state1,tspan,dt,"OP",cb,coes,deg,perts);

  //-----------------------pintar las orbitas--------------------

  
//Plot_orbit(el vector de propagadores,cuerpo central ,"título", save);
  Plot_orbit(OP,cb,"N orbits", save);
  
  return 0;
}
