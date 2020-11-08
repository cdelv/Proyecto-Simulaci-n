//main.cpp

#include "master.h"

int main(void)
{
  //----------------parametros del propagador--------------------
  
  int N=2;                    //numero de orbitas
  earth cb;                  //cuerpo central
  perturbations perts;      //diccionario de perturbaciones
  double tspan=3600*24*5;  //tmax en segundos
  double dt=10;            // paso de tiempo en segundos
  bool coes=true;        //si condiciones iniciales son coes
  bool deg=true;        //si están en grados
  double masa=50;      //kg
  int tcuadro=10;    //cada cuanto imprime datos
 
  //------------------definir las perturbaciones----------------
  
  //perts.J2=true;
  //perts.aero=true; perts.Cd=2.2; perts.A=std::pow(0.001,2); //Cd=coeficiente de fricción A=area en km²
  perts.thrust=0.327; perts.isp=4300; perts.thrust_direction=1;
  
  //------------------definir el propagador---------------------
  
  std::vector <OrbitPropagator> OP(N);

  //-----------------definir condiciones iniciales-------------
  
  // si coes=false state0{rx,ry,rx,vx,vy,vz}
  //si coes=true   stateo{a,e,i,ta,aop,raan}

  double r0=cb.radius+414;
  std::vector <double> state0{r0,0.0006189,51.6393,0.0,234.1955,105.6372};
 std::vector <double>state1= tlecoes("HJ-2A.txt",cb); //sacar condiciones iniciales de los TLE
  
 //------------------parametros de Plot_orbit--------------------
  
  bool save=false; //guarda la gráfica
  bool show=true; //muestra la grafica al final de correr
  std::vector<std::string> labels {"iss","HJ-2A"}; //labels de los propagadores
  std::string title_orbit ="Many Orbits"; //titulo de la gráfica
  
 //---------------------- Propagar orbitas ----------------------

  //si hay varios propagadores el nombre "OP" debe ser distinto
  //es el archivo .dat donde se guardan los datos
  
  OP[0].inicie(state0,tspan,dt,"OP",cb,coes,deg,perts,masa,tcuadro);
  OP[1].inicie(state1,tspan,dt,"OP1",cb,coes,deg,perts,masa,tcuadro);

  //-----------------------pintar las orbitas--------------------

  
//Plot_orbit(el vector de propagadores,cuerpo central ,"título", save);
  Plot_orbit_gnuplot(OP,cb,"N orbits", save); //grafica bonita pero sin ejes ni labels
  
//Plot_orbit(OP,cb,show,save,title_orbit,labels); //grafica en python con labeles y ejes
  
  //  Plot_coes();
  
  return 0;
}
