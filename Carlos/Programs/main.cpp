//main.cpp

#include "master.h"



int main(void)
{
  int N=2;
  sun cb;
  std::vector <OrbitPropagator> OP(N);
  double r0=cb.radius+500;
  double omega=std::sqrt(cb.mu/(r0*r0*r0)),T=2*M_PI/omega, V0=r0*omega;
  std::vector <double> state0{r0,0,0,0,V0,0};
  std::vector <double> state1{r0,0,0,0,0,V0};
  double tspan=T;
  double dt=1;
  
  OP[0].inicie(state0,tspan,dt,"OP",cb);
  OP[1].inicie(state1,tspan,dt,"OP1",cb);

  Plot_orbit(OP,cb,"N orbits");
  
  return 0;
}
