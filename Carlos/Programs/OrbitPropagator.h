//OrbitPropagator.h
#pragma once
#include "master.h"


class OrbitPropagator{
 private:
  vector3D r, v, a;
  double Tmax=0, dt=0;
  bool coes, deg;

 public:
  std::string file="OP";
  
  template <typename T>
  void inicie(std::vector <double> &state0, double tspan, double DT,std::string name,const  T &CB, bool COES, bool DEG, perturbations &perts);
  template <typename T>
  void propagate_orbit(const T &CB,perturbations &perts);
  template <typename T>
  void CalculeAceleracion(const T &CB,perturbations &perts);
  void Mueva_r(double t, double coeficiente);
  void Mueva_v(double t, double coeficiente);
  void BorreAceleracion(void){a.cargue(0,0,0);};
  double Getx(void){return r.x();};                                                                                           
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
};

//-------------------------Implementar funciones------------------------------

template <typename T>
void OrbitPropagator::inicie(std::vector <double> &state0, double tspan, double DT,std::string name,const T &CB, bool COES, bool DEG,perturbations &perts)
{
  coes=COES;
  deg=DEG;
  
  if (coes)
    {
      std::vector <double> state(6,0);
      state=coes2rv(state0,deg,CB);
      r.cargue(state[0],state[1],state[2]);
      v.cargue(state[3],state[4],state[5]);
    }
  else{
    r.cargue(state0[0],state0[1],state0[2]);
    v.cargue(state0[3],state0[4],state0[5]);
  }

  Tmax=tspan; dt=DT;
  file=name;
  propagate_orbit(CB,perts);
}
template <typename T>
void OrbitPropagator::CalculeAceleracion(const T &CB,perturbations &perts){
 BorreAceleracion();
 double aux=-CB.mu*std::pow(norma2(r),-1.5);
 a+=aux*r;
 if(perts.J2)
   {
     double z2=std::pow(r.z(),2);
     double r2=norma2(r);
     double norm_r=std::sqrt(r2);
     double tx=r.x()/norm_r*(5*z2/r2-1);
     double ty=r.y()/norm_r*(5*z2/r2-1);
     double tz=r.z()/norm_r*(5*z2/r2-3);
     aux=1.5*CB.J2*CB.mu*std::pow(CB.radius,2)/std::pow(norm_r,4);
     vector3D a_j2; a_j2.cargue(tx,ty,tz);
     a+=aux*a_j2;
   }
}
void OrbitPropagator::Mueva_r(double t, double coeficiente){
  r+=v*t*coeficiente;
}

void OrbitPropagator::Mueva_v(double t, double coeficiente){
  v+=a*t*coeficiente;
}
template <typename T>
void OrbitPropagator::propagate_orbit(const T &CB,perturbations &perts)
{
double E=0.1786178958448091e00;
double L=-0.2123418310626054e0;
double X=-0.6626458266981849e-1;
double coeficiente1=(1-2*L)/2;
double coeficiente2=(1-2*(X+E));
 
 std::ofstream fout;   //Salida de la orbita.                                                                                                                   
 fout.open(file+".dat");

  for(double t=0; t<=Tmax; t+=dt)
    {
    Mueva_r(dt,E);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,coeficiente1);
    Mueva_r(dt,X);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,L);
    Mueva_r(dt,coeficiente2);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,L);
    Mueva_r(dt,X);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,coeficiente1);
    Mueva_r(dt,E);
    
    fout << r.x() <<"\t"<< r.y() <<"\t"<< r.z() << std::endl;
    }
  fout.close();
}
