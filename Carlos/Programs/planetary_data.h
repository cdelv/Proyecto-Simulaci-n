//planetary_data.h
#pragma once
#include <cmath>

const double G_meters=6.67408e-11;
const double G=G_meters*std::pow(10,-9);

class sun {
 private:

 public:
  std::string name="SUN";
  double mass=1.989e30;
  double mu=1.32712e11;
  double radius=695700.0;
};

class earth {
 private:

 public:
  std::string name="EARTH";
  double mass=5.972e24;
  double mu=5.972e24*G;
  double radius=6378.0;
  double J2=1.082635854e-3;
  //  'zs':atm[:,0], #km
  // 'rhos': atm[:,1]*10**8, # kg/km³
  // 'atm_rot_vector':np.array([0.0,0.0,72.9211e-6]), #rad/s
  double deorbit_altitude=100;
  //'spice_file':A+'de432s.bsp'
};

class perturbations {
 private:

 public:
  bool J2=false;
  bool aero=false;
  bool srp=false;
  //n_bodies=[];
  double Cd=0;
  double A=0;
  double trhust=0;
  double thrust_direction=0;
  double isp=0;
  double CR=0;
  double A_srp=0;
  double C20=0;
  double mu=0;
};

