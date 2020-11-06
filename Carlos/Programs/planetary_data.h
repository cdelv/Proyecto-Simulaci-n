//planetary_data.h
#pragma once

const double G_meters=6.67408e-11;
const double G=G_meters*10e-9;

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
};
