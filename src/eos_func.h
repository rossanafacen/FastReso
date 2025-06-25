/*
 * Copyright (c) 2018-2021 Aleksas Mazeliauskas, Stefan Floerchinger, 
 *                    Eduardo Grossi, and Derek Teaney
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */
#ifndef FASTRESO_eos_func_h
#define FASTRESO_eos_func_h
#include "gsl/gsl_sf_bessel.h"
double pfunc(double mOverT, double QmuOverT, int sign=1);
double efunc(double mOverT, double QmuOverT, int sign=1);
double sfunc(double mOverT, double QmuOverT, int sign=1);
double nfunc(double mOverT, double QmuOverT, int sign=1);

namespace eos_params{
 const int nmax = 105;
 const double rel = 1e-10;

}




// return pressure/T^4/nu
double pfunc(double mOverT, double QmuOverT, int sign)
{
  double x = mOverT;
  double pressure = 0;
  int A =sign;
  bool donext = true;
  using namespace eos_params;
  int i=1;
  //for (int i =1; i <= nmax; i++){
  while (donext and i<nmax ){
    A *=sign; 
    double dp = A/pow(i,4)*( x>0 ? gsl_sf_bessel_Kn_scaled(2,x*i)*(x*i)*(x*i)/2.0: 1)*exp(-i*(x-QmuOverT));
    pressure += dp;
    i++;
    if ( abs(dp/pressure) < rel){
      donext=false;
    }
  }
  return pressure/M_PI/M_PI;
}
// return energy/T^4/nu
double efunc(double mOverT, double QmuOverT, int sign){
  double x = mOverT;
  double energy = 0;
  int A =sign;
  bool donext = true;
  using namespace eos_params;
  int i=1;
  while (donext and i<nmax ){
    A *=sign; 
    double de =   A/pow(i,4)*(x>0 ? (gsl_sf_bessel_Kn_scaled(2,x*i) + 1.0/3.0*(x*i)*gsl_sf_bessel_Kn_scaled(1,x*i)   )*(x*i)*(x*i)/2.0 : 1 )*exp(-i*(x-QmuOverT));
    energy += de;
    i++;
    if ( abs(de/energy) < rel){
      donext=false;
    }
  }
  return 3*energy/M_PI/M_PI;
}



// return entropy/T^3
double sfunc(double mOverT, double QmuOverT, int sign){
  double x = mOverT;
  double entropy = 0;
  int A =sign;


  bool donext = true;
  using namespace eos_params;
  int i=1;
  while (donext and i<nmax ){
    A *=sign; 
    double dent = A/pow(i,4)*( x >0 ? ( (4 - QmuOverT*i)*gsl_sf_bessel_Kn_scaled(2,x*i) + (x*i)*gsl_sf_bessel_Kn_scaled(1,x*i)   )*(x*i)*(x*i)/2.0 : 4 - QmuOverT*i )*exp(-i*(x-QmuOverT));
    entropy += dent;
    i++;
    if ( abs(dent/entropy) < rel){
      donext=false;
    }
  }
  return entropy/M_PI/M_PI;
}
// return number density n/T^3
double nfunc(double mOverT, double QmuOverT, int sign){
  double x = mOverT;
  double number = 0;
  int A =sign;
    
  bool donext = true;
  using namespace eos_params;
  int i=1;
  while (donext and i<nmax ){
    A *=sign; 
    double dn = A/pow(i,3)*( x>0 ? (gsl_sf_bessel_Kn_scaled(2,x*i) )*(x*i)*(x*i)/2. : 1)*exp(-i*(x-QmuOverT));
    number += dn;
    i++;
    if ( abs(dn/number) < rel){
      donext=false;
    }

  }
  return number/M_PI/M_PI;
}
#endif
