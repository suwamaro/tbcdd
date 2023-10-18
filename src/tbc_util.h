/***************************************************************************
* Functions used in the tight-binding calculation.
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _TBC_UTIL_H_
#define _TBC_UTIL_H_

#include <cmath>

namespace TBC
{
  double round0(double v){
    if (std::abs(v) < 1e-12) { return 0.; }
    else { return v; }
  }
  
  double fermi_energy(double x, double kT, double mu) {
    double alpha = (x-mu)/std::abs(kT);
    if (kT < 1e-15 || std::abs(alpha) > 20) {
      return (x < mu) ? (x-mu) : 0.0;
    }
    else {
      if ( alpha >= 0 ) {
	return - kT * log(1. + exp(-alpha));
      } else {
	return x - mu - kT * log(1. + exp(alpha));
      }
    }
  }
  
  double fermi_density(double x, double kT, double mu) {
    double alpha = (x-mu)/std::abs(kT);
    if (kT < 1e-15 || std::abs(alpha) > 20) {
      return (x < mu) ? 1.0 : 0.0;
    } else {
      if ( alpha >= 0 ) {
	return exp(-alpha) / (exp(-alpha) + 1.);
      } else {
	return 1. / (exp(alpha) + 1.);
      }
    }
  }
  
}

#endif // _TBC_UTIL_H_
