/***************************************************************************
* Electron d-orbitals
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _ORBITALS_H_
#define _ORBITALS_H_

#include "tbc.h"

using namespace std::complex_literals;

namespace TBC
{
  class dOrbital {
  public:
    dOrbital();
    ~dOrbital(){}
    cx_mat_type psi() const { return psi_; }

  private:
    cx_mat_type psi_;
  };

  dOrbital::dOrbital(){
    double a = 1. / sqrt(2.);
    cx_double ai = 1i * a;    
    psi_ = {{0,0,-ai},{ai,-a,0},{0,0,0},{ai,a,0},{0,0,ai}};
  }

}

#endif // _ORBITALS_H_
