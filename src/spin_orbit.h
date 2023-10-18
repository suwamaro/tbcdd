/***************************************************************************
* Spin-orbit interaction
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _SPIN_ORBIT_H_
#define _SPIN_ORBIT_H_

#include "tbc.h"

using namespace std::complex_literals;

namespace TBC
{
  class SpinOrbit {
  public:
    SpinOrbit();
    explicit SpinOrbit(Parameters const& pr);
    ~SpinOrbit(){}
    cx_mat_type H_SO(double lambda, double sigma) const;
    double E_J05() const;
    void output_eigenvalue(ofstream& out) const;
    cx_mat_type H_SO_up, H_SO_down;
    mat_type H_TD;
    cx_mat_type H_loc_up, H_loc_down;
    vec_type Edown, Eup;
    cx_mat_type Vdown, Vup;
  };
  
  SpinOrbit::SpinOrbit(Parameters const& pr){
    H_SO_up = H_SO(pr.lambda, 1.0);
    H_SO_down = H_SO(pr.lambda, -1.0);
    H_TD = {{pr.Delta, 0, 0},{0, pr.Delta, 0}, {0,0,0}};
    H_loc_up = H_SO_up  + H_TD;
    H_loc_down = H_SO_down + H_TD;

    // Diagonalization
    arma::eig_sym(Edown, Vdown, H_loc_down);
    arma::eig_sym(Eup, Vup, H_loc_up);
  }

  cx_mat_type SpinOrbit::H_SO(double lambda, double sigma) const {
    cx_double b = lambda/2.;  
    cx_mat_type H = {{0, 1i * sigma * b, - sigma * b}, {- 1i * sigma * b, 0, 1i * b}, {- sigma * b, - 1i * b, 0}};
    return H;
  }

  double SpinOrbit::E_J05() const {
    return Edown[2];
  }
  
  void SpinOrbit::output_eigenvalue(ofstream& out) const {
    out << "Edown" << std::endl;
    out << Edown << std::endl;
    out << "Vdown" << std::endl;
    out << Vdown << std::endl;

    out << "Eup" << std::endl;
    out << Eup << std::endl;
    out << "Vup" << std::endl;
    out << Vup << std::endl;    
  }
}

#endif // _SPIN_ORBIT_H_
