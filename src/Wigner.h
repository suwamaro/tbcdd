/***************************************************************************
* Wigner's matrix
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/
#ifndef _WIGNER_H_
#define _WIGNER_H_

#include "tbc.h"

using namespace std::complex_literals;

namespace TBC
{
  // Wigner's D matrix for J=1/2  
  class Wigner05 {
  public:
    Wigner05(){};
    ~Wigner05(){}
    cx_mat_type Dmat(double alpha, double beta, double gamma) const;
    
  private:
    cx_mat_type diag(double alpha) const;
    mat_type dmat(double beta) const;    
  };

  cx_mat_type Wigner05::Dmat(double alpha, double beta, double gamma) const {
    // Rotation operator in the z-y-z convention    
    return diag(alpha) * dmat(beta) * diag(gamma);
  }

  cx_mat_type Wigner05::diag(double alpha) const {
    cx_double z = std::exp(0.5 * 1i * alpha);
    cx_vec_type v({1./z, z});
    return arma::diagmat(v);
  }
  
  mat_type Wigner05::dmat(double beta) const {
    mat_type d(2,2, arma::fill::zeros);
    d(0,0) = cos(0.5*beta);
    d(0,1) = - sin(0.5*beta);
    d(1,0) = sin(0.5*beta);
    d(1,1) = cos(0.5*beta);
    return d;
  }

  
  // Wigner's D matrix for J=2
  class Wigner2 {
  public:
    Wigner2();
    ~Wigner2(){}
    cx_mat_type Dmat(double alpha, double beta, double gamma) const;
    
  private:
    cx_mat_type diag(double alpha) const;
    mat_type dmat(double beta) const;    
    mat_type JpJm;

  };

  Wigner2::Wigner2():JpJm(5,5, arma::fill::zeros){
    double a = 2.;
    double b = sqrt(6.);
  
    JpJm(0,1) = JpJm(3,4) = a;
    JpJm(1,0) = JpJm(4,3) = - a;
    JpJm(1,2) = JpJm(2,3) = b;
    JpJm(2,1) = JpJm(3,2) = - b;
  }

  cx_mat_type Wigner2::Dmat(double alpha, double beta, double gamma) const {
    // Rotation operator in the z-y-z convention    
    return diag(alpha) * dmat(beta) * diag(gamma);
  }

  cx_mat_type Wigner2::diag(double alpha) const {
    cx_double z = std::exp(1i * alpha);
    cx_double z2 = z * z;
    cx_vec_type v({1./z2, 1./z, 1.0, z, z2});
    return arma::diagmat(v);
  }
  
  mat_type Wigner2::dmat(double beta) const {
    return arma::expmat(-0.5*beta*JpJm);
  }  
  
}

#endif // _WIGNER_H_
