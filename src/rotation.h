/***************************************************************************
* Rotation matrices
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _ROTATION_H_
#define _ROTATION_H_

#include "Wigner.h"
#include "slater_koster.h"
#include "orbitals.h"
#include "spin_orbit.h"
#include "parameters_tbc.h"

namespace TBC
{
  class Rotations {
  public:
    Rotations();
    explicit Rotations(Parameters const& pr);
    ~Rotations(){}
    std::tuple<cx_mat_type, cx_mat_type> rotation_and_tilt(double theta, double phi) const;
    mat_type Rz(double theta) const;
    mat_type Ry(double theta) const;    
    mat_type O3_rotation(double alpha, double beta, double gamma) const;
    cx_mat_type spin_rotation() const;
    mat_type O3_spin_rotation(int i) const;    
    cx_mat_type J_basis() const;    
    cx_mat_type elems(int i, int j, cx_mat_type const& mat) const;
    std::tuple<int, int> get_index_pair(int i) const;
    
    int n_sublattices;
    std::vector<double> thetas, phis;
    Wigner05 wigner05;    
    Wigner2 wigner2;
    SlaterKoster sk;
    dOrbital dorb;
    SpinOrbit so;
    std::vector<cx_mat_type> R_locs;
    std::vector<cx_mat_type> spin_locs;
    cx_mat_type U_spin, U_J;
    cx_mat_type R_x, R_y, R_xy, R_xyb;
    cx_mat_type R_z, R_xz, R_xbz, R_yz, R_ybz;
    cx_mat_type R_xyz, R_xbyz, R_xybz, R_xbybz;

    /* Assume the unite cell of the two-layer structure */
    /* Layer1: 0  1  */
    /* Layer2: 3  2  */    
    cx_mat_type t_x01, t_x10, t_y01, t_y10, t_xy0, t_xyb0, t_xy1, t_xyb1, t_x0, t_x1, t_y0, t_y1;  // Layer 1
    cx_mat_type t_x23, t_x32, t_y23, t_y32, t_xy2, t_xyb2, t_xy3, t_xyb3, t_x2, t_x3, t_y2, t_y3;  // Layer 2
    cx_mat_type t_z03, t_z12, t_xz02, t_xbz02, t_xz13, t_xbz13, t_yz02, t_ybz02, t_yz13, t_ybz13;  // Interlayer
    cx_mat_type t_xyz03, t_xbyz03, t_xybz03, t_xbybz03, t_xyz12, t_xbyz12, t_xybz12, t_xbybz12;        // Interlayer
  };

  std::tuple<cx_mat_type, cx_mat_type> Rotations::rotation_and_tilt(double theta, double phi) const {
    /* Converting from angle to radian */
    theta *= M_PI / 180.;
    phi *= M_PI / 180.;
    cx_mat_type D05 = wigner05.Dmat(- M_PI/4. - theta, - phi, M_PI/4.);  // Inverse transformation
    cx_mat_type D2 = wigner2.Dmat(- M_PI/4. - theta, - phi, M_PI/4.);  // Inverse transformation
    return std::make_tuple(D05, D2);
  }

  mat_type Rotations::Rz(double theta) const {
    double ct = cos(theta);
    double st = sin(theta);
    mat_type R = {{ct, - st, 0},{st, ct, 0},{0, 0, 1.}};
    return R;
  }

  mat_type Rotations::Ry(double theta) const {
    double ct = cos(theta);
    double st = sin(theta);
    mat_type R = {{ct, 0, st},{0., 1., 0},{- st, 0, ct}};
    return R;
  }  
  
  mat_type Rotations::O3_rotation(double alpha, double beta, double gamma) const {
    return Rz(alpha) * Ry(beta) * Rz(gamma);
  }  
  
  cx_mat_type Rotations::elems(int i, int j, cx_mat_type const& mat) const {
    mat_type pos(2*n_sublattices, 2*n_sublattices, arma::fill::zeros);      
    pos(i,j) = 1;
    return arma::kron(pos, mat);
  }

  std::tuple<int, int> Rotations::get_index_pair(int i) const {
    // Indices for spin up and spin down
    int i1 = ((i >> 1) << 2) | (i & 1);
    int i2 = i1 + 2;
    return std::make_tuple(i1, i2);
  }
  
  cx_mat_type Rotations::spin_rotation() const{
    /* For t2g */
    int mat_size = 2 * n_sublattices * dorb.psi().n_cols;
    cx_mat_type U(mat_size, mat_size, arma::fill::zeros);
    
    for(int idx=0; idx < n_sublattices; ++idx){
      cx_double z = spin_locs[idx](0,0);
      cx_double w = spin_locs[idx](1,0);
      cx_mat_type U_diag = arma::diagmat(cx_vec_type{std::conj(z), std::conj(z), z});
      cx_mat_type U_off = arma::diagmat(cx_vec_type{w, w, - std::conj(w)});

      int idx1, idx2;
      std::tie(idx1, idx2) = get_index_pair(idx);
      U += elems(idx1, idx1, U_diag);
      U += elems(idx1, idx2, U_off);
      U += elems(idx2, idx1, - U_off.t());
      U += elems(idx2, idx2, U_diag.t());    
    }

    // // for check
    // U = arma::eye<arma::cx_mat>(24,24);
    
    return U;
  }

  mat_type Rotations::O3_spin_rotation(int i) const {
    /* For O(3) classical spin in sublattice i */
    double theta = thetas[i];
    double phi = phis[i];
    return O3_rotation(- M_PI/4., phi, M_PI/4. + theta);
  }  

  cx_mat_type Rotations::J_basis() const {
    int mat_size = 2 * n_sublattices * so.Vdown.n_cols;
    cx_mat_type U(mat_size, mat_size, arma::fill::zeros);
    
    for(int idx=0; idx < n_sublattices; ++idx){
      int idx1, idx2;
      std::tie(idx1, idx2) = get_index_pair(idx);
      U += elems(idx1, idx1, so.Vdown.t());
      U += elems(idx2, idx2, so.Vup.t());
    }
    
    return U;
  }  
  
  Rotations::Rotations(Parameters const& pr):wigner05(),wigner2(), sk(), dorb(), so(pr){
    n_sublattices = pr.n_sublattices;
    thetas = pr.thetas;
    phis = pr.phis;
    R_locs.clear();
    R_locs.resize(n_sublattices);
    spin_locs.resize(n_sublattices);
    for(int i=0; i < n_sublattices; ++i){
      std::tie(spin_locs[i], R_locs[i]) = rotation_and_tilt(thetas[i], phis[i]);
    }
    
    R_x = wigner2.Dmat(-pr.alpha_mono, M_PI/2, 0);  // z -> x
    R_y = wigner2.Dmat(M_PI/2 + pr.alpha_mono, M_PI/2, - M_PI/2);  // z -> y: The last -pi/2 does not matter.
    R_xy = wigner2.Dmat(M_PI/4, M_PI/2, - M_PI/4);  // z -> x y
    R_xyb = wigner2.Dmat(- M_PI/4, M_PI/2, M_PI/4);  // z -> x -y
    
    R_z = arma::eye<arma::cx_mat>(5,5);
    R_xz = wigner2.Dmat(- pr.alpha_mono, M_PI/4, 0);
    R_xbz = wigner2.Dmat(- pr.alpha_mono, - M_PI/4, 0);
    R_yz = wigner2.Dmat(M_PI/2 + pr.alpha_mono, M_PI/4, - M_PI/2);
    R_ybz = wigner2.Dmat(- M_PI/2 + pr.alpha_mono, M_PI/4, M_PI/2);

    double angle3 = atan(sqrt(2.));
    R_xyz = wigner2.Dmat(M_PI/4, angle3, - M_PI/4);  // z -> x y z
    R_xbyz = wigner2.Dmat(- M_PI/4, - angle3, M_PI/4);  // z -> -x y z
    R_xybz = wigner2.Dmat(- M_PI/4, angle3, M_PI/4);  // z -> x -y z
    R_xbybz = wigner2.Dmat(M_PI/4, - angle3, - M_PI/4);  // z -> x y z    
    
    // Hopping amplitudes
    t_x01 = dorb.psi().t() * R_locs[0] * R_x * sk.t_SK() * R_x.t() * R_locs[1].t() * dorb.psi();
    t_x10 = dorb.psi().t() * R_locs[1] * R_x * sk.t_SK() * R_x.t() * R_locs[0].t() * dorb.psi();
    t_y01 = dorb.psi().t() * R_locs[0] * R_y * sk.t_SK() * R_y.t() * R_locs[1].t() * dorb.psi();
    t_y10 = dorb.psi().t() * R_locs[1] * R_y * sk.t_SK() * R_y.t() * R_locs[0].t() * dorb.psi();    
    t_xy0 = dorb.psi().t() * R_locs[0] * R_xy * sk.t_SK() * R_xy.t() * R_locs[0].t() * dorb.psi();
    t_xyb0 = dorb.psi().t() * R_locs[0] * R_xyb * sk.t_SK() * R_xyb.t() * R_locs[0].t() * dorb.psi();
    t_xy1 = dorb.psi().t() * R_locs[1] * R_xy * sk.t_SK() * R_xy.t() * R_locs[1].t() * dorb.psi();
    t_xyb1 = dorb.psi().t() * R_locs[1] * R_xyb * sk.t_SK() * R_xyb.t() * R_locs[1].t() * dorb.psi();
    t_x0 = dorb.psi().t() * R_locs[0] * R_x * sk.t_SK() * R_x.t() * R_locs[0].t() * dorb.psi();
    t_x1 = dorb.psi().t() * R_locs[1] * R_x * sk.t_SK() * R_x.t() * R_locs[1].t() * dorb.psi();
    t_y0 = dorb.psi().t() * R_locs[0] * R_y * sk.t_SK() * R_y.t() * R_locs[0].t() * dorb.psi();
    t_y1 = dorb.psi().t() * R_locs[1] * R_y * sk.t_SK() * R_y.t() * R_locs[1].t() * dorb.psi();
    
    t_x23 = dorb.psi().t() * R_locs[2] * R_x * sk.t_SK() * R_x.t() * R_locs[3].t() * dorb.psi();
    t_x32 = dorb.psi().t() * R_locs[3] * R_x * sk.t_SK() * R_x.t() * R_locs[2].t() * dorb.psi();
    t_y23 = dorb.psi().t() * R_locs[2] * R_y * sk.t_SK() * R_y.t() * R_locs[3].t() * dorb.psi();
    t_y32 = dorb.psi().t() * R_locs[3] * R_y * sk.t_SK() * R_y.t() * R_locs[2].t() * dorb.psi();    
    t_xy2 = dorb.psi().t() * R_locs[2] * R_xy * sk.t_SK() * R_xy.t() * R_locs[2].t() * dorb.psi();
    t_xyb2 = dorb.psi().t() * R_locs[2] * R_xyb * sk.t_SK() * R_xyb.t() * R_locs[2].t() * dorb.psi();
    t_xy3 = dorb.psi().t() * R_locs[3] * R_xy * sk.t_SK() * R_xy.t() * R_locs[3].t() * dorb.psi();
    t_xyb3 = dorb.psi().t() * R_locs[3] * R_xyb * sk.t_SK() * R_xyb.t() * R_locs[3].t() * dorb.psi();
    t_x2 = dorb.psi().t() * R_locs[2] * R_x * sk.t_SK() * R_x.t() * R_locs[2].t() * dorb.psi();
    t_x3 = dorb.psi().t() * R_locs[3] * R_x * sk.t_SK() * R_x.t() * R_locs[3].t() * dorb.psi();
    t_y2 = dorb.psi().t() * R_locs[2] * R_y * sk.t_SK() * R_y.t() * R_locs[2].t() * dorb.psi();
    t_y3 = dorb.psi().t() * R_locs[3] * R_y * sk.t_SK() * R_y.t() * R_locs[3].t() * dorb.psi();    
    
    t_z03 = dorb.psi().t() * R_locs[0] * R_z * sk.t_SK() * R_z.t() * R_locs[3].t() * dorb.psi();
    t_z12 = dorb.psi().t() * R_locs[1] * R_z * sk.t_SK() * R_z.t() * R_locs[2].t() * dorb.psi();
    t_xz02 = dorb.psi().t() * R_locs[0] * R_xz * sk.t_SK() * R_xz.t() * R_locs[2].t() * dorb.psi();
    t_xbz02 = dorb.psi().t() * R_locs[0] * R_xbz * sk.t_SK() * R_xbz.t() * R_locs[2].t() * dorb.psi();
    t_xz13 = dorb.psi().t() * R_locs[1] * R_xz * sk.t_SK() * R_xz.t() * R_locs[3].t() * dorb.psi();
    t_xbz13 = dorb.psi().t() * R_locs[1] * R_xbz * sk.t_SK() * R_xbz.t() * R_locs[3].t() * dorb.psi();
    t_yz02 = dorb.psi().t() * R_locs[0] * R_yz * sk.t_SK() * R_yz.t() * R_locs[2].t() * dorb.psi();
    t_ybz02 = dorb.psi().t() * R_locs[0] * R_ybz * sk.t_SK() * R_ybz.t() * R_locs[2].t() * dorb.psi();
    t_yz13 = dorb.psi().t() * R_locs[1] * R_yz * sk.t_SK() * R_yz.t() * R_locs[3].t() * dorb.psi();
    t_ybz13 = dorb.psi().t() * R_locs[1] * R_ybz * sk.t_SK() * R_ybz.t() * R_locs[3].t() * dorb.psi();

    t_xyz03 = dorb.psi().t() * R_locs[0] * R_xyz * sk.t_SK() * R_xyz.t() * R_locs[3].t() * dorb.psi();
    t_xbyz03 = dorb.psi().t() * R_locs[0] * R_xbyz * sk.t_SK() * R_xbyz.t() * R_locs[3].t() * dorb.psi();
    t_xybz03 = dorb.psi().t() * R_locs[0] * R_xybz * sk.t_SK() * R_xybz.t() * R_locs[3].t() * dorb.psi();
    t_xbybz03 = dorb.psi().t() * R_locs[0] * R_xbybz * sk.t_SK() * R_xbybz.t() * R_locs[3].t() * dorb.psi();
    t_xyz12 = dorb.psi().t() * R_locs[1] * R_xyz * sk.t_SK() * R_xyz.t() * R_locs[2].t() * dorb.psi();
    t_xbyz12 = dorb.psi().t() * R_locs[1] * R_xbyz * sk.t_SK() * R_xbyz.t() * R_locs[2].t() * dorb.psi();
    t_xybz12 = dorb.psi().t() * R_locs[1] * R_xybz * sk.t_SK() * R_xybz.t() * R_locs[2].t() * dorb.psi();
    t_xbybz12 = dorb.psi().t() * R_locs[1] * R_xbybz * sk.t_SK() * R_xbybz.t() * R_locs[2].t() * dorb.psi();    

    // Spin rotation
    U_spin = spin_rotation();

    // J basis
    U_J = J_basis();
  }
}

#endif // _ROTATION_H_
