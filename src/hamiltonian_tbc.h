/***************************************************************************
* Tight-binding hamiltonian
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _HAMILTONIAN_TBC_H_
#define _HAMILTONIAN_TBC_H_

#include "Wigner.h"
#include "slater_koster.h"
#include "parameters_tbc.h"
#include "rotation.h"
#include "mean_field.h"
#include "tbc_util.h"

namespace TBC
{
  class Hoppings {
  public:
    Hoppings(){};
    ~Hoppings(){}
    void output(ofstream &out) const;
    double tx, tx_bar, ty, ty_bar, txy, txyb, tz, tz_bar;
    double tx2_0, tx2_1, ty2_0, ty2_1;  // Depend on sublattices.
    double tfx, tfx_bar, tfy, tfy_bar, tfz, tfz_bar;
  };

  /* Member functions of Hoppings */
  void Hoppings::output(ofstream &out) const {
    out << "tx = " << tx << std::endl;
    out << "tx_bar = " << tx_bar << std::endl;
    out << "ty = " << ty << std::endl;
    out << "ty_bar = " << ty_bar << std::endl;
    out << "txy = " << txy << std::endl;
    out << "txyb = " << txyb << std::endl;
    out << "tx2_0 = " << tx2_0 << std::endl;
    out << "tx2_1 = " << tx2_1 << std::endl;    
    out << "ty2_0 = " << ty2_0 << std::endl;
    out << "ty2_1 = " << ty2_1 << std::endl;    
    
    out << "tz = " << tz << std::endl;
    out << "tz_bar = " << tz_bar << std::endl;

    out << "tfx = " << tfx << std::endl;
    out << "tfx_bar = " << tfx_bar << std::endl;
    out << "tfy = " << tfy << std::endl;
    out << "tfy_bar = " << tfy_bar << std::endl;
    out << "tfz = " << tfz << std::endl;
    out << "tfz_bar = " << tfz_bar << std::endl;    
  }
  
  class Hamiltonian {
  public:
    Hamiltonian();
    explicit Hamiltonian(Parameters const& pr);
    std::tuple<cx_mat_type, cx_mat_type> separate_mat(cx_mat_type const& mat) const;
    cx_mat_type H_elem(int i, int j, cx_mat_type const& mat) const;
    double hopping_amplitude(std::string const& bond) const;
    cx_mat_type H_k(double kx, double ky, double) const;
    cx_mat_type H_k_J05(double kx, double ky, double) const;
    cx_mat_type extract_J05(cx_mat_type const& h) const;
    void get_hopping_parameters();    
    void output_hopping_parameters(ofstream &out) const;
    void output_magnetization(ofstream &out) const;
    void output_magnetization_in_basis(ofstream &out, mat_type const& R) const;
    cx_mat_type onsite_U_mf(int i) const;    
    void set_U_mf();
    void set_U(double U);
    void set_mf(MeanField const& mf);
    void print_H_J05_hopping() const;
    int n_sublattices() const { return rot.n_sublattices; }
    ~Hamiltonian(){}
    
    int n;  // Matrix dimension
    bool single_layer, bilayer;
    double t;
    double U;
    double eta2, eta_a, eta3, etaz, etaz2, etaz3;  // Distance factors
    Rotations rot;
    cx_mat_type H_x, H_y, H_xy, H_xyb, H_x2, H_y2, H_z, H_xz, H_xbz, H_yz, H_ybz, H_xyz, H_xbyz, H_xybz, H_xbybz;  // Electron hopping
    cx_mat_type H_x_J05, H_y_J05, H_xy_J05, H_xyb_J05, H_x2_J05, H_y2_J05, H_z_J05, H_xz_J05, H_xbz_J05, H_yz_J05, H_ybz_J05, H_xyz_J05, H_xbyz_J05, H_xybz_J05, H_xbybz_J05;  // Electron hopping between Jeff=1/2 orbitals    
    cx_mat_type H_local, H_local_J05;  // Local part
    arma::uvec index_J05;
    Hoppings hop;
    MeanField mf;
    cx_mat_type Umf;
  };

  /* Member functions of Hamiltonian */
  Hamiltonian::Hamiltonian(Parameters const& pr)
    :n(24),rot(pr),
     H_x(n,n),H_y(n,n),H_xy(n,n),H_xyb(n,n),H_x2(n,n),H_y2(n,n),H_z(n,n),H_xz(n,n),H_xbz(n,n),H_yz(n,n),H_ybz(n,n),H_xyz(n,n),H_xbyz(n,n),H_xybz(n,n),H_xbybz(n,n),
     H_x_J05(n/3,n/3),H_y_J05(n/3,n/3),H_xy_J05(n/3,n/3),H_xyb_J05(n/3,n/3),H_x2_J05(n/3,n/3),H_y2_J05(n/3,n/3),H_z_J05(n/3,n/3),H_xz_J05(n/3,n/3),H_xbz_J05(n/3,n/3),H_yz_J05(n/3,n/3),H_ybz_J05(n/3,n/3),H_xyz_J05(n/3,n/3),H_xbyz_J05(n/3,n/3),H_xybz_J05(n/3,n/3),H_xbybz_J05(n/3,n/3),
     H_local(n,n),H_local_J05(n/3,n/3),Umf(n/3,n/3){
    /* Parameters */
    single_layer = pr.single_layer;    
    bilayer = pr.bilayer;
    t = pr.t;
    U = pr.U;
    eta2 = pr.eta2;
    eta_a = pr.eta_a;
    eta3 = pr.eta3;
    etaz = pr.etaz;
    etaz2 = pr.etaz2;
    etaz3 = pr.etaz3;    
    if (etaz2 != 0 || etaz3 != 0) {
      std::cout << "Warning: Contributions from the interlayer hoppings are mixed in the four-sublattice formulation. Interlayer hopping amplitudes for shifting 2D positions are not available." << std::endl;
    }
    
    /* Indices of Jeff=1/2 orbitals, ad-hoc. */
    index_J05 = {2,5,8,11,14,17,20,23};
    
    /* Setting mean fields */
    mf.init(n_sublattices());
    double ne = 2. * pr.filling;
    mf.ns = {ne, ne, ne, ne};
    mf.mags[0] = pr.m1;
    mf.mags[1] = pr.m2;
    mf.mags[2] = pr.m3;
    mf.mags[3] = pr.m4;    
    set_U_mf();
    
    // Matrix dimension is 4*3*2=24.
    // (0,yz,down) (0,xz,down) (0,xy,up) (1,yz,down) (1,xz,down) (1,xy,up)
    // (0,yz,up) (0,xz,up) (0,xy,down) (1,yz,up) (1,xz,up) (1,xy,down)
    // (2,yz,down) (2,xz,down) (2,xy,up) (3,yz,down) (3,xz,down) (3,xy,up)
    // (2,yz,up) (2,xz,up) (2,xy,down) (3,yz,up) (3,xz,up) (3,xy,down)
    
    // x
    H_x.zeros();
    H_x += H_elem(0, 1, rot.t_x01);
    H_x += H_elem(1, 0, rot.t_x10);    
    H_x += H_elem(2, 3, rot.t_x23);
    H_x += H_elem(3, 2, rot.t_x32);    
    H_x = rot.U_spin * H_x * rot.U_spin.t();  // Spin rotation
    H_x_J05 = extract_J05(H_x);
    
    // y
    H_y.zeros();
    H_y += H_elem(0, 1, rot.t_y01);
    H_y += H_elem(1, 0, rot.t_y10);    
    H_y += H_elem(2, 3, rot.t_y23);
    H_y += H_elem(3, 2, rot.t_y32);    
    H_y = rot.U_spin * H_y * rot.U_spin.t();  // Spin rotation
    H_y_J05 = extract_J05(H_y);
    
    // xy
    H_xy.zeros();
    H_xy += H_elem(0, 0, rot.t_xy0);
    H_xy += H_elem(1, 1, rot.t_xy1);
    H_xy += H_elem(2, 2, rot.t_xy2);
    H_xy += H_elem(3, 3, rot.t_xy3);
    H_xy = rot.U_spin * H_xy * rot.U_spin.t();  // Spin rotation
    H_xy_J05 = extract_J05(H_xy);
    
    // xyb
    H_xyb.zeros();
    H_xyb += H_elem(0, 0, rot.t_xyb0);
    H_xyb += H_elem(1, 1, rot.t_xyb1);
    H_xyb += H_elem(2, 2, rot.t_xyb2);
    H_xyb += H_elem(3, 3, rot.t_xyb3);    
    H_xyb = rot.U_spin * H_xyb * rot.U_spin.t();  // Spin rotation
    H_xyb_J05 = extract_J05(H_xyb);
    
    // x2
    H_x2.zeros();
    H_x2 += H_elem(0, 0, rot.t_x0);
    H_x2 += H_elem(1, 1, rot.t_x1);
    H_x2 += H_elem(2, 2, rot.t_x2);
    H_x2 += H_elem(3, 3, rot.t_x3);    
    H_x2 = rot.U_spin * H_x2 * rot.U_spin.t();  // Spin rotation
    H_x2_J05 = extract_J05(H_x2);
    
    // y2
    H_y2.zeros();
    H_y2 += H_elem(0, 0, rot.t_y0);
    H_y2 += H_elem(1, 1, rot.t_y1);
    H_y2 += H_elem(2, 2, rot.t_y2);
    H_y2 += H_elem(3, 3, rot.t_y3);    
    H_y2 = rot.U_spin * H_y2 * rot.U_spin.t();  // Spin rotation    
    H_y2_J05 = extract_J05(H_y2);
    
    // z
    H_z.zeros();
    H_z += H_elem(0, 3, rot.t_z03);
    H_z += H_elem(3, 0, rot.t_z03.t());    
    H_z += H_elem(1, 2, rot.t_z12);
    H_z += H_elem(2, 1, rot.t_z12.t());    
    H_z = rot.U_spin * H_z * rot.U_spin.t();  // Spin rotation
    H_z_J05 = extract_J05(H_z);
    
    // xz
    H_xz.zeros();
    H_xz += H_elem(0, 2, rot.t_xz02);
    H_xz += H_elem(2, 0, rot.t_xz02.t());
    H_xz += H_elem(1, 3, rot.t_xz13);
    H_xz += H_elem(3, 1, rot.t_xz13.t());
    H_xz = rot.U_spin * H_xz * rot.U_spin.t();  // Spin rotation
    H_xz_J05 = extract_J05(H_xz);
    
    // xbz
    H_xbz.zeros();
    H_xbz += H_elem(0, 2, rot.t_xbz02);
    H_xbz += H_elem(2, 0, rot.t_xbz02.t());
    H_xbz += H_elem(1, 3, rot.t_xbz13);
    H_xbz += H_elem(3, 1, rot.t_xbz13.t());
    H_xbz = rot.U_spin * H_xbz * rot.U_spin.t();  // Spin rotation
    H_xbz_J05 = extract_J05(H_xbz);
    
    // yz
    H_yz.zeros();
    H_yz += H_elem(0, 2, rot.t_yz02);
    H_yz += H_elem(2, 0, rot.t_yz02.t());
    H_yz += H_elem(1, 3, rot.t_yz13);
    H_yz += H_elem(3, 1, rot.t_yz13.t());
    H_yz = rot.U_spin * H_yz * rot.U_spin.t();  // Spin rotation
    H_yz_J05 = extract_J05(H_yz);
    
    // ybz
    H_ybz.zeros();
    H_ybz += H_elem(0, 2, rot.t_ybz02);
    H_ybz += H_elem(2, 0, rot.t_ybz02.t());
    H_ybz += H_elem(1, 3, rot.t_ybz13);
    H_ybz += H_elem(3, 1, rot.t_ybz13.t());
    H_ybz = rot.U_spin * H_ybz * rot.U_spin.t();  // Spin rotation
    H_ybz_J05 = extract_J05(H_ybz);
    
    // xyz
    H_xyz.zeros();
    H_xyz += H_elem(0, 3, rot.t_xyz03);
    H_xyz += H_elem(3, 0, rot.t_xyz03.t());    
    H_xyz += H_elem(1, 2, rot.t_xyz12);
    H_xyz += H_elem(2, 1, rot.t_xyz12.t());    
    H_xyz = rot.U_spin * H_xyz * rot.U_spin.t();  // Spin rotation
    H_xyz_J05 = extract_J05(H_xyz);
    
    // xbyz
    H_xbyz.zeros();
    H_xbyz += H_elem(0, 3, rot.t_xbyz03);
    H_xbyz += H_elem(3, 0, rot.t_xbyz03.t());    
    H_xbyz += H_elem(1, 2, rot.t_xbyz12);
    H_xbyz += H_elem(2, 1, rot.t_xbyz12.t());    
    H_xbyz = rot.U_spin * H_xbyz * rot.U_spin.t();  // Spin rotation
    H_xbyz_J05 = extract_J05(H_xbyz);
    
    // xybz
    H_xybz.zeros();
    H_xybz += H_elem(0, 3, rot.t_xybz03);
    H_xybz += H_elem(3, 0, rot.t_xybz03.t());    
    H_xybz += H_elem(1, 2, rot.t_xybz12);
    H_xybz += H_elem(2, 1, rot.t_xybz12.t());    
    H_xybz = rot.U_spin * H_xybz * rot.U_spin.t();  // Spin rotation
    H_xybz_J05 = extract_J05(H_xybz);
    
    // xbybz
    H_xbybz.zeros();
    H_xbybz += H_elem(0, 3, rot.t_xbybz03);
    H_xbybz += H_elem(3, 0, rot.t_xbybz03.t());    
    H_xbybz += H_elem(1, 2, rot.t_xbybz12);
    H_xbybz += H_elem(2, 1, rot.t_xbybz12.t());    
    H_xbybz = rot.U_spin * H_xbybz * rot.U_spin.t();  // Spin rotation            
    H_xbybz_J05 = extract_J05(H_xbybz);
    
    // Local
    H_local.zeros();
    for(int idx=0; idx < n_sublattices(); ++idx){
      int idx1, idx2;
      std::tie(idx1, idx2) = rot.get_index_pair(idx);
      H_local += rot.elems(idx1, idx1, rot.so.H_loc_down);
      H_local += rot.elems(idx2, idx2, rot.so.H_loc_up);
    }

    /* Energy of the Jeff=1/2 orbital */
    H_local_J05 = rot.so.E_J05() * arma::eye<cx_mat_type>(n/3, n/3);
  }
  
  std::tuple<cx_mat_type, cx_mat_type> Hamiltonian::separate_mat(cx_mat_type const& mat) const {
    int size = mat.n_rows;

    // Spin nonflip
    cx_mat_type mat_nonflip(size, size, arma::fill::zeros);        
    mat_nonflip(0,0) = mat(0,0);
    mat_nonflip(0,1) = mat(0,1);
    mat_nonflip(1,0) = mat(1,0);
    mat_nonflip(1,1) = mat(1,1);
    mat_nonflip(2,2) = mat(2,2);    

    // Spin flip    
    cx_mat_type mat_flip(size, size, arma::fill::zeros);
    mat_flip(0,2) = mat(0,2);
    mat_flip(1,2) = mat(1,2);
    mat_flip(2,0) = mat(2,0);
    mat_flip(2,1) = mat(2,1);

    return std::make_tuple(mat_nonflip, mat_flip);
  }

  cx_mat_type Hamiltonian::H_elem(int i, int j, cx_mat_type const& mat_ij) const {
    // Getting the indices
    int i1, i2;
    std::tie(i1, i2) = rot.get_index_pair(i);
    int j1, j2;
    std::tie(j1, j2) = rot.get_index_pair(j);

    // Separating the matrices
    cx_mat_type tij_nonflip, tij_flip;
    std::tie(tij_nonflip, tij_flip) = separate_mat(mat_ij);

    // Adding the matrix elements    
    cx_mat_type Hij(n, n, arma::fill::zeros);
    Hij += rot.elems(i1, j1, tij_nonflip);
    Hij += rot.elems(i2, j2, tij_nonflip);    
    Hij += rot.elems(i1, j2, tij_flip);
    Hij += rot.elems(i2, j1, tij_flip);

    return Hij;
  }

  double Hamiltonian::hopping_amplitude(std::string const& bond) const {
    double tb = t;
    if (bond == "x") {
    } else if (bond == "y") {
    } else if (bond == "xy") {
      tb *= eta2;
    } else if (bond == "xyb") {
      tb *= eta2 * eta_a;
    } else if (bond == "x2") {
      tb *= eta3;
    } else if (bond == "y2") {
      tb *= eta3;
    } else if (bond == "z") {
      tb *= etaz;
    } else if (bond == "xz") {
      tb *= etaz2;
    } else if (bond == "xbz") {
      tb *= etaz2;
    } else if (bond == "yz") {
      tb *= etaz2;
    } else if (bond == "ybz") {
      tb *= etaz2;
    } else if (bond == "xyz") {
      tb *= etaz3;
    } else if (bond == "xbyz") {
      tb *= etaz3 * eta_a;
    } else if (bond == "xybz") {
      tb *= etaz3 * eta_a;
    } else if (bond == "xbybz") {
      tb *= etaz3;                        
    } else {
      std::cerr << "bond " << bond << " is not supported in " << __func__ << ".\n";
      std::exit(EXIT_FAILURE);
    }
    return tb;
  }
  
  cx_mat_type Hamiltonian::H_k(double kx, double ky, double kz) const {
    /* */
    cx_mat_type Hk(n, n, arma::fill::zeros);
    Hk += H_local;    
    Hk += hopping_amplitude("x") * 2.0 * cos(kx) * H_x;
    Hk += hopping_amplitude("y") * 2.0 * cos(ky) * H_y;
    Hk += hopping_amplitude("xy") * 2.0 * cos(kx + ky) * H_xy;
    Hk += hopping_amplitude("xyb") * 2.0 * cos(kx - ky) * H_xyb;
    Hk += hopping_amplitude("x2") * 2.0 * cos(2.*kx) * H_x2;
    Hk += hopping_amplitude("y2") * 2.0 * cos(2.*ky) * H_y2;
    if (single_layer) { return Hk; }
    
    double z_factor = 2.0;
    if (bilayer) { z_factor = 1.0; }
    Hk += hopping_amplitude("z") * z_factor * cos(kz) * H_z;
    Hk += hopping_amplitude("xz") * z_factor * cos(kx + kz) * H_xz;
    Hk += hopping_amplitude("xbz") * z_factor * cos(- kx + kz) * H_xbz;
    Hk += hopping_amplitude("yz") * z_factor * cos(ky + kz) * H_yz;
    Hk += hopping_amplitude("ybz") * z_factor * cos(- ky + kz) * H_ybz;

    Hk += hopping_amplitude("xyz") * z_factor * cos(kx + ky + kz) * H_xyz;
    Hk += hopping_amplitude("xbyz") * z_factor * cos(- kx + ky + kz) * H_xbyz;
    Hk += hopping_amplitude("xybz") * z_factor * cos(kx - ky + kz) * H_xybz;
    Hk += hopping_amplitude("xbybz") * z_factor * cos(- kx - ky + kz) * H_xbybz;    
    
    return Hk;
  }

  cx_mat_type Hamiltonian::H_k_J05(double kx, double ky, double kz) const {
    cx_mat_type Hk(n/3, n/3, arma::fill::zeros);
    Hk += H_local_J05;    
    Hk += hopping_amplitude("x") * 2.0 * cos(kx) * H_x_J05;
    Hk += hopping_amplitude("y") * 2.0 * cos(ky) * H_y_J05;
    Hk += hopping_amplitude("xy") * 2.0 * cos(kx + ky) * H_xy_J05;
    Hk += hopping_amplitude("xyb") * 2.0 * cos(kx - ky) * H_xyb_J05;
    Hk += hopping_amplitude("x2") * 2.0 * cos(2.*kx) * H_x2_J05;
    Hk += hopping_amplitude("y2") * 2.0 * cos(2.*ky) * H_y2_J05;
    if (single_layer) { return Hk; }
    
    double z_factor = 2.0;
    if (bilayer) { z_factor = 1.0; }
    Hk += hopping_amplitude("z") * z_factor * cos(kz) * H_z_J05;
    Hk += hopping_amplitude("xz") * z_factor * cos(kx + kz) * H_xz_J05;
    Hk += hopping_amplitude("xbz") * z_factor * cos(- kx + kz) * H_xbz_J05;
    Hk += hopping_amplitude("yz") * z_factor * cos(ky + kz) * H_yz_J05;
    Hk += hopping_amplitude("ybz") * z_factor * cos(- ky + kz) * H_ybz_J05;

    Hk += hopping_amplitude("xyz") * z_factor * cos(kx + ky + kz) * H_xyz_J05;
    Hk += hopping_amplitude("xbyz") * z_factor * cos(- kx + ky + kz) * H_xbyz_J05;
    Hk += hopping_amplitude("xybz") * z_factor * cos(kx - ky + kz) * H_xybz_J05;
    Hk += hopping_amplitude("xbybz") * z_factor * cos(- kx - ky + kz) * H_xbybz_J05;        
    
    /* Adding the mean field potential */
    Hk += Umf;
    
    return Hk;
  }  

  cx_mat_type Hamiltonian::extract_J05(cx_mat_type const& h) const {
    cx_mat_type h2 = rot.U_J * h * rot.U_J.t();    
    return h2.submat(index_J05, index_J05);
  }
  
  void Hamiltonian::get_hopping_parameters(){
    /* x */
    hop.tx = hopping_amplitude("x") * round0(std::real(H_x_J05(0,1)));
    hop.tx_bar = hopping_amplitude("x") * round0(std::imag(H_x_J05(0,1)));
    hop.tfx = hopping_amplitude("x") * round0(std::real(H_x_J05(0,3)));
    hop.tfx_bar = hopping_amplitude("x") * round0(std::imag(H_x_J05(0,3)));
    
    /* y */
    hop.ty = hopping_amplitude("y") * round0(std::real(H_y_J05(0,1)));
    hop.ty_bar = hopping_amplitude("y") * round0(std::imag(H_y_J05(0,1)));
    hop.tfy = hopping_amplitude("y") * round0(std::real(H_y_J05(0,3)));
    hop.tfy_bar = hopping_amplitude("y") * round0(std::imag(H_y_J05(0,3)));        

    /* xy */
    hop.txy = hopping_amplitude("xy") * round0(std::real(H_xy_J05(0,0)));

    /* xyb */
    hop.txyb = hopping_amplitude("xyb") * round0(std::real(H_xyb_J05(0,0)));

    /* x2 depends on sublattices. */
    hop.tx2_0 = hopping_amplitude("x2") * round0(std::real(H_x2_J05(0,0)));
    hop.tx2_1 = hopping_amplitude("x2") * round0(std::real(H_x2_J05(1,1)));    

    /* y2 depends on sublattices. */    
    hop.ty2_0 = hopping_amplitude("y2") * round0(std::real(H_y2_J05(0,0)));
    hop.ty2_1 = hopping_amplitude("y2") * round0(std::real(H_y2_J05(1,1)));        

    /* z */
    hop.tz = hopping_amplitude("z") * round0(std::real(H_z_J05(0,5)));
    hop.tz_bar = hopping_amplitude("z") * round0(std::imag(H_z_J05(0,5)));
    hop.tfz = hopping_amplitude("z") * round0(std::real(H_z_J05(0,7)));
    hop.tfz_bar = hopping_amplitude("z") * round0(std::imag(H_z_J05(0,7)));

    /* Interlayer hopping contributions are mixed in this four-sublattice formulation. */
  }

  void Hamiltonian::print_H_J05_hopping() const {
    /* Energy of the Jeff=1/2 orbital. */
    std::cout << "E_J05 = " << rot.so.E_J05() << std::endl;

    /* Tight-binding Hamiltonian for each bond */
    std::cout << "H_x_J05" << std::endl;
    std::cout << hopping_amplitude("x") * H_x_J05 << std::endl;

    std::cout << "H_y_J05" << std::endl;
    std::cout << hopping_amplitude("y") * H_y_J05 << std::endl;

    std::cout << "H_xy_J05" << std::endl;
    std::cout << hopping_amplitude("xy") * H_xy_J05 << std::endl;

    std::cout << "H_xyb_J05" << std::endl;
    std::cout << hopping_amplitude("xyb") * H_xyb_J05 << std::endl;

    std::cout << "H_x2_J05" << std::endl;
    std::cout << hopping_amplitude("x2") * H_x2_J05 << std::endl;

    std::cout << "H_y2_J05" << std::endl;
    std::cout << hopping_amplitude("y2") * H_y2_J05 << std::endl;

    std::cout << "H_z_J05" << std::endl;
    std::cout << hopping_amplitude("z") * H_z_J05 << std::endl;

    std::cout << "H_xz_J05" << std::endl;
    std::cout << hopping_amplitude("xz") * H_xz_J05 << std::endl;

    std::cout << "H_xbz_J05" << std::endl;
    std::cout << hopping_amplitude("xbz") * H_xbz_J05 << std::endl;

    std::cout << "H_yz_J05" << std::endl;
    std::cout << hopping_amplitude("yz") * H_yz_J05 << std::endl;

    std::cout << "H_ybz_J05" << std::endl;
    std::cout << hopping_amplitude("ybz") * H_ybz_J05 << std::endl;

    std::cout << "H_xyz_J05" << std::endl;
    std::cout << hopping_amplitude("xyz") * H_xyz_J05 << std::endl;

    std::cout << "H_xbyz_J05" << std::endl;
    std::cout << hopping_amplitude("xbyz") * H_xbyz_J05 << std::endl;

    std::cout << "H_xybz_J05" << std::endl;
    std::cout << hopping_amplitude("xybz") * H_xybz_J05 << std::endl;

    std::cout << "H_xbybz_J05" << std::endl;
    std::cout << hopping_amplitude("xbybz") * H_xbybz_J05 << std::endl;
  }

  void Hamiltonian::output_hopping_parameters(ofstream &out) const {
    hop.output(out);
  }

  void Hamiltonian::output_magnetization(ofstream &out) const {
    /* O(3) rotation is applied for each sublattice. */    
    out << std::setprecision(prec);
    out << "# Sublattice   n   mx   my   mz   |m|" << std::endl;
    for(int i=0; i < n_sublattices(); ++i){
      vec_type mag = rot.O3_spin_rotation(i) * mf.mags[i].to_arma();
      out << i << std::setw(spc) << round0(mf.ns[i]) << std::setw(spc) << round0(mag[0]) << std::setw(spc) << round0(mag[1]) << std::setw(spc) << round0(mag[2]) << std::setw(spc) << round0(arma::norm(mag)) << std::endl;
    }
  }

  void Hamiltonian::output_magnetization_in_basis(ofstream &out, mat_type const& R) const {
    /* New bases: (b1 b2 b3)^t = R (ex ey ez)^t. */
    assert(R.n_rows == 3 && R.n_cols == 3);    
    mat_type R_inv_t = arma::inv(R).t();
    
    /* O(3) rotation is applied for each sublattice. */    
    out << std::setprecision(prec);
    out << "# Sublattice   n   m_b1   m_b2   m_b3   |m|" << std::endl;
    for(int i=0; i < n_sublattices(); ++i){
      vec_type mag = R_inv_t * rot.O3_spin_rotation(i) * mf.mags[i].to_arma();  // Multiplied by R_inv_t
      out << i << std::setw(spc) << round0(mf.ns[i]) << std::setw(spc) << round0(mag[0]) << std::setw(spc) << round0(mag[1]) << std::setw(spc) << round0(mag[2]) << std::setw(spc) << round0(arma::norm(mag)) << std::endl;
    }
  }
  
  cx_mat_type Hamiltonian::onsite_U_mf(int i) const {
    double n = mf.ns[i];
    double mx = mf.mags[i].x;
    double my = mf.mags[i].y;
    double mz = mf.mags[i].z;
    cx_mat_type U_mf = {{0.5 * n - mz, - (mx - 1i*my)}, {- (mx + 1i*my), 0.5 * n + mz}};
    U_mf *= U;
    return U_mf;
  }

  void Hamiltonian::set_U_mf(){
    Umf.zeros();
    for(int i=0; i < n_sublattices(); ++i){
      cx_mat_type Umfi = onsite_U_mf(i);
      int i1, i2;
      std::tie(i1, i2) = rot.get_index_pair(i);
      Umf(i1,i1) += Umfi(0,0);
      Umf(i1,i2) += Umfi(0,1);
      Umf(i2,i1) += Umfi(1,0);
      Umf(i2,i2) += Umfi(1,1);
    }
  }

  void Hamiltonian::set_U(double _U){
    U = _U;
    set_U_mf();
  }
  
  void Hamiltonian::set_mf(MeanField const& mf2){
    mf = mf2;
    set_U_mf();
  }
}

#endif // _HAMILTONIAN_TBC_H_
