/***************************************************************************
* Wave vectors
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _WAVE_VECTOR_H_
#define _WAVE_VECTOR_H_

// #include "tbc.h"
#include "parameters_tbc.h"

namespace TBC
{
  /* WaveVectorPath */  
  class WaveVectorPath {
  public:
    WaveVectorPath();
    explicit WaveVectorPath(Parameters const& pr);
    bool unfinished() const;
    void get_wave_vector(double *kx, double *ky, double *kz);    
    void init(Parameters const& pr);
    void reset();
    ~WaveVectorPath(){}
    
    std::string lattice;
    std::size_t idx;
    std::size_t n_kpoints;
    std::size_t n_tot_kpoints;

    /* route_type for square lattice */
    /* 0: (0,0) -> (pi/2,pi/2) -> (pi,0) -> (0,0) */
    /* 1: (0,0) -> (pi,0) -> (pi/2,-pi/2) -> (pi/2,pi/2) -> (0,0) -> (pi/2,-pi/2) */    
    int route_type;
  };

  /* Member functions of WaveVectorPath */
  WaveVectorPath::WaveVectorPath(Parameters const& pr){
    init(pr);
  }

  void WaveVectorPath::init(Parameters const& pr){
    lattice = pr.lattice;
    idx = 0;
    n_kpoints = pr.n_kpoints;
    route_type = pr.k_route_type;

    if (lattice == "square") {
      if (route_type == 0) {
	n_tot_kpoints = 4 * n_kpoints + 1;
      } else if (route_type == 1) {
	n_tot_kpoints = 7 * n_kpoints + 1;
      } else {
	std::cerr << "\"route_type\" " << route_type << " is not supported.\n";
	std::exit(EXIT_FAILURE);
      }
    } else {
      std::cerr << "\"lattice\" " << lattice << " is not supported.\n";
      std::exit(EXIT_FAILURE);
    }
  }

  void WaveVectorPath::reset(){
    idx = 0;
  }
  
  bool WaveVectorPath::unfinished() const {
    return idx < n_tot_kpoints;
  }
  
  void WaveVectorPath::get_wave_vector(double *kx, double *ky, double *kz){
    if (lattice == "square") {    
      double dk = 0.5 * M_PI / n_kpoints;
      if (route_type == 0) {
	*kz = 0;
	if (idx < n_kpoints) {
	  *kx = idx * dk;
	  *ky = *kx;
	} else if (idx < 2 * n_kpoints) {
	  *kx = 0.5 * M_PI + (idx - n_kpoints) * dk;
	  *ky = M_PI - *kx;
	} else {
	  *kx = M_PI - (idx - 2 * n_kpoints) * dk;
	  if (std::abs(*kx) < 1e-12) { *kx = 0.0; }
	  *ky = 0.0;
	}
      } else if (route_type == 1) {
	*kz = 0;
	if (idx < 2 * n_kpoints) {
	  *kx = idx * dk;
	  *ky = 0;
	} else if (idx < 3 * n_kpoints) {
	  int di = idx - 2 * n_kpoints;
	  *kx = M_PI - di * dk;
	  *ky = - di * dk;
	} else if (idx < 5 * n_kpoints) {
	  int di = idx - 3 * n_kpoints;
	  *kx = 0.5 * M_PI;
	  *ky = - 0.5 * M_PI + di * dk;
	} else if (idx < 6 * n_kpoints) {
	  int di = idx - 5 * n_kpoints;
	  *kx = 0.5 * M_PI - di * dk;
	  *ky = *kx;
	} else {
	  int di = idx - 6 * n_kpoints;
	  *kx = di * dk;
	  *ky = - *kx;
	}
      }
    }
    ++idx;
  }

  /* WaveVector */
  class WaveVector {
  public:
    WaveVector(){};
    virtual ~WaveVector(){}
    
    std::size_t Nu;  // Number of unit cells
    double eps = 1e-12;
    std::size_t n_tot_kpoints;
    std::vector<double> factor_k;
    std::vector<vec3> ks;
  };

  
  /* WaveVectorSquare */
  class WaveVectorSquare : public WaveVector {
  public:
    WaveVectorSquare(){};
    explicit WaveVectorSquare(Parameters const& pr);
    ~WaveVectorSquare(){}    
    void init(Parameters const& pr);
    bool in_BZ(double kx, double ky) const;
    bool on_BZ(double kx, double ky) const;
    bool on_BZ_corner(double kx, double ky) const;    
    std::tuple<double, double> index_to_vectors(std::size_t idx) const;
    void set_vectors();
    
    int L;
    double k1;

    /* Instantiation */
    static std::unique_ptr<WaveVectorSquare> mk_wavevector_square(Parameters const& pr);
  };

  std::unique_ptr<WaveVectorSquare> WaveVectorSquare::mk_wavevector_square(Parameters const& pr) {
    return std::make_unique<WaveVectorSquare>(pr);
  }
  
  /* Member functions of WaveVectorSquare */
  WaveVectorSquare::WaveVectorSquare(Parameters const& pr):WaveVector(){
    init(pr);
  }

  void WaveVectorSquare::init(Parameters const& pr){
    L = pr.L;
    Nu = L * L / 2;  // Number of unit cells    
    n_tot_kpoints = (L+1) * (L+1);
    k1 = 2.0 * M_PI / (double)L;
    set_vectors();
  }  

  bool WaveVectorSquare::in_BZ(double kx, double ky) const {
    return std::abs(kx) + std::abs(ky) < M_PI + eps;
  }

  bool WaveVectorSquare::on_BZ(double kx, double ky) const {
    double k = std::abs(kx) + std::abs(ky);
    return M_PI - eps < k && k < M_PI + eps;
  }

  bool WaveVectorSquare::on_BZ_corner(double kx, double ky) const {
    bool kx_0 = std::abs(kx) < eps;
    bool ky_0 = std::abs(ky) < eps;    
    bool kx_pi = std::abs(std::abs(kx) - M_PI) < eps;
    bool ky_pi = std::abs(std::abs(ky) - M_PI) < eps;
    return (kx_0 && ky_pi) || (kx_pi && ky_0);
  }  

  std::tuple<double, double> WaveVectorSquare::index_to_vectors(std::size_t idx) const {
    int x = idx % (L+1);
    int y = idx / (L+1);
    double kx = k1 * (x - L/2);
    double ky = k1 * (y - L/2);
    return std::make_tuple(kx, ky);
  }
  
  void WaveVectorSquare::set_vectors(){
    factor_k.clear();
    ks.clear();    
    for(std::size_t idx = 0; idx < n_tot_kpoints; ++idx){
      double kx, ky;
      std::tie(kx, ky) = index_to_vectors(idx);
      if (in_BZ(kx, ky)) {
	ks.push_back(vec3(kx, ky, 0));
	if (on_BZ(kx, ky)) {
	  if (on_BZ_corner(kx, ky)) {
	    factor_k.push_back(0.25);
	  } else {
	    factor_k.push_back(0.5);
	  }
	} else {
	  factor_k.push_back(1.0);
	}
      }
    }
  }
}

#endif // _WAVE_VECTOR_H_
