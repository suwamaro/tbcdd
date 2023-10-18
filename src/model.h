/***************************************************************************
* Mean-field models
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _MODEL_TBC_H_
#define _MODEL_TBC_H_

#include "hamiltonian_tbc.h"
#include "wavevector.h"
#include "vec.h"
#include "tbc_util.h"

namespace TBC
{
  /* Model */  
  class Model {
  public:
    Model();
    explicit Model(Parameters const& pr);
    virtual ~Model(){}
    bool solve_self_consistent_equation();
    std::size_t n_elems() const;    
    void randomize_mean_field();
    void set_eigenvectors();
    void print_n_elems() const;
    bool self_consistent(MeanField const& mean) const;
    double calc_chemical_potential(double fraction) const;
    MeanField calc_mean();
    double calc_elec_energy() const;
    double calc_classical_energy() const;        
    void output_result(ofstream &out) const;
    void output_energy(ofstream &out) const;

    std::unique_ptr<WaveVector> wv;
    bool verbose;
    bool random_initialization;
    RNG rng;
    std::size_t max_iter;
    double tol;
    Hamiltonian H;
    double filling;
    double mu = 0;  // Chemical potential    
    mat_type Eks;  // Eigenvalues
    cx_cube_type Uks;  // Eigenmatrices
  };

  /* Member functions of Model */
  Model::Model(Parameters const& pr):rng(pr.seed),H(pr){
    verbose = pr.verbose_sc;
    random_initialization = pr.random_initialization;
    max_iter = pr.max_iter;
    tol = pr.tol_sc;
    filling = pr.filling;
    if (verbose) {
      std::cout << "max_iter = " << max_iter << "  tol = " << tol << "  filling = " << filling << std::endl;
    }
  }

  std::size_t Model::n_elems() const {
    std::size_t nk = wv->ks.size();
    std::size_t Hdim = 2 * H.n_sublattices();
    return Hdim * Hdim * nk;
  }
  
  void Model::randomize_mean_field(){
    MeanField mf(H.n_sublattices());
    for(int i=0; i < H.n_sublattices(); ++i){
      mf.ns[i] = 2. * filling;
      std::uniform_real_distribution<> dist(0.0, 1.0);      
      double theta = dist(rng) * M_PI;
      double phi = dist(rng) * 2.0 * M_PI;
      double S = 0.5;      
      mf.mags[i].x = S * sin(theta) * cos(phi);
      mf.mags[i].y = S * sin(theta) * sin(phi);
      mf.mags[i].z = S * cos(theta);      
    }
    H.set_mf(mf);    
  }

  void Model::set_eigenvectors(){
    std::size_t nk = wv->ks.size();
    std::size_t Hdim = 2 * H.n_sublattices();
    Eks.resize(Hdim, nk);
    Uks.resize(Hdim, Hdim, nk);
      
    vec_type E(Hdim);
    cx_mat_type U(Hdim, Hdim);
    for(std::size_t ki=0; ki < nk; ++ki){
      vec3 k = wv->ks[ki];
      cx_mat_type Hk = H.H_k_J05(k.x, k.y, k.z);      
      arma::eig_sym(E, U, Hk);  
      Eks.col(ki) = E;
      Uks.slice(ki) = U;
    }
  }

  void Model::print_n_elems() const {
    std::cout << "Number of elements: " << n_elems() << std::endl;
  }
  
  bool Model::solve_self_consistent_equation(){
    if (verbose) {
      std::cout << "Solving the self-consistent equation..." << std::endl;
    }
    
    /* Initialization */
    std::size_t t = 0;
    MeanField mf_mean(H.n_sublattices());
    if (random_initialization) {
      randomize_mean_field();
    }
    
    do {      
      if (verbose) {
	std::cout << "t = " << t << std::endl;
	std::cout << H.mf << std::endl;	
      }
      mf_mean = calc_mean();
      ++t;
      
      if (self_consistent(mf_mean)){ return true; }
      if (t == max_iter) { return false; }
      H.set_mf(mf_mean);
    } while (true);
  }
  
  bool Model::self_consistent(MeanField const& mean) const {
    MeanField mf_diff = mean - H.mf;
    double diff2 = mf_diff.norm2();
    
    if (verbose) {
      std::cerr << "diff = " << diff2 << std::endl;
    }
    
    return diff2 < tol;
  }

  double Model::calc_chemical_potential(double fraction) const {
    /* Half-filling */
    if (std::abs(fraction - 0.5) < 1e-12) {
      return arma::median(arma::vectorise(Eks));
    }
    
    /* Starting with the previous chemical potential. */
    double _mu = mu;    
    double dmu_init_abs = 0.1;
    double dmu = 0;
    double factor = 0.51;  // Set to a value a little greater than 0.5.
    double n_elec_target = 2 * H.n_sublattices() * wv->Nu * fraction;
    while(true){
      double n_elec = 0;
      for(std::size_t ki=0; ki < wv->ks.size(); ++ki){
	arma::uword count = arma::accu(Eks.col(ki) < _mu);  // Counting the number of elements smaller than mu.
	n_elec += wv->factor_k[ki] * (double)count;
      }
      
      if (n_elec < n_elec_target - 1e-12) {
	if (dmu > 0) {
	} else if (dmu < 0) {
	  dmu *= - factor;
	} else {
	  dmu = dmu_init_abs;
	}
      } else if (n_elec > n_elec_target + 1e-12) {
	if (dmu < 0) {
	} else if (dmu > 0) {
	  dmu *= - factor;
	} else {
	  dmu = - dmu_init_abs;
	}
      } else {
	break;
      }
      _mu += dmu;
    };
    
    return _mu;
  }
  
  MeanField Model::calc_mean(){
    /* Calculating the expectation values of the charge density and spin operators. */
    set_eigenvectors();
    mu = calc_chemical_potential(filling);
    
    double kT = 0.0;  // Assume zero temperature
    MeanField mf(H.n_sublattices());
    int Hdim = 2 * H.n_sublattices();    
    cx_mat_type D(Hdim, Hdim);
    
    /* Summation over the BZ */
    for(std::size_t ki=0; ki < wv->ks.size(); ++ki){
      vec_type fdens(Hdim);
      for(int l=0; l < Hdim; ++l){
	fdens(l) = fermi_density(Eks(l,ki), kT, mu);
      }
      sp_mat_type F(arma::diagmat(fdens));
      D = Uks.slice(ki) * F * Uks.slice(ki).t();

      /* Calculating the mean value for each sublattice. */
      for(int gamma=0; gamma < H.n_sublattices(); ++gamma){
	int i1, i2;
	std::tie(i1, i2) = H.rot.get_index_pair(gamma);

	double factor = wv->factor_k[ki];
	mf.ns[gamma] += factor * std::real(D(i1, i1) + D(i2, i2));
	mf.mags[gamma].x += factor * std::real(D(i2, i1));
	mf.mags[gamma].y += factor * std::imag(D(i2, i1));
	mf.mags[gamma].z += factor * 0.5 * std::real(D(i1, i1) - D(i2, i2));
      }
    }

    /* Normalization */
    for(int gamma=0; gamma < H.n_sublattices(); ++gamma){    
      mf.ns[gamma] /= (double)wv->Nu;
      mf.mags[gamma] /= (double)wv->Nu;
    }
    
    return mf;
  }

  double Model::calc_elec_energy() const {
    double g = 0;
    double kT = 0;
  
    /* Summation over the BZ */
    std::size_t Hdim = 2 * H.n_sublattices();    
    for(std::size_t ki=0; ki < wv->ks.size(); ++ki){
      double factor = wv->factor_k[ki];
      for(int l=0; l < Hdim; ++l){
	g += factor * fermi_energy(Eks(l,ki), kT, mu);
      }
    }

    /* Normalization */
    g /= (wv->Nu * H.n_sublattices());

    return g;
  }

  double Model::calc_classical_energy() const {
    double e_cl_charge = 0;
    double e_cl_spin = 0;
    for(int i=0; i < H.n_sublattices(); ++i){
      e_cl_charge += H.mf.ns[i] * H.mf.ns[i];
      e_cl_spin += H.mf.mags[i].norm2();
    }
    e_cl_charge *= - 0.25 * H.U;
    e_cl_spin *= H.U;

    /* Normalization */
    e_cl_charge /= (double)H.n_sublattices();
    e_cl_spin /= (double)H.n_sublattices();

    return e_cl_charge + e_cl_spin;
  }
  
  void Model::output_result(ofstream &out) const {
    out << std::setprecision(prec);
    out << "# mu = " << round0(mu) << std::endl;
    out << "# Sublattice   n   mx   my   mz   |m|" << std::endl;
    for(int i=0; i < H.n_sublattices(); ++i){
      out << i << std::setw(spc) << round0(H.mf.ns[i]) << std::setw(spc) << round0(H.mf.mags[i].x) << std::setw(spc) << round0(H.mf.mags[i].y) << std::setw(spc) << round0(H.mf.mags[i].z) << std::setw(spc) << round0(H.mf.mags[i].norm()) << std::endl;
    }    
  }

  void Model::output_energy(ofstream &out) const {
    out << std::setprecision(prec);
    double g_e = calc_elec_energy();
    double g_cl = calc_classical_energy();    
    double mu_eff = mu - H.rot.so.E_J05();
    out << "Free_energy_plus_mu_eff = " << g_e + g_cl + mu_eff << std::endl;
    out << "Free_energy = " << g_e + g_cl << std::endl;
    out << "Electronic_energy = " << g_e << std::endl;
    out << "Classical_energy = " << g_cl << std::endl;        
    out << "mu_eff = " << mu_eff << std::endl;    
    out << "mu = " << mu << std::endl;
    out << "E_J05 = " << H.rot.so.E_J05() << std::endl;    
  }

  
  /* ModelSquare */
  class ModelSquare: public Model {
  public:
    ModelSquare();
    explicit ModelSquare(Parameters const& pr);
    ~ModelSquare(){}

    /* Instantiation */
    static std::unique_ptr<ModelSquare> mk_model_square(Parameters const& pr);
  };

  std::unique_ptr<ModelSquare> ModelSquare::mk_model_square(Parameters const& pr){
    return std::make_unique<ModelSquare>(pr);
  }
  
  /* Member functions of ModelSquare */
  ModelSquare::ModelSquare(Parameters const& pr):Model(pr){
    wv = WaveVectorSquare::mk_wavevector_square(pr);
    set_eigenvectors();
    if (verbose) { print_n_elems(); }
  }
  
}

#endif // _MODEL_TBC_H_
