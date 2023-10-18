/*****************************************************************************
*
* Parameters for tight-binding calculations.
*
* Copyright (C) 2023 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _PARAMETERS_TBC_H_
#define _PARAMETERS_TBC_H_

#include <string>
#include "cpptoml.h"
#include "cpptoml_helper.h"

namespace TBC {
  class Parameters {    
  public:
    Parameters(){};
    ~Parameters(){};
    explicit Parameters(std::string const& ifn);
    int n_sublattices = 4;
    std::string lattice;
    int n_kpoints;
    int k_route_type;
    bool single_layer, bilayer;
    double t;  // Energy scale of the hopping amplitude
    double U;  // Onsite Coulomb interaction in the mean field approximation
    std::vector<double> thetas;
    std::vector<double> phis;    
    double alpha_mono;
    double lambda, Delta;  // Spin-orbit coupling and the tetragonal distotion of octahedra
    double eta2, eta_a, eta3, etaz, etaz2, etaz3;  // Distance factors for hopping amplitudes
    double eta2_J05, eta3_J05;  // Distance factors for Jeff=1/2
    std::vector<double> m1, m2, m3, m4;  // Magnetization for mean field
    bool solve_self_consistent_equation;
    bool verbose_sc;
    bool random_initialization;
    unsigned int seed;
    double tol_sc;
    int L;
    std::size_t max_iter;
    double filling;    
  };

  Parameters::Parameters(std::string const& ifn){
    auto config = cpptoml::parse_file(ifn);
    lattice = config->get_as<std::string>("lattice").value_or("square");    
    n_kpoints = config->get_as<int>("n_kpoints").value_or(8);
    k_route_type = config->get_as<int>("k_route_type").value_or(0);    
    single_layer = config->get_as<bool>("single_layer").value_or(false);    
    bilayer = config->get_as<bool>("bilayer").value_or(true);
    t = config->get_as<double>("t").value_or(1.0);        
    U = config->get_as<double>("U").value_or(0.0);
    if (config->contains("thetas")){ thetas = get_string<double>(*config, "thetas"); }
    else { thetas.resize(4, 0.0); }
    if (config->contains("phis")){ phis = get_string<double>(*config, "phis"); }
    else { phis.resize(4, 0.0); }
    assert(thetas.size() == phis.size());
    alpha_mono = config->get_as<double>("alpha_mono").value_or(0.0);
    lambda = config->get_as<double>("lambda").value_or(0.0);
    Delta = config->get_as<double>("Delta").value_or(0.0);
    eta2 = config->get_as<double>("eta2").value_or(0.0);
    eta_a = config->get_as<double>("eta_a").value_or(0.0);
    eta3 = config->get_as<double>("eta3").value_or(0.0);
    etaz = config->get_as<double>("etaz").value_or(0.0);
    etaz2 = config->get_as<double>("etaz2").value_or(0.0);
    etaz3 = config->get_as<double>("etaz3").value_or(0.0);
    eta2_J05 = config->get_as<double>("eta2_J05").value_or(0.0);
    eta3_J05 = config->get_as<double>("eta3_J05").value_or(0.0);        
    std::vector<double> mzero = {0,0,0};
    if (config->contains("mag1")){ m1 = get_string<double>(*config, "mag1"); }
    else { m1 = mzero; }
    if (config->contains("mag2")){ m2 = get_string<double>(*config, "mag2"); }
    else { m2 = mzero; }    
    if (config->contains("mag3")){ m3 = get_string<double>(*config, "mag3"); }
    else { m3 = mzero; }    
    if (config->contains("mag4")){ m4 = get_string<double>(*config, "mag4"); }
    else { m4 = mzero; }
    
    solve_self_consistent_equation = config->get_as<bool>("solve_self_consistent_equation").value_or(true);
    verbose_sc = config->get_as<bool>("verbose_sc").value_or(true);
    random_initialization = config->get_as<bool>("random_initialization").value_or(true);
    seed = config->get_as<unsigned int>("seed").value_or(0);            
    tol_sc = config->get_as<double>("tol_sc").value_or(1e-20);            
    L = config->get_as<int>("L").value_or(8);
    max_iter = config->get_as<std::size_t>("max_iter").value_or(100000);
    filling = config->get_as<double>("filling").value_or(0.5);        
  }
}

#endif // _PARAMETERS_TBC_H_
