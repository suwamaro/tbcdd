/***************************************************************************
* Tight-binding calculation
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#include "hamiltonian_tbc.h"
#include "wavevector.h"
#include "model.h"

std::unique_ptr<TBC::Model> mk_model(TBC::Parameters const& pr){
  std::unique_ptr<TBC::Model> m;
  if (pr.lattice == "square") {    
    std::unique_ptr<TBC::ModelSquare> m2 = TBC::ModelSquare::mk_model_square(pr);
    m = std::move(m2);
  } else {
    std::cerr << "\"lattice\" " << pr.lattice << " is not supported.\n";
    std::exit(EXIT_FAILURE);
  }
  return m;
}

int main(int argc, char *argv[]){
  /* Parameter-file name */
  std::string input_filename = "config.toml";  
  
  /* Simulation directory */
  if ( argc < 2 ) {
    std::cerr << "Input a path to " << input_filename << ".\n";
    std::exit(EXIT_FAILURE);
  }

  /* Input parameters */  
  path base_dir(argv[1]);
  auto input_file = base_dir / input_filename;
  if ( !exists(input_file) ) {
    std::cerr << "Required input file " << input_file << " does not exist.\n";
    std::exit(EXIT_FAILURE);
  }
  TBC::Parameters pr(input_file);
  
  /* Output file */
  ofstream out_bands;
  out_bands.open(base_dir / "bands.text");
  out_bands << "# kidx  kx  ky  kz   E0  E1  E2  ..." << std::endl;
  ofstream out_bands_J05;
  out_bands_J05.open(base_dir / "bands_J05.text");
  out_bands_J05 << "# kidx  kx  ky  kz   E0  E1  E2  ..." << std::endl;
  ofstream out_bands_J05_U;  // With U
  out_bands_J05_U.open(base_dir / "bands_J05_U.text");
  out_bands_J05_U << "# kidx  kx  ky  kz   E0  E1  E2  ..." << std::endl;  
  
  /* Hamiltonian */
  TBC::Hamiltonian H(pr);

  /* Output the eigenvalues of the spin-orbit coupling. */
  ofstream out_so;
  out_so.open(base_dir / "eigenvalue_of_spin_orbit_coupling.text");
  H.rot.so.output_eigenvalue(out_so);
  out_so.close();
      
  /* Calculating bands */
  TBC::WaveVectorPath wv(pr);
  double kx, ky, kz;
  
  auto diag_output = [&](cx_mat_type const& H, ofstream& out){
    vec_type E;
    cx_mat_type V;
    arma::eig_sym(E, V, H);
    
    /* Output */
    out << wv.idx - 1 << " " << kx << " " << ky << " " << kz;
    for(double e: E){
      out << " " << e;
    }
    out << std::endl;
  };

  /* t2g */
  while (wv.unfinished()) {
    /* Wave vector */
    wv.get_wave_vector(&kx, &ky, &kz);
    
    /* Diagonalization for each wave vector */
    cx_mat_type Hk = H.H_k(kx, ky, kz);
    diag_output(Hk, out_bands);
  };

  /* Jeff=1/2 */
  /* Resetting some parameters to reproduce the t2g bands. */
  H.eta2 = pr.eta2_J05;
  H.eta3 = pr.eta3_J05;  

  /* Hopping parameters of Jeff=1/2 */
  ofstream out_hop_param_J05;
  out_hop_param_J05.open(base_dir / "hopping_parameters_J05.text");
  H.get_hopping_parameters();
  H.output_hopping_parameters(out_hop_param_J05);
  H.print_H_J05_hopping();
  out_hop_param_J05.close();
  
  /* Solving the self-consistent equation */
  bool solved = false;
  if (pr.solve_self_consistent_equation) {
    std::unique_ptr<TBC::Model> model = mk_model(pr);
    solved = model->solve_self_consistent_equation();
    if (solved) {
      /* Setting the solution. */      
      std::cout << "A self-consistent solution was found.\n";
      H.set_mf(model->H.mf);

      /* Output the result. */
      ofstream out_sc;
      out_sc.open(base_dir / "self_consistent_solution.text");
      model->output_result(out_sc);
      out_sc.close();

      /* Output the free energy. */
      ofstream out_ene;
      out_ene.open(base_dir / "free_energy.text");
      model->output_energy(out_ene);
      out_ene.close();            
    } else {
      std::cout << "No self-consistent solution was found.\n";
    }
  }
      
  /* Output magnetization */
  ofstream out_mag;
  out_mag.open(base_dir / "magnetization.text");
  H.output_magnetization(out_mag);
  out_mag.close();

  /* Output magnetization in a different basis set. */
  ofstream out_mag_b;
  out_mag_b.open(base_dir / "magnetization_in_basis.text");
  double u = 1. / sqrt(2.);
  mat_type Rb = {{u, -u, 0}, {u, u, 0}, {0, 0, 1.}};  // a b c direction
  H.output_magnetization_in_basis(out_mag_b, Rb);
  out_mag_b.close();  

  /* Basis set */
  ofstream out_basis;
  out_basis.open(base_dir / "basis_set.text");  
  out_basis << "Basis set:" << std::endl;
  out_basis << Rb << std::endl;
  out_basis.close();
  
  /* Resetting WaveVectorPath */
  wv.reset();
  
  while (wv.unfinished()) {
    /* Wave vector */
    wv.get_wave_vector(&kx, &ky, &kz);
    
    /* Jeff=1/2 without U */
    H.set_U(0);
    cx_mat_type Hk_J05 = H.H_k_J05(kx, ky, kz);    
    diag_output(Hk_J05, out_bands_J05);

    /* Jeff=1/2 with U */
    H.set_U(pr.U);    
    cx_mat_type Hk_J05_U = H.H_k_J05(kx, ky, kz);    
    diag_output(Hk_J05_U, out_bands_J05_U);    
  };
  
  /* Closing files */
  out_bands.close();
  out_bands_J05.close();
  out_bands_J05_U.close();    
  
  return EXIT_SUCCESS;
}
