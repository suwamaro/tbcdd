tbcdd
==========================

This code produces a tight-binding model of d-orbitals and calculates the electron bands of the four-sublattice bilayer system. The model Hamiltonian includes spin-orbit coupling and tetragonal distortion as local terms. It also projects out an effective Jeff=1/2 model from the t2g orbital model, targeting iridium oxides. In the effective model, the on-site Coulomb repulsion is considered, constructing a Hubbard model. The magnetic order of the Hubbard model is obtained by solving the self-consistent equation in the Hartree-Fock approximation.

The main parameters of this code are the distance factors of electron hopping amplitudes and the octahedral rotation (theta) and tilt (phi) angles of oxygen anions for each sublattice. The input parameters must be set in a toml file (config.toml). The parameter details can be found in parameters_tbc.h.
	
## Prerequisites

 - CMake - http://www.cmake.org 
 - Armadillo - http://arma.sourceforge.net
 
## Compilation

    cmake [path-to-source-directory]
    make
    
## Developer

 - Hidemaro Suwa (University of Tokyo) suwamaro@phys.s.u-tokyo.ac.jp
