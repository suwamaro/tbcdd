tbcdd
==========================

This code produces a tight-binding model of d-orbitals and calculates the electron bands of the four-sublattice bilayer system. The model Hamiltonian includes spin-orbit coupling and tetragonal distortion as local terms. It projects out an effective Jeff=1/2 model from the t2g orbital model and calculates the magnetic order in the Hartree-Fock approximation. One of the target systems is iridium oxide. The main parameters of this code are the distance factors of electron hopping amplitudes and the octahedral rotation (theta) and tilt (phi) angles of oxygen anions for each sublattice. The input parameters must be set in a toml file (config.toml).

## Prerequisites

 - CMake - http://www.cmake.org 
 - Armadillo - http://arma.sourceforge.net
 
## Compilation

    cmake [path-to-source-directory]
    make
    
## Developer

 - Hidemaro Suwa (University of Tokyo) suwamaro@phys.s.u-tokyo.ac.jp
