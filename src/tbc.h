/***************************************************************************
* Tight-binding calculation
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _TBC_H_
#define _TBC_H_

#include <iostream>
#include <complex>
#include <random>
#include <armadillo>

/* Using std::filesystem */
typedef std::ifstream ifstream;
typedef std::ofstream ofstream;
#include <filesystem>
typedef std::filesystem::path path;
typedef std::filesystem::directory_iterator directory_iterator;
using std::filesystem::exists;
using std::filesystem::is_directory;
using std::filesystem::is_regular_file;
using std::filesystem::create_directories;

/* Value types */
typedef std::complex<double> cx_double;
typedef arma::vec vec_type;
typedef arma::cx_vec cx_vec_type;
typedef arma::mat mat_type;
typedef arma::sp_mat sp_mat_type;
typedef arma::cx_mat cx_mat_type;
typedef arma::sp_cx_mat sp_cx_mat_type;
typedef arma::cube cube_type;
typedef arma::cx_cube cx_cube_type;

typedef std::mt19937 RNG;

constexpr int prec = 10;
constexpr int spc = 18;

#endif // _TBC_H_
