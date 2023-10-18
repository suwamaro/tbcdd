/*****************************************************************************
*
* Functions for input parameters.
*
* Copyright (C) 2020 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __PARAMETERS_UTIL_
#define __PARAMETERS_UTIL_

#include <string>
#include <vector>

/* Find */
template<class T> bool is_elem_in_vector(std::vector<T> const& v, T elem);
bool is_elem_in_vector(std::vector<std::string> const& v, const char* elem);

/* Split */
std::vector<std::string> split(std::string const& s, char delim);
std::vector<std::string> split_string(std::string str, char sep);

#endif // __PARAMETERS_UTIL_
