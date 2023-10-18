/*****************************************************************************
*
* Functions for input parameters.
*
* Copyright (C) 2020 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <iostream>
#include <sstream>
#include <algorithm>
#include "parameters_util.h"

/* Find */
template<class T> bool is_elem_in_vector(std::vector<T> const& v, T elem){
  return std::find(v.begin(), v.end(), elem) != v.end();
}

/* Instantiation */
template<> bool is_elem_in_vector<int64_t>(std::vector<int64_t> const& v, int64_t elem);
template<> bool is_elem_in_vector<double>(std::vector<double> const& v, double elem);
  
bool is_elem_in_vector(std::vector<std::string> const& v, const char* elem){
  return is_elem_in_vector<std::string>(v, std::string(elem));
}


/* Split */
std::vector<std::string> split(std::string const& s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    if (!item.empty()) {
      elems.push_back(item);
    }
  }
  return elems;
}

std::vector<std::string> split_string(std::string str, char sep) {
    int first = 0;
    int last = str.find_first_of(sep);
    std::vector<std::string> result;    
    if ( last == std::string::npos ) {
      result.push_back(str);
    } else {
      while (first < str.size()) {
	std::string subStr(str, first, last - first);
	result.push_back(subStr);
	first = last + 1;
	last = str.find_first_of(sep, first);
	if (last == std::string::npos) {
	  last = str.size();
	}
      }
    }
    return result;
}
