/*****************************************************************************
*
* Helper functions for cpptoml.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/
#ifndef _CPPTOML_HELPER_
#define _CPPTOML_HELPER_

#include "cpptoml.h"

template <class T> T get_from_table(const std::shared_ptr<cpptoml::table> table, const std::string& key);

template<class T> std::vector<T> get_string_impl(cpptoml::table const& c, std::string const& pname);
template<class T> std::vector<T> get_string(cpptoml::table const& c, std::string const& pname);
  
#endif // _CPPTOML_HELPER_
