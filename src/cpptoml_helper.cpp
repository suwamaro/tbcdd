/*****************************************************************************
*
* Helper functions for cpptoml.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <iostream>
#include "cpptoml_helper.h"
#include "parameters_util.h"

template <class T> T get_from_table(const std::shared_ptr<cpptoml::table> table, const std::string& key){
  auto v = table->get_qualified_as<T>(key);
  if (v) {
    return *v;
  }
  std::cerr << "No key '" << key << "' was found in the table.\n" << "Table list:\n" << *table << std::endl;
  std::exit(EXIT_FAILURE);
}

// Instantiation
template int64_t get_from_table<int64_t>(const std::shared_ptr<cpptoml::table> table, const std::string& key);
template double get_from_table<double>(const std::shared_ptr<cpptoml::table> table, const std::string& key);
template bool get_from_table<bool>(const std::shared_ptr<cpptoml::table> table, const std::string& key);


template<class T> std::vector<T> get_string_impl(cpptoml::table const& c, std::string const& pname){
  std::string p_string = *(c.get_as<std::string>(pname));
  std::vector<std::string> p_strings = split_string(p_string, ' ');  // Separated by ' '
  std::vector<T> p_values;
  for(int i=0; i < p_strings.size(); i++){
    p_values.push_back(static_cast<T>(std::stod(p_strings[i])));    
  }
  return p_values;
}

template<class T> std::vector<T> get_string(cpptoml::table const& c, std::string const& pname){
  try {
    return get_string_impl<T>(c, pname);
  } catch (...){
    std::cout << "Function get_string threw an exception." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

/* Instantiation */
template std::vector<int> get_string<int>(cpptoml::table const& c, std::string const& pname);
template std::vector<double> get_string<double>(cpptoml::table const& c, std::string const& pname);
