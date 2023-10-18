/***************************************************************************
* Mean fields
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _MEAN_FIELD_H_
#define _MEAN_FIELD_H_

#include <vector>
#include "vec.h"

namespace TBC
{
  /* MeanField */
  class MeanField {
  public:
    MeanField(){}
    explicit MeanField(int n_sublattices);
    MeanField(MeanField const& mf);
    MeanField(MeanField&& mf) noexcept;
    MeanField& operator=(MeanField const& mf);
    MeanField& operator=(MeanField&& mf) noexcept;
    MeanField operator+(MeanField const& other) const;
    MeanField operator-(MeanField const& other) const;
    double norm2() const;
    double norm() const;    
    void init(int n_sublattices);
    void print() const;
    ~MeanField(){}

    std::vector<double> ns;  // Deviation of charge density
    std::vector<vec3> mags;  // Magnetization
  };

  MeanField::MeanField(int n_sublattices){
    init(n_sublattices);
  }

  MeanField::MeanField(MeanField const& mf){
    ns = mf.ns;
    mags = mf.mags;
  }

  MeanField::MeanField(MeanField&& mf) noexcept {
    ns = std::move(mf.ns);
    mags = std::move(mf.mags);
  }
  
  MeanField& MeanField::operator=(MeanField const& mf){
    if (this != &mf) {
      ns = mf.ns;
      mags = mf.mags;
    }
    return *this;
  }

  MeanField& MeanField::operator=(MeanField&& mf) noexcept {
    if (this != &mf) {
      ns = std::move(mf.ns);
      mags = std::move(mf.mags);
    }
    return *this;
  }  

  MeanField MeanField::operator+(MeanField const& other) const {
    int n_sub = mags.size();
    MeanField mf(n_sub);
    for(int i=0; i < n_sub; ++i){
      mf.ns[i] = ns[i] + other.ns[i];
    }
    for(int i=0; i < n_sub; ++i){
      mf.mags[i] = mags[i] + other.mags[i];
    }    
    
    return mf;
  }

  MeanField MeanField::operator-(MeanField const& other) const {
    int n_sub = mags.size();    
    MeanField mf(n_sub);
    for(int i=0; i < n_sub; ++i){
      mf.ns[i] = ns[i] - other.ns[i];
    }
    for(int i=0; i < n_sub; ++i){
      mf.mags[i] = mags[i] - other.mags[i];
    }    
    
    return mf;
  }  

  double MeanField::norm2() const {
    double sum = 0;
    for(int i=0; i < ns.size(); ++i){
      sum += ns[i] * ns[i];
    }
    for(int i=0; i < mags.size(); ++i){
      sum += mags[i].norm2();
    }    
    
    return sum;
  }

  double MeanField::norm() const {
    return sqrt(norm2());
  }
  
  void MeanField::init(int n_sublattices){    
    ns.resize(n_sublattices, 0.0);
    vec3 mag(0,0,0);
    mags.resize(n_sublattices, mag);
  }

  void MeanField::print() const {
    std::cout << "ns" << std::endl;
    for(int i=0; i < ns.size(); ++i){      
      std::cout << ns[i] << std::endl;
    }

    std::cout << "mags" << std::endl;    
    for(int i=0; i < mags.size(); ++i){
      std::cout << mags[i] << std::endl;
    }
  }

  /* Operator overloading */
  std::ostream& operator<<(std::ostream& os, MeanField const& mf){
    os << "ns:";
    for(int i=0; i < mf.ns.size(); ++i){
      os << "  " << mf.ns[i];
    }
    os << std::endl;
    
    os << "mags:";
    for(int i=0; i < mf.mags.size(); ++i){
      os << "  " << mf.mags[i];
    }
    
    return os;
  }
}

#endif // _MEAN_FIELD_H_
