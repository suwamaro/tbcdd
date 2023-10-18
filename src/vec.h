/***************************************************************************
* Vector class
*
* Copyright (C) 2023 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/

#ifndef _VEC_H_
#define _VEC_H_

#include <vector>
#include <cassert>
#include <armadillo>

/* Vec3 */
template<typename T> class Vec3 {
public:
  Vec3();
  Vec3(T x, T y, T z);
  Vec3(Vec3 const& v);
  Vec3(Vec3&& v) noexcept;
  Vec3(std::vector<T> const& v);
  Vec3(arma::Col<T> const& v);
  Vec3<T>& operator=(Vec3<T> const& v);
  Vec3<T>& operator=(Vec3<T>&& v) noexcept;
  Vec3<T>& operator=(std::vector<T> const& v);
  Vec3<T>& operator=(arma::Col<T> const& v);  
  Vec3<T> operator-() const;
  Vec3<T> operator+(Vec3<T> const& other) const;
  Vec3<T> operator-(Vec3<T> const& other) const;
  Vec3<T> operator*(T a) const;
  Vec3<T> operator/(T a) const;
  void operator+=(Vec3<T> const& other);
  void operator-=(Vec3<T> const& other);
  void operator*=(T a);
  void operator/=(T a);      
  T dot(Vec3<T> const& other) const;
  Vec3<T> cross(Vec3<T> const& other) const;
  T norm2() const;  
  T norm() const;
  Vec3<T> normalized() const;
  arma::Col<T> to_arma() const;  
  
  T x;
  T y;
  T z;
};

template<typename T> Vec3<T>::Vec3():x(0),y(0),z(0){}

template<typename T> Vec3<T>::Vec3(T x, T y, T z):x(x),y(y),z(z){}

template<typename T> Vec3<T>::Vec3(Vec3 const& v){
  x = v.x;
  y = v.y;
  z = v.z;
}

template<typename T> Vec3<T>::Vec3(Vec3&& v) noexcept {
  x = std::move(v.x);
  y = std::move(v.y);
  z = std::move(v.z);
}

template<typename T> Vec3<T>::Vec3(std::vector<T> const& v){
  assert(v.size() == 3);
  x = v[0];
  y = v[1];
  z = v[2];
}

template<typename T> Vec3<T>::Vec3(arma::Col<T> const& v){
  assert(v.n_elem == 3);
  x = v[0];
  y = v[1];
  z = v[2];
}

template<typename T> Vec3<T>& Vec3<T>::operator=(Vec3<T> const& v){
  if (this != &v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  return *this;
}

template<typename T> Vec3<T>& Vec3<T>::operator=(Vec3<T>&& v) noexcept {
  if (this != &v) {
    x = std::move(v.x);
    y = std::move(v.y);
    z = std::move(v.z);
  }
  return *this;
}

template<typename T> Vec3<T>& Vec3<T>::operator=(std::vector<T> const& v){
  assert(v.size() == 3);
  x = v[0];
  y = v[1];
  z = v[2];
  return *this;
}

template<typename T> Vec3<T>& Vec3<T>::operator=(arma::Col<T> const& v){
  assert(v.n_elem == 3);
  x = v[0];
  y = v[1];
  z = v[2];
  return *this;
}

template<typename T> Vec3<T> Vec3<T>::operator-() const {
  return Vec3(-x, -y, -z);
}

template<typename T> Vec3<T> Vec3<T>::operator+(Vec3<T> const& other) const {
  return Vec3(x+other.x, y+other.y, z+other.z);
}

template<typename T> Vec3<T> Vec3<T>::operator-(Vec3<T> const& other) const {
  return Vec3(x-other.x, y-other.y, z-other.z);
}

template<typename T> Vec3<T> Vec3<T>::operator*(T a) const {
  return Vec3(x*a, y*a, z*a);
}

template<typename T> Vec3<T> Vec3<T>::operator/(T a) const {
  return Vec3(x/a, y/a, z/a);
}

template<typename T> void Vec3<T>::operator+=(Vec3<T> const& other){
  x += other.x;
  y += other.y;
  z += other.z;
}
    
template<typename T> void Vec3<T>::operator-=(Vec3<T> const& other){
  x -= other.x;
  y -= other.y;
  z -= other.z;
}
    
template<typename T> void Vec3<T>::operator*=(T a){
  x *= a;
  y *= a;
  z *= a;
}
    
template<typename T> void Vec3<T>::operator/=(T a){
  x /= a;
  y /= a;
  z /= a;
}

template<typename T> T Vec3<T>::dot(Vec3<T> const& other) const {
  return x*other.x + y*other.y + z*other.z;
}
    
template<typename T> Vec3<T> Vec3<T>::cross(Vec3<T> const& other) const {
  return Vec3(y*other.z - z*other.y, z*other.x - x*other.z, x*other.y - y*other.x);
}

template<typename T> T Vec3<T>::norm2() const {
  return dot(*this);
}

template<typename T> T Vec3<T>::norm() const {
  return sqrt(norm2());
}

template<typename T> Vec3<T> Vec3<T>::normalized() const {
  return *this / norm();
}

template<typename T> arma::Col<T> Vec3<T>::to_arma() const {
  arma::Col<T> v = {x, y, z};
  return v;
}

/* Operator overloading */
template<typename T> std::ostream& operator<<(std::ostream& os, Vec3<T> const& v){
  return os << "(x=" << v.x << ", y=" << v.y << ", z=" << v.z << ")";
}

template<typename T> Vec3<T> operator*(T a, Vec3<T> x) {
  return x * a;
}

template<typename T> Vec3<T> operator*(Vec3<T> x, T a) {
  return x * a;
}


/* Typedef */
typedef Vec3<double> vec3;

#endif
