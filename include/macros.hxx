#pragma once

#include <cstdint> // uint64_t, int64_t unit32_t, int32_t, ...
#include <cmath> // std::round

  template<int nbits> size_t inline align(int64_t const in) 
          { return ((((in - 1) >> nbits) + 1) << nbits); }
  // template<int nbits> size_t inline align(int64_t const in) { return in; } // alignment switched off

  int inline required_bits(uint64_t const in) {
      int nbits = 0;
      auto n = in - 1;
      while (n) { 
          ++nbits;
          n >>= 1;
      } // while
      return nbits;
  } // required_bits

  template <typename real_t> inline char const * real_t_name(); // provide no implementation for the general case
  template <> inline char const * real_t_name<double>() { return "double"; }
  template <> inline char const * real_t_name<float> () { return "float"; }


  template<typename T> inline int constexpr sgn(T const val) {
      return (T(0) < val) - (val < T(0));
  } // sgn

  template<typename T> inline T constexpr pow2(T const x) { return x*x; } // pow2
  template<typename T> inline T constexpr pow3(T const x) { return x*pow2(x); } // pow3
  template<typename T> inline T constexpr pow4(T const x) { return pow2(pow2(x)); }
  template<typename T> inline T constexpr pow8(T const x) { return pow2(pow4(x)); }

  template<typename real_t> inline 
  real_t intpow(real_t const x, unsigned const nexp) {
      // power function using recursive doubling, only non-negative powers possible
      unsigned n = nexp;
      real_t xbin = x, xpow = (real_t)1;
      while (n) {
          if (n & 1) xpow *= xbin; // if n modulo 2 == 1
          n >>= 1; // divide n by 2
          xbin *= xbin; // square x
      } // while n is nonzero
      return xpow;
  } // intpow


  template<unsigned Step=1> inline
  double constexpr factorial(unsigned const n) {
      return (n > 1)? factorial<Step>(n - Step)*((double)n) : 1;
  } // factorial n! for Step=1 and double_factorial n!! for Step=2

  template<typename real_t> inline
  void set(real_t y[], size_t const n, real_t const a) {
      for(size_t i = 0; i < n; ++i) y[i] = a;
  } // set

  template<typename real_t, typename real_a_t> inline
  void set(real_t y[], size_t const n, real_a_t const a[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*f;
  } // set

  template<typename real_t, typename real_a_t> inline
  void scale(real_t y[], size_t const n, real_a_t const a[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] *= a[i]*f;
  } // scale

  template<typename real_t> inline
  void scale(real_t y[], size_t const n, real_t const f) {
      // no default value for f here, as scaling with 1.0 has no effect
      for(size_t i = 0; i < n; ++i) y[i] *= f;
  } // scale

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  void product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*b[i]*f;
  } // product

  template<typename real_t, typename real_a_t, typename real_b_t, typename real_c_t> inline
  void product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_c_t const c[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*b[i]*c[i]*f;
  } // product

  template<typename real_t, typename real_a_t> inline
  void add_product(real_t y[], size_t const n, real_a_t const a[], real_t const f) {
      // no default value for f here, otherwise the name add_product is missleading!
      for(size_t i = 0; i < n; ++i) y[i] += a[i]*f;
  } // add_product == axpy-type

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  void add_product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] += a[i]*b[i]*f;
  } // add_product

  template<typename real_t, typename real_a_t> inline
  double dot_product(size_t const n, real_t const bra[], real_a_t const ket[]) {
      double dot = 0;
      for(size_t i = 0; i < n; ++i) {
          dot += bra[i]*ket[i];
      } // i
      return dot;
  } // dot_product

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  double dot_product(size_t const n, real_t const bra[], real_a_t const ket[], real_b_t const metric[]) {
      double dot = 0;
      for(size_t i = 0; i < n; ++i) {
          dot += bra[i]*metric[i]*ket[i];
      } // i
      return dot;
  } // dot_product
   
  inline bool is_integer(double const f) { return (f == std::round(f)); }
