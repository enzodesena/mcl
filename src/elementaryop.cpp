/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "elementaryop.h"
#include "comparisonop.h"
#include "vectorop.h"
#include "vector.h"
#include <limits>

#if MCL_LOAD_BOOST
  #include "boost/math/special_functions/spherical_harmonic.hpp"
#endif

namespace mcl {



  
Complex Conj(Complex scalar) {
  return Complex(scalar.real(), -scalar.imag());
}

Real RealPart(Complex scalar) {
  return scalar.real();
}
  
Real ImagPart(Complex scalar) {
  return scalar.imag();
}
  
Int NextPow2(Real input) { return (UInt) std::ceil(log2(std::fabs((double) input))); }

  
double StringToDouble(const std::string& s) {
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
}
  
Real Abs(Real input) {
  return (Real) std::fabs((double) input);
}
  
Real Abs(Complex input) {
  return (Real) std::abs(input);
}
  
Int Next2(Int input) {
  return (UInt) pow(2, NextPow2(input));
}
  
template<typename T, size_t length>
T Entropy(Vector<T, length> pdf, Real base) {
  pdf = Multiply(pdf, 1.0/Sum(pdf));
  return -Sum(Multiply(pdf, Log(pdf)))/log(base);
}
  
  
Int Factorial(const Int input) {
  if(input <= 1) return 1;
  
  Int temp;
  temp = input * Factorial(input - 1);
  return temp;
}
  
Real LinearInterpolation(Real x0, Real y0, Real x1, Real y1, Real x) {
  Real m = (y1-y0)/(x1-x0);
  return y0+(x-x0)*m;
}

#if MCL_LOAD_BOOST
Real AssociatedLegendreP(Int n, Int m, Real x) {
  ASSERT(n >= 0); // As in Matlab we don't accept n<0
  ASSERT(m >= 0); // As in Matlab we don't accept m<0
  ASSERT(n >= m); // As in Matlab we don't accept n<m
  return boost::math::legendre_p<Real>((int) n, (int) m, x);
}
  
Complex SphericalHarmonic(Int n, Int m, Real theta, Real phi) {
  return boost::math::spherical_harmonic<Real, Real>((int) n, (int) m,
                                                     theta, phi);
}
#endif
  
  
  
} // namespace mcl
