/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include "comparisonop.h"

namespace mcl {
  
template<class T>
T Max(
  const T scalar_a,
  const T scalar_b) noexcept
{
  if (scalar_a >= scalar_b)
  {
    return scalar_a;
  }
  else
  {
    return scalar_b;
  }
}
  
template<class T>
T Min(
  const T scalar_a,
  const T scalar_b) noexcept
{
  if (scalar_a < scalar_b)
  {
    return scalar_a;
  }
  else
  {
    return scalar_b;
  }
}


/** Equivalent to Matlab's fix(scalar) */
template<typename T>
Int Fix(const T scalar)
{
  if (scalar >= 0.0)
  {
    return (Int) floor((double) scalar);
  }
  else
  {
    return (Int) ceil((double) scalar);
  }
}

/**
 Equivalent to Matlab's sign. Returns 1 if the element
 is greater than zero, 0 if it equals zero and -1 if it is less than zero.
 */
template<typename T>
Int Sign(const T scalar)
{
  if (IsApproximatelyEqual(scalar, 0.0, std::numeric_limits<T>::epsilon()))
  {
    return 0;
  }
  else if (scalar > 0.0)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}


/** Equivalent to Matlab's rem(scalar_a,scalar_b) */
template<typename T>
T Rem(
  const T x,
  const T y)
{
  if (IsApproximatelyEqual(y, (T) 0.0, std::numeric_limits<T>::epsilon()))
  {
    return x;
  }
  if (IsApproximatelyEqual(x, y, std::numeric_limits<T>::epsilon()))
  {
    return static_cast<T>(0.0);
  }
  Int n = Fix(x/y);
  return x - ((T) n)*y;
}


template<typename T>
T Floor(const T input)
{
  // TODO: verify for all types.
  return std::floor(input);
}


/** Equivalent to Matlab's mod(scalar_a,scalar_b) */
template<typename T>
T Mod(
  const T x,
  const T y)
{
  static_assert(! std::is_integral<T>::value);
  if (IsApproximatelyEqual(y, 0, std::numeric_limits<T>::epsilon()))
  {
    return x;
  }
  if (IsApproximatelyEqual(x, y, std::numeric_limits<T>::epsilon()))
  {
    return static_cast<T>(0.0);
  }
  Int signum(Sign(x/y));
  if (signum == 1 || signum == 0)
  {
    return Rem(x, y);
  }
  else
  {
    return Rem(x, y);
  }
}

template<>
int Mod<int>(
  const int x,
  const int y)
{
  if (y == 0)
  {
    return x;
  }
  if (x == y || x == -y) {
    return 0;
  }
  int n(static_cast<int>(Fix(((Real) x)/((Real) y))));
  int signum(static_cast<int>(Sign(((Real) x)/((Real) y))));
  if (signum == 1 || signum == 0) {
    return x - n*y;
  }
  else {
    return x - n*y + y;
  }
}

template<>
Int Mod<Int>(
  const Int x,
  const Int y)
{
  if (y == 0)
  {
    return x;
  }
  if (x == y || x == -y) {
    return 0;
  }
  Int n = Fix(((Real) x)/((Real) y));
  Int signum(Sign(((Real) x)/((Real) y)));
  if (signum == 1 || signum == 0) {
    return x - n*y;
  }
  else {
    return x - n*y + y;
  }
}


/** Equivalent to Matlab's abs(scalar) */
Real Abs(Real input);
  
/** Equivalent to Matlab's abs(scalar) */
Real Abs(Complex<Real> input);
  
/** Power function. Equivalent to Matlab's input^exponent. */
Real Pow(Real input, Real exponent);
  
/** Square root function. Equivalent to Matlab's sqrt(input) */
Real Sqrt(Real input);
  
/** Equivalent to Matlab's round(input). This is faster than the standard
 C++ round function, especially on Windows. Returns an integer. */
inline Int RoundToInt(Real input) {
  Int output = static_cast<int>(input);
  output += (input-output >= 0.5) - (input-output <= -0.5);
  return output;
}
  
/** Returns the conjugate of the element. Equivalent to Matlab's conj(scalar). */
Complex<Real> Conj(Complex<Real> scalar);

/** Returns the real part of a complex scalar. Equivalent to Matlab's 
 real(scalar). I am calling it `RealPart' since `Real' denotes the number type */
Real RealPart(Complex<Real> scalar);
  
/** Returns the imaginary part of a complex scalar. Equivalent to Matlab's
 imag(scalar). I am calling it `ImagPart' for consistency with `RealPart' */
Real ImagPart(Complex<Real> scalar);
  
/** Equivalent to Matlab's nextpow2(input) */
Int NextPow2(Real input);
  
/** This returns the next power of 2. For instance 5=>8, 12=>16, 16=>16. */
Int Next2(Int input);
  
/** Converts a string to a double */
double StringToDouble(const std::string& s);
  
/** Equivalent to Matlab's factorial(input) */
Int Factorial(const Int input);
  
/** Linear interpolation between two values */
Real LinearInterpolation(Real x1, Real y1, Real x2, Real y2, Real x);
  
/** 
 Returns true if the imaginary part is approximately zero. The precision used
 is VERY_SMALL in equality operations, hence use only for testing.
 */
bool IsReal(const Vector<Complex<Real>>& input);
  
/** 
 Calculates the entropy of a discreate random variable with given `pdf'.
 It normalises the pdf if its sum is not 1.
 Note: this function is identical to Matlab's only for uint8 values.
 */
Real Entropy(Vector<Real> pdf, Real base);
  
bool ElementaryOpTest();
  
#ifdef MCL_LOAD_BOOST
/** 
 Returns the value of the associated Legendre polynomial of degree `n' and
 order `m' of the values x. Equivalent to the m-th value of the vector
 returned by Matlab's legendre(n, x) 
 */
Real AssociatedLegendreP(Int n, Int m, Real x);
  
/** 
 Returns the value of the spherical harmonic of degree n and order m,
 and where theta is elevation, measured as the angle formed with the z-axis,
 and phi is the aximuth, measured as the angle formed with the x-axis.
 The harmonics are defined as:
 Y_\ell^m( \theta , \varphi ) =
 \sqrt{{(2\ell+1)\over 4\pi}{(\ell-m)!\over (\ell+m)!}}  \,
 P_\ell^m ( \cos{\theta} ) \, e^{i m \varphi }
 which are orthonormal
 */
Complex<Real> SphericalHarmonic(Int n, Int m, Real theta, Real phi);
#endif
  
  
} /**< namespace mcl  */
