/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include "comparisonop.h"
#include "vectorop.h"

namespace mcl {
  
template<class T>
inline T Max(
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
inline T Min(
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
inline Int Fix(const T scalar)
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
inline Int Sign(const T scalar)
{
  if (scalar == T(0.0))
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
inline T Rem(
  const T x,
  const T y)
{
  if (IsApproximatelyEqual(y, (T) 0.0, std::numeric_limits<T>::epsilon()))
  {
    return std::numeric_limits<T>::quiet_NaN();
  }
  if (IsApproximatelyEqual(x, y, std::numeric_limits<T>::epsilon()))
  {
    return static_cast<T>(0.0);
  }
  Int n = Fix(x/y);
  return x - ((T) n)*y;
}


template<typename T>
inline T Floor(const T input)
{
  // TODO: verify for all types.
  return std::floor(input);
}


/** Equivalent to Matlab's mod(scalar_a,scalar_b) */
template<typename T>
inline T Mod(
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
    return Rem(x, y) + y;
  }
}

template<>
inline int Mod<int>(
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
  int n(static_cast<int>(Fix(((double) x)/((double) y))));
  int signum(static_cast<int>(Sign(((double) x)/((double) y))));
  if (signum == 1 || signum == 0) {
    return x - n*y;
  }
  else {
    return x - n*y + y;
  }
}

template<>
inline Int Mod<Int>(
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
  Int n = Fix(((double) x)/((double) y));
  Int signum(Sign(((double) x)/((double) y)));
  if (signum == 1 || signum == 0) {
    return x - n*y;
  }
  else {
    return x - n*y + y;
  }
}


/** Equivalent to Matlab's abs(scalar) */
inline double Abs(double input)
{
  return std::fabs(input);
}
  
  
/** Equivalent to Matlab's abs(scalar) */
inline float Abs(float input)
{
  return std::abs(input);
}


/** Equivalent to Matlab's abs(scalar) */
template<typename T>
inline T Abs(Complex<T> input)
{
  return static_cast<T>(std::abs(input));
}
  
  
/** Power function. Equivalent to Matlab's input^exponent. */
template<typename T>
inline T Pow(
  T input,
  T exponent)
{
  return (T) pow((double) input, (double) exponent);
}

/** Square root function. Equivalent to Matlab's sqrt(input) */
template<typename T>
inline T Sqrt(
  T input)
{
  return (T) sqrt((double) input);
}
  
  
/** Equivalent to Matlab's round(input). This is faster than the standard
 C++ round function, especially on Windows. Returns an integer. */
template<typename T>
inline Int RoundToInt(
  T input) noexcept
{
  Int output = static_cast<Int>(input);
  output += (input-output >= 0.5) - (input-output <= -0.5);
  return output;
}
  
/** Returns the conjugate of the element. Equivalent to Matlab's conj(scalar). */
template<typename T>
inline Complex<T> Conj(
  const Complex<T> scalar) noexcept
{
  return Complex(scalar.real(), -scalar.imag());
}

/** Returns the real part of a complex scalar. Equivalent to Matlab's 
 real(scalar). I am calling it `RealPart' since `T' denotes the number type */
template<typename T>
inline T RealPart(
  const Complex<T> scalar) noexcept
{
  return scalar.real();
}

/** Returns the imaginary part of a complex scalar. Equivalent to Matlab's
 imag(scalar). I am calling it `ImagPart' for consistency with `RealPart' */
template<typename T>
inline T ImagPart(
  const Complex<T> scalar) noexcept
{
  return scalar.imag();
}
  
/** Equivalent to Matlab's nextpow2(input) */
template<typename T>
inline Int NextPow2(
  T input) noexcept
{
  return static_cast<int>(std::ceil(log2(std::fabs((double) input))));
}
  
/** This returns the next power of 2. For instance 5=>8, 12=>16, 16=>16. */
//template<typename T>
inline Int Next2(
  Int input) noexcept
{
  return pow(2, NextPow2(input));
}
  
/** Converts a string to a double */
inline double StringToDouble(
  const std::string& s) noexcept
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
  {
    return 0;
  }
  else
  {
    return x;
  }
}
  
/** Equivalent to Matlab's factorial(input) */
template<typename T>
inline T Factorial(const T input) // TODO: inlining this is a bit dangerous
{
  if(input <= 1)
  {
    return 1;
  }
  else
  {
    return input * Factorial(input - 1);
  }
}
  
/** Linear interpolation between two values */
template<typename T>
inline T LinearInterpolation(
  T x0,
  T y0,
  T x1,
  T y1,
  T x) noexcept
{
  T m = (y1-y0)/(x1-x0);
  return y0+(x-x0)*m;
}
  
  
/**
 Returns true if the imaginary part is exactly zero
 (tests equality with the type's 0).
 */
template<typename T>
bool IsReal(
  const Complex<T>& input)
{
  return input.imag() == T(0.0);
}

/** 
 Returns true if the imaginary part is approximately zero. The precision used
 is VERY_SMALL in equality operations, hence use only for testing.
 */
template<typename T>
bool IsApproximatelyReal(
  const Complex<T>& input,
  const T precision = VERY_SMALL)
{
  return IsApproximatelyEqual(input.imag(), T(0.0), precision);
}

template<typename T>
bool IsReal(
  const Vector<T>& vector)
{
  for (auto& element : vector)
  {
    if (! IsReal(element))
    {
      return false;
    }
  }
  return true;
}

template<typename T>
bool IsApproximatelyReal(
  const Vector<T>& vector,
  const T precision = VERY_SMALL)
{
  for (auto& element : vector)
  {
    if (! IsApproximatelyReal(element, precision))
    {
      return false;
    }
  }
  return true;
}

  
/** 
 Calculates the entropy of a discreate random variable with given `pdf'.
 It normalises the pdf if its sum is not 1.
 Note: this function is identical to Matlab's only for uint8 values.
 */
template<typename T, size_t length>
T Entropy(
  const Vector<T,length> pdf,
  T base)
{
  Vector<T,length> normalised_pdf(Multiply(pdf, 1.0/Sum(pdf)));
  return -Sum(Multiply(pdf, Log(pdf)))/log(base);
}
  

  
  

#if MCL_LOAD_BOOST
/**
 Returns the value of the associated Legendre polynomial of degree `n' and
 order `m' of the values x. Equivalent to the m-th value of the vector
 returned by Matlab's legendre(n, x)
 */
template<typename T>
inline T AssociatedLegendreP(
  Int n,
  Int m,
  T x)
{
  ASSERT(n >= 0); // As in Matlab we don't accept n<0
  ASSERT(m >= 0); // As in Matlab we don't accept m<0
  ASSERT(n >= m); // As in Matlab we don't accept n<m
  return boost::math::legendre_p<T>((int) n, (int) m, x);
}
  
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
template<typename T>
inline Complex<T> SphericalHarmonic(
  Int n,
  Int m,
  T theta,
  T phi)
{
  return boost::math::spherical_harmonic<T,T>((int) n, (int) m, theta, phi);
}
#endif
  
  
bool ElementaryOpTest();
  
  
} /**< namespace mcl  */
