/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once


//#include <cassert>
//#include "mcltypes.h"
#include "vector.h"
//#include <limits>

namespace mcl {

  
//template<typename T, typename U>
//void ForEach(
//  Vector<T>& vector,
//  T (*operation)(U))
//{
//  for (auto& element : vector) {
//    element = operation(element);
//  }
//}
//
//template<typename T, typename U>
//void ForEach(
//  const Vector<T>& input_vector,
//  T (*operation)(U),
//  Vector<T>& output_vector)
//{
//  ASSERT(input_vector.length() == output_vector.length());
//  auto input_iter = input_vector.begin();
//  auto output_iter = output_vector.begin();
//  while (input_iter != input_vector.end())
//  {
//    *output_iter = operation(*input_iter);
//  }
//}

/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<class T>
inline Vector<T> Add(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  ASSERT(vector_a.length() == vector_b.length());
  
  Vector<T> output((Int)vector_a.length());
  for (Int i=0; i<(Int)vector_a.length(); ++i)
  {
    output[i] = vector_a[i]+vector_b[i];
  }
  return output;
}
  
////void Add(
////  const Real* input_data_a,
////  const Real* input_data_b,
////  const size_t num_samples,
////  Real* output_data) noexcept;
//
//
////template<>
////inline Vector<Real> Add(
////  const Vector<Real>& vector_a,
////  const Vector<Real>& vector_b) noexcept {
////  ASSERT(vector_a.length() == vector_b.length());
////
////  Vector<Real> output((Int)vector_a.length());
////  Add(vector_a.data(), vector_b.data(), (Int)vector_a.length(),
////      output.data());
////  return output;
////}
//  
/** Returns the opposite vector.Equivalent to Matlab's -vector. */
template<class T>
Vector<T> Opposite(
  const Vector<T>& vector) noexcept
{
  // Checking we are not dealing with unsigned types.
  // The assert below responds false to complex. TODO: fix this
//  ASSERT(std::numeric_limits<T>::is_signed);
  Vector<T> output(vector.length());
  for (size_t i=0; i<vector.length(); ++i)
  {
    output[i] = -vector[i];
  }
  return output;
}
//
//  
///** Returns the inverse vector.Equivalent to Matlab's 1./vector. */
//Vector<Real> Inverse(const Vector<Real>& vector) noexcept;
//
///** 
// Returns the point by point subtraction of the two vectors.
// Equivalent to Matlab's vector_a-vector_b.
// */
//template<class T> 
//Vector<T> Subtract(const Vector<T>& vector_a,
//                        const Vector<T>& vector_b) noexcept {
//  return Add(vector_a, Opposite(vector_b));
//}
//
//
/**
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T>
Vector<T> Multiply(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  ASSERT(vector_a.length() == vector_b.length());
  Vector<T> output(vector_a.length());
  for (size_t i=0; i<vector_a.length(); ++i)
  {
    output[i] = vector_a[i]*vector_b[i];
  }
  return output;
}
//
//void Multiply(const Real* input_data_a, const Real* input_data_b,
//              Int num_samples, Real* output_data) noexcept;
//  
//  
//template<>
//inline Vector<Real> Multiply(const Vector<Real>& vector_a,
//                                  const Vector<Real>& vector_b) noexcept {
//  ASSERT(vector_a.length() == vector_b.length());
//  
//  Vector<Real> output((Int)vector_a.length());
//  Multiply(vector_a.data(), vector_b.data(), (Int)vector_a.length(),
//           output.data());
//  return output;
//}
//  
///** 
// Returns the point by point multiplication of the two vectors.
// Equivalent to Matlab's vector_a.*vector_b.
// */
//template<class T>
//Vector<T> Divide(const Vector<T>& vector_a,
//                      const Vector<T>& vector_b) noexcept {
//  ASSERT(vector_a.length() == vector_b.length());
//  Vector<T> output((Int)vector_a.length());
//  for (Int i=0; i<(Int)vector_a.length(); ++i) {
//    output[i] = vector_a[i]/vector_b[i];
//  }
//  return output;
//}
//
//
///** Equivalent to Matlab's exp(vector). */
//template<class T>
//Vector<T> Exp(const Vector<T>& vector) noexcept {
//  Int n(vector.length());
//  Vector<T> output(vector.length());
//  for (Int i=0; i<n; ++i) { output[i] = exp(vector[i]); }
//  return output;
//}
//  
///** 
// Returns the vector with conjugate entries.
// Equivalent to Matlab's conj(vector).
// */
//template<typename T>
//Vector<Complex<T>> Conj(
//  const Vector<Complex<T>>& input) noexcept
//{
//  Vector<Complex<T>> output(input.length());
//  ForEach(input, &Conj, output);
//  return std::move(output);
//}
//
/** Transform real vector into complex vector with null imaginary part */
template<typename T>
Vector<Complex<T>>
ComplexVector(const Vector<T>& input) noexcept
{
  Vector<Complex<T>> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i)
  {
    output[i] = Complex<T>(input[i], 0.0);
  }
  return std::move(output);
}
//
///** Equivalent to Matlab's real(input). */
//template<typename T>
//Vector<T> RealPart(
//  const Vector<Complex<T>>& input) noexcept
//{
//  Vector<T> output(input.length());
//  ForEach(input, &std::real, output);
//  return std::move(output);
//}
//
///** Equivalent to Matlab's imag(input). */
//template<typename T>
//Vector<T> Imag(
//  const Vector<Complex<T>>& input) noexcept
//{
//  Vector<T> output(input.length());
//  ForEach(input, &std::imag, output);
//  return std::move(output);
//}
//
//
///** 
// Returns the point-wise poser to exponent.
// Equivalent to Matlab's vector.^exponent
// */
//template<typename T>
//Vector<T> Pow(
//  const Vector<T>& input,
//  T exponent) noexcept
//{
//  Vector<T> output(input.length());
//  for (size_t i=0; i<input.length(); ++i)
//  {
//    output[i] = Pow(input[i], exponent);
//  }
//  return std::move(output);
//}
//
//
///** Equivalent to Matlab's abs(vector) */
//template<typename T>
//Vector<T> Abs(
//  const Vector<T>& input) noexcept
//{
//  Vector<T> output(input.length());
//  ForEach(input, &std::fabs, output);
//  return output;
//}
//
///** Equivalent to Matlab's abs(vector) */
//template<typename T>
//Vector<double> Abs(
//  const Vector<Complex<T>>& input) noexcept
//{
//  Vector<T> output(input.length());
//  ForEach(input, &std::abs, output);
//  return std::move(output);
//}
//
///** Equivalent to Matlab's vector.*(vector>0) */
//Vector<Real> HalfWave(const Vector<Real>& vector) noexcept;
//
///** Equivalent to Matlab's cos(vector) */
//Vector<Real> Cos(const Vector<Real>& vector) noexcept;
//
///** Equivalent to Matlab's sin(vector) */
//Vector<Real> Sin(const Vector<Real>& vector) noexcept;
//  
///** 
// Returns the natural logarithm of the elements of vector.
// Equivalent to Matlab's log(vector).
// */
//Vector<Real> Log(const Vector<Real>& vector) noexcept;
//  
///**
// Returns the 10-base logarithm of the elements of vector.
// Equivalent to Matlab's log10(vector).
// */
//Vector<Real> Log10(const Vector<Real>& vector) noexcept;
//
template<typename TOrigin, typename TDestination>
Vector<TDestination> Convert(
  const Vector<TOrigin>& vector) noexcept
{
  const size_t length = vector.length();
  Vector<TDestination> output(length);
  for (size_t i=0; i<length; ++i)
  {
    output[i] = static_cast<TDestination>(vector[i]);
  }
  return output;
}

bool PointWiseOpTest();

} /**< namespace mcl */

