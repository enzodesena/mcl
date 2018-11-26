/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "elementaryop.hpp"
#include <functional>

namespace mcl
{


template<typename T, typename U>
inline void ForEach(
  const Vector<T>& input_vector,
  std::function<U(T)> operation,
  Vector<U>& output_vector) noexcept
{
  ASSERT(input_vector.length() == output_vector.length());
  auto input_iter = input_vector.begin();
  auto output_iter = output_vector.begin();
  while (input_iter != input_vector.end())
  {
    *(output_iter++) = operation(*(input_iter++));
  }
}

template<typename T>
inline void ForEach(
  Vector<T>& vector,
  std::function<T(T)> operation) noexcept
{
  ForEach(vector, operation, vector);
}


template<typename T>
inline void ForEach(
  const Vector<T>& input_a,
  const Vector<T>& input_b,
  std::function<T(T,T)> pointwise_operation,
  Vector<T>& output) noexcept
{
  ASSERT(input_a.length() == input_b.length());
  ASSERT(input_a.length() == output.length());
  auto input_a_iter = input_a.begin();
  auto input_b_iter = input_b.begin();
  auto output_iter = output.begin();
  while (output_iter != output.end())
  {
    *(output_iter++) = pointwise_operation(*(input_a_iter++), *(input_b_iter++));
  }
}




/** Fall-back multiply by a constant in case of no available optimisations. */
template<typename T>
inline void MultiplySerial(
  const Vector<T>& input,
  const T gain,
  Vector<T>& output) noexcept
{
  ASSERT(input.length() == output.length());
//  auto input_iter(input.begin());
//  auto output_iter(input.begin());
  for (size_t i = 0; i<output.length(); ++i)
  {
    output[i] = input[i] * gain;
//    *(output_iter++) = *(input_iter++) * gain;
  }
}

//template<typename T>
//inline void ForEach(
//  const Vector<T>& input_a,
//  const Vector<T>& input_b,
//  std::function<T(T,T)> pointwise_operation,
//  Vector<T>& output) noexcept

/** Fall-back multiply vectors in case of no available optimisations. */
template<typename T>
inline void MultiplySerial(
  const Vector<T>& input_a,
  const Vector<T>& input_b,
  Vector<T>& output) noexcept
{
  ForEach<T>(input_a, input_b, [] (T a, T b) { return a * b; }, output);
}

/**
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T>
inline Vector<T> MultiplySerial(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  ASSERT(vector_a.length() == vector_b.length());
  Vector<T> output(vector_a.length());
  MultiplySerial(vector_a, vector_b, output);
  
  return output;
}


template<typename T>
inline void MultiplyAddSerial(
  const Vector<T>& input_to_multiply,
  const T gain,
  const Vector<T>& input_to_add,
  Vector<T>& output) noexcept
{
  ASSERT(input_to_multiply.length() == input_to_add.length());
  ASSERT(input_to_add.length() == output.length());
  auto input_to_multiply_iter(input_to_multiply.begin());
  auto input_to_add_iter(input_to_add.begin());
  auto output_iter(output.begin());
  while (output_iter != output.end())
  {
    *(output_iter++) =
      *(input_to_multiply_iter++) * gain + *(input_to_add_iter++);
  }
}




/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<typename T>
inline void AddSerial(
  const Vector<T>& input_a,
  const Vector<T>& input_b,
  Vector<T>& output) noexcept
{
  ASSERT(input_a.length() == input_a.length());
  ASSERT(input_a.length() == output.length());
  auto input_a_iter = input_a.begin();
  auto input_b_iter = input_b.begin();
  auto output_iter = output.begin();
  while (output_iter != output.end())
  {
    *(output_iter++) = *(input_a_iter++) + *(input_b_iter++);
  }
}


  
template<typename T>
inline void ForEach(
  const Vector<T>& input_vector,
  T (*operation)(T),
  Vector<T>& output_vector)
{
  ASSERT(input_vector.length() == output_vector.length());
  auto input_iter = input_vector.begin();
  auto output_iter = output_vector.begin();
  while (input_iter != input_vector.end())
  {
    *(output_iter++) = operation(*(input_iter++));
  }
}
  
/** Returns the opposite vector.Equivalent to Matlab's -vector. */
template<typename T>
inline Vector<T> Opposite(
  const Vector<T>& input) noexcept
{
  // Checking we are not dealing with unsigned types.
  static_assert(
    std::is_same<T, Complex<double>>::value ||
    std::is_same<T, Complex<float>>::value ||
    std::is_signed<T>::value, "");
  Vector<T> output(input.length());
  T (*operation)(T) = [] (T value) -> T { return -value; };
  ForEach(input, operation, output);
  return std::move(output);
}

  
/** Returns the inverse vector.Equivalent to Matlab's 1./vector. */
template<typename T>
inline Vector<T> Inverse(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (T value) = [] (T value) { return T(1.0) / value; };
  ForEach(input, operation, output);
  return std::move(output);
}

template<typename T>
inline Vector<T> HalfWave(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (T value) = [] (T value) { return mcl::Max(T(0.0), value); };
  ForEach(input, operation, output);
  return output;
}

template<typename T>
inline Vector<T> Cos(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (T value) = [] (T value) { return cos(value); };
  ForEach(input, operation, output);
  return output;
}

template<typename T>
inline Vector<T> Sin(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (T value) = [] (T value) { return sin(value); };
  ForEach(input, operation, output);
  return output;
}


  
/**
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T>
Vector<T> Divide(const Vector<T>& vector_a,
                      const Vector<T>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  Vector<T> output((Int)vector_a.length());
  for (Int i=0; i<(Int)vector_a.length(); ++i) {
    output[i] = vector_a[i]/vector_b[i];
  }
  return output;
}

//
/** Equivalent to Matlab's exp(vector). */
template<class T>
Vector<T> Exp(const Vector<T>& vector) noexcept {
  Int n(vector.length());
  Vector<T> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = exp(vector[i]); }
  return output;
}

/**
 Returns the vector with conjugate entries.
 Equivalent to Matlab's conj(vector).
 */
template<typename T>
Vector<Complex<T>> Conj(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<Complex<T>> output(input.length());
  Complex<T> (*operation) (Complex<T>) =
    [] (Complex<T> value) { return Conj(value); };
  ForEach(input, operation, output);
  return std::move(output);
}

/** Equivalent to Matlab's real(input). */
template<typename T>
Vector<T> RealPart(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (Complex<T>) = [] (Complex<T> value) { return value.real(); };
  ForEach<Complex<T>,T>(input, operation, output);
  return std::move(output);
}

/** Equivalent to Matlab's imag(input). */
template<typename T>
Vector<T> Imag(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.length());
  ForEach(input, &std::imag, output);
  return std::move(output);
}


/**
 Returns the point-wise poser to exponent.
 Equivalent to Matlab's vector.^exponent
 */
template<typename T>
inline Vector<T> Pow(
  const Vector<T>& input,
  const T exponent) noexcept
{
  Vector<T> output(input.length());
  for (size_t i=0; i<input.length(); ++i)
  {
    output[i] = Pow(input[i], exponent);
  }
  return std::move(output);
}


/** Equivalent to Matlab's abs(vector) */

inline Vector<double> Abs(
  const Vector<double>& input) noexcept
{
  Vector<double> output(input.length());
  ForEach(input, &std::fabs, output);
  return output;
}

/** Equivalent to Matlab's abs(vector) */

inline Vector<float> Abs(
  const Vector<float>& input) noexcept
{
  Vector<float> output(input.length());
  ForEach(input, &std::abs, output);
  return output;
}

/** Equivalent to Matlab's abs(vector) */
template<typename T>
inline Vector<double> Abs(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.length());
  T (*operation) (Complex<T>) = [] (Complex<T> value) { return mcl::Abs(value); };
  ForEach<Complex<T>,T>(input, operation, output);
  return std::move(output);
}
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
/**
 Returns the natural logarithm of the elements of vector.
 Equivalent to Matlab's log(vector).
 */
template<typename T>
inline Vector<T> Log(
  const Vector<T>& vector) noexcept
{
  Vector<T> output(vector.length());
  for (size_t i=0; i<vector.length(); ++i)
  {
    output[i] = log(vector[i]);
  }
  return output;
}
  
  
/**
 Returns the 10-base logarithm of the elements of vector.
 Equivalent to Matlab's log10(vector).
 */
template<typename T>
inline Vector<T> Log10(
  const Vector<T>& vector) noexcept
{
  Vector<T> output(vector.length());
  for (size_t i=0; i<vector.length(); ++i)
  {
    output[i] = log10(vector[i]);
  }
  return output;
}

} /**< namespace mcl */

