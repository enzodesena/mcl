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
template<
  typename ConstForwardIterator,
  typename ForwardIterator,
  typename ForwardIteratorT = typename std::iterator_traits<ForwardIterator>::value_type,
  typename ConstForwardIteratorT = typename std::iterator_traits<ConstForwardIterator>::value_type>
void ForEach(
  ConstForwardIterator input_begin,
  const ConstForwardIterator input_end,
  std::function<ForwardIteratorT(ConstForwardIteratorT)> operation,
  ForwardIterator output_begin) noexcept
{
  while (input_begin != input_end)
  {
    *(output_begin++) = operation(*(input_begin++));
  }
}


template<
  typename ForwardIteratorInputA,
  typename ForwardIteratorInputB,
  typename ForwardIteratorOutput,
  typename T = typename std::iterator_traits<ForwardIteratorOutput>::value_type>
void ForEach(
  ForwardIteratorInputA input_a_begin,
  const ForwardIteratorInputA input_a_end,
  ForwardIteratorInputB input_b_begin,
  std::function<T(T,T)> pointwise_operation,
  ForwardIteratorOutput output_iter) noexcept
{
  while (input_a_begin != input_a_end)
  {
    *(output_iter++) = pointwise_operation(
      *(input_a_begin++), *(input_b_begin++));
  }
}

template<typename T>
void ForEach(
  const Vector<T>& input_vector,
  T (*operation)(T),
  Vector<T>& output_vector)
{
  ASSERT(input_vector.size() == output_vector.size());
  auto input_iter = input_vector.begin();
  auto output_iter = output_vector.begin();
  while (input_iter != input_vector.end())
  {
    *(output_iter++) = operation(*(input_iter++));
  }
}


template<typename T, typename U>
void ForEach(
  const Vector<T>& input_vector,
  std::function<U(
    T)> operation,
  Vector<U>& output_vector) noexcept
{
  ASSERT(input_vector.size() == output_vector.size());
  auto input_iter = input_vector.begin();
  auto output_iter = output_vector.begin();
  while (input_iter != input_vector.end())
  {
    *(output_iter++) = operation(*(input_iter++));
  }
}


template<typename T>
void ForEach(
  Vector<T>& vector,
  std::function<T(
    T)> operation) noexcept
{
  ForEach(vector, operation, vector);
}


template<typename T>
void ForEach(
  const Vector<T>& input_a,
  const Vector<T>& input_b,
  std::function<T(
    T,
    T)> pointwise_operation,
  Vector<T>& output) noexcept
{
  ASSERT(input_a.size() == input_b.size());
  ASSERT(input_a.size() == output.size());
  auto input_a_iter = input_a.begin();
  auto input_b_iter = input_b.begin();
  auto output_iter = output.begin();
  while (output_iter != output.end())
  {
    *(output_iter++) = pointwise_operation(*(input_a_iter++), *(input_b_iter++));
  }
}


template<
  typename ConstForwardIterator,
  typename ForwardIterator,
  typename T = typename std::iterator_traits<ForwardIterator>::value_type>
void Multiply(
  ConstForwardIterator input_begin,
  const ConstForwardIterator input_end,
  const T gain,
  ForwardIterator output_begin) noexcept
{
  std::function<T(T)> operation = [gain](T a) -> T { return  a * gain; };
  ForEach
  (
    input_begin,
    input_end,
    operation,
    output_begin);
}


template<
  typename ConstForwardIterator,
  typename ForwardIterator,
  typename T = typename std::iterator_traits<ForwardIterator>::value_type>
void Multiply(
  ConstForwardIterator input_a_begin,
  const ConstForwardIterator input_a_end,
  ConstForwardIterator input_b_begin,
  ForwardIterator output_begin) noexcept
{
  std::function<T(T,T)> operation = [](T a, T b) -> T { return  a * b; };
  ForEach
  (
    input_a_begin,
    input_a_end,
    input_b_begin,
    operation,
    output_begin);
}


template<
  typename ForwardIteratorInputA,
  typename ForwardIteratorInputB,
  typename ForwardIteratorOutput,
  typename T = typename std::iterator_traits<ForwardIteratorOutput>::value_type>
void MultiplyAdd(
  ForwardIteratorInputA input_to_multiply_begin,
  const ForwardIteratorInputA input_to_multiply_end,
  const T gain,
  ForwardIteratorInputB input_to_add_begin,
  ForwardIteratorOutput output_begin) noexcept
{
  std::function<T(T,T)> operation = [gain](T a, T b) -> T
  {
    return  a * gain + b;
  };
  ForEach
  (
    input_to_multiply_begin,
    input_to_multiply_end,
    input_to_add_begin,
    operation,
    output_begin);
}


template<
  typename ForwardIteratorInputA,
  typename ForwardIteratorInputB,
  typename ForwardIteratorOutput,
  typename T = typename std::iterator_traits<ForwardIteratorOutput>::value_type>
void Add(
  ForwardIteratorInputA input_a_begin,
  const ForwardIteratorInputA input_a_end,
  ForwardIteratorInputB input_b_begin,
  ForwardIteratorOutput output_begin) noexcept
{
  std::function<T(T,T)> operation = [](T a, T b) -> T { return  a + b; };
  ForEach
  (
    input_a_begin,
    input_a_end,
    input_b_begin,
    operation,
    output_begin);
}





/** Returns the opposite vector.Equivalent to Matlab's -vector. */
template<typename T>
Vector<T> Opposite(
  const Vector<T>& input) noexcept
{
  // Checking we are not dealing with unsigned types.
  static_assert(
    std::is_same<T,Complex<double>>::value ||
    std::is_same<T,Complex<float>>::value ||
    std::is_signed<T>::value, "");
  Vector<T> output(input.size());
  T (*operation)(
      T) = [](
    T value) -> T
  {
    return -value;
  };
  ForEach(input, operation, output);
  return std::move(output);
}


/** Returns the inverse vector.Equivalent to Matlab's 1./vector. */
template<typename T>
Vector<T> Inverse(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      T value) = [](
    T value)
  {
    return T(1.0) / value;
  };
  ForEach(input, operation, output);
  return std::move(output);
}


template<typename T>
Vector<T> HalfWave(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      T value) = [](
    T value)
  {
    return mcl::Max(T(0.0), value);
  };
  ForEach(input, operation, output);
  return output;
}


template<typename T>
Vector<T> Cos(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      T value) = [](
    T value)
  {
    return cos(value);
  };
  ForEach(input, operation, output);
  return output;
}


template<typename T>
Vector<T> Sin(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      T value) = [](
    T value)
  {
    return sin(value);
  };
  ForEach(input, operation, output);
  return output;
}


/**
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T>
Vector<T> Divide(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  ASSERT(vector_a.size() == vector_b.size());
  Vector<T> output((Int)vector_a.size());
  for (Int i = 0; i < (Int)vector_a.size(); ++i)
  {
    output[i] = vector_a[i] / vector_b[i];
  }
  return output;
}


//
/** Equivalent to Matlab's exp(vector). */
template<class T>
Vector<T> Exp(
  const Vector<T>& vector) noexcept
{
  Int n(vector.size());
  Vector<T> output(vector.size());
  for (Int i = 0; i < n; ++i)
  {
    output[i] = exp(vector[i]);
  }
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
  Vector<Complex<T>> output(input.size());
  Complex<T> (*operation)(
      Complex<T>) =
    [](
    Complex<T> value)
  {
    return Conj(value);
  };
  ForEach(input, operation, output);
  return std::move(output);
}


/** Equivalent to Matlab's real(input). */
template<typename T>
Vector<T> RealPart(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      Complex<T>) = [](
    Complex<T> value)
  {
    return value.real();
  };
  ForEach<Complex<T>,T>(input, operation, output);
  return std::move(output);
}


/** Equivalent to Matlab's imag(input). */
template<typename T>
Vector<T> Imag(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.size());
  ForEach(input, &std::imag, output);
  return std::move(output);
}


/**
 Returns the point-wise poser to exponent.
 Equivalent to Matlab's vector.^exponent
 */
template<typename T>
Vector<T> Pow(
  const Vector<T>& input,
  const T exponent) noexcept
{
  Vector<T> output(input.size());
  for (size_t i = 0; i < input.size(); ++i)
  {
    output[i] = Pow(input[i], exponent);
  }
  return std::move(output);
}


/** Equivalent to Matlab's abs(vector) */

inline Vector<double> Abs(
  const Vector<double>& input) noexcept
{
  Vector<double> output(input.size());
  ForEach(input, &std::fabs, output);
  return output;
}


/** Equivalent to Matlab's abs(vector) */

inline Vector<float> Abs(
  const Vector<float>& input) noexcept
{
  Vector<float> output(input.size());
  ForEach(input, &std::abs, output);
  return output;
}


/** Equivalent to Matlab's abs(vector) */
template<typename T>
Vector<double> Abs(
  const Vector<Complex<T>>& input) noexcept
{
  Vector<T> output(input.size());
  T (*operation)(
      Complex<T>) = [](
    Complex<T> value)
  {
    return mcl::Abs(value);
  };
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
Vector<T> Log(
  const Vector<T>& vector) noexcept
{
  Vector<T> output(vector.size());
  for (size_t i = 0; i < vector.size(); ++i)
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
Vector<T> Log10(
  const Vector<T>& vector) noexcept
{
  Vector<T> output(vector.size());
  for (size_t i = 0; i < vector.size(); ++i)
  {
    output[i] = log10(vector[i]);
  }
  return output;
}
} /**< namespace mcl */
