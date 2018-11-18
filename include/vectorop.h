/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once


//#include "mcltypes.h"
#include "elementaryop.h"
//#include "basicop.h"
//#include "matrixop.h"
#include "vector.h"
//#include <vector>
//#include <iostream>

namespace mcl {

  
  
/** Equivalent to Matlab's length(input). */
template<class T, int length>
int Length(const Vector<T,length>& input) noexcept {
  return input.length();
}


template <class T, int length>
void SetToZero(Vector<T, length>& vector)
{
  for (auto& element : vector)
  {
    element = 0.0;
  }
}


/** Returns a vector of zeros */
template <class T, int length>
Vector<T, length> Zeros() noexcept {
  Vector<T, length> vector(length);
  SetToZero(vector);
  return std::move(vector);
}


template <class T>
Vector<T, kDynamicLength> Zeros(int length) noexcept {
  Vector<T, kDynamicLength> vector(length);
  SetToZero(vector);
  return std::move(vector);
}

template <class T>
  Vector<T, kDynamicLength> EmptyVector() noexcept {
  return Vector<T, kDynamicLength>(0);
}
//
//
//template<typename T, int length>
//void Multiply(
//    const Vector<T, length>& input,
//    const T gain,
//    Vector<T, length>& output) noexcept {
//  ASSERT(input.length() == output.length());
//  for (Int i=0; i<output.length(); ++i) {
//    output.At(i) = input.At(i)*gain;
//  }
//}
//
///**
// Returns the point by point multiplication of the vector with the gain.
// Equivalent to Matlab's vector_a.*gain.
// */
//template<typename T, int length>
//Vector<T, length> Multiply(
//    const Vector<T, length>& vector,
//    const T gain) noexcept {
//  Vector<T,length> output;
//  Multiply(vector, gain, output);
//  return std::move(output);
//}
//
//
//void Multiply(
//    const double* input_data,
//    const int num_samples,
//    const double gain,
//    double* output_data) noexcept;
//
//template<int length>
//Vector<double, length> Multiply(
//    const Vector<double, length>& input,
//    const double gain) noexcept {
//  Vector<double,length> output;
//  Multiply(input.data(), input.length(), gain, output.data());
//  return std::move(output);
//}
//
///** This calculates the multiplication of a vector (`input_data_mult`)
// by a constant (`gain`), and then adds the resulting vector to
// another vector (`input_data_add'). */
//void MultiplyAdd(
//    const double* input_data_mult,
//    const double gain,
//    const double* input_data_add,
//    const int num_samples,
//    double* output_data) noexcept;
//
//
///**
// Returns the point by point addition of the two vectors.
// Equivalent to Matlab's vector_a+vector_b.
// */
//template<class T, int length>
//void Add(
//    const Vector<T, length>& vector,
//    const T scalar,
//    Vector<T, length>& output_vector) noexcept {
//  ASSERT(vector.length() == output_vector.length());
//  for (int i=0; i<vector.length(); ++i) { output_vector[i] = vector[i]+scalar; }
//}
//
//template<class T, int length>
//Vector<T, length> Add(
//    const Vector<T, length>& vector,
//    const T scalar) noexcept {
//  Vector<T, length> output(vector.length());
//  Add(vector, scalar, output);
//  return std::move(output);
//}
//
//template<class T>
//T Add(const T* input_data, const Int num_samples) noexcept;
//
//
/**
 Returns the subset of elements with indexes from_index and to_index.
 Equivalent to Matlab's vector(from_index:to_index). (Careful about the
 different indexes convention between here and Matlab.
 */
template<class T>
Vector<T,kDynamicLength> Subset(
    const Vector<T>& input,
    const int from_index,
    const int to_index) noexcept {
  ASSERT(from_index < input.length());
  ASSERT(to_index < input.length());
  ASSERT(from_index <= to_index);
  Vector<T,kDynamicLength> output(to_index-from_index+1);
  for (int i = from_index; i<=to_index; ++i) {
    output[i-from_index] = input[i];
  }
  return std::move(output);
}


/**
 Returns the concatenation of vector_a and vector_b. Equivalent to Matlab's
 [vector_a; vector_b].
 */
template<class T>
Vector<T,kDynamicLength> Concatenate(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  Vector<T,kDynamicLength> output(vector_a.length() + vector_b.length());
  auto output_iter = output.begin();
  for (auto& element : vector_a)
  {
    *(output_iter++) = element;
  }
  for (auto& element : vector_b)
  {
    *(output_iter++) = element;
  }
  return std::move(output);
}


/** Returns a vector with only one element. */
template<class T>
Vector<T,kDynamicLength> UnaryVector(
  const T& element) noexcept
{
  Vector<T,kDynamicLength> output(1);
  output[0] = element;
  return std::move(output);
}

//
///** Returns a vector with two elements. */
//template<class T>
//Vector<T> BinaryVector(const T& element_a, const T& element_b) noexcept {
//  Vector<T> output(2);
//  output[0] = element_a;
//  output[1] = element_b;
//  return output;
//}
//
/**
 Flips the vector. Equivalent to matlab's flipud or fliplr (which for vectors
 are equivalent).
 */
template<class T, int length>
Vector<T,length> Flip(
  const Vector<T,length>& input) noexcept
{
  if (input.length() <= 1) { return input; }
  Vector<T,length> output(input.length());
  auto input_iter = input.end()-1; // Start from last element of vector
  auto output_iter = output.begin();
  while (output_iter != output.end()) {
    *(output_iter++) = *(input_iter--);
  }
  return std::move(output);
}

/**
 Equivalent to Matlab's circshift(vector, num_positions). A positive
 num_positions corresponds to a forward shift.
 */
template<class T, int length>
Vector<T,length> CircShift(
  const Vector<T,length>& vector,
  int num_positions) noexcept
{
  int N(vector.length());
  Vector<T> output(N);
  for (int i=0; i<N; ++i) {
    int index = static_cast<int>(Mod(i - num_positions, N));
    output[i] = vector[index];
  }
  return output;
}

///** Equivalent to Matlab's conv(vector_a, vector_b). */
//template<class T>
//Vector<T> Conv(const Vector<T>& vector_a,
//                    const Vector<T>& vector_b) noexcept {
//  Int N_a = (Int)vector_a.length();
//  Int N_b = (Int)vector_b.length();
//  Int out_length = N_a+N_b-1;
//
//  Vector<T> moving_vector_temp = Concatenate(Zeros<T>(N_b-1),
//                                                  Flip(vector_a));
//  Vector<T> moving_vector_a = Concatenate(moving_vector_temp,
//                                               Zeros<T>(N_b-1));
//
//  Vector<T> output = Zeros<T>(out_length);
//  for (Int n=0; n<out_length; ++n) {
//    for (Int m=0; m<N_b; ++m) {
//      output[out_length-n-1] += moving_vector_a[n+m]*vector_b[m];
//    }
//  }
//  return output;
//}
//
//
///**
// Adds all the vectors and zero-pads short vectors if they have different
// lengths.
// */
//template<class T>
//Vector<T>
//AddVectors(const Vector<Vector<T> >& vectors) noexcept {
//  // Get maximum length
//  Vector<Int> vector_lengths(vectors.length());
//  for (Int i=0; i<(Int)vectors.length(); ++i) {
//    vector_lengths[i] = (Int)vectors[i].length();
//  }
//  Int max_length(Max(vector_lengths));
//
//  Vector<T> output = Zeros<T>(max_length);
//  for (Int i=0; i<(Int)vectors.length(); ++i) {
//    output = Add(output, ZeroPad(vectors[i], max_length));
//  }
//
//  return output;
//}
//
///**
// Adds two vectors and zero-pads the shorter one if they have different
// lengths.
// */
//template<class T>
//Vector<T> AddVectors(const Vector<T>& vector_a,
//                          const Vector<T>& vector_b) noexcept {
//  // Get maximum length
//  Int max_length(Max((Int)vector_a.length(), (Int)vector_b.length()));
//
//  Vector<T> output = Zeros<T>(max_length);
//  output = Add(output, ZeroPad(vector_a, max_length));
//  output = Add(output, ZeroPad(vector_b, max_length));
//
//  return output;
//}
//
//
///** Interleaves two vectors, with the first element of `vector_a` going
// first.*/
//template<class T>
//Vector<T> Interleave(const Vector<T>& vector_a,
//                          const Vector<T>& vector_b) noexcept {
//  ASSERT(vector_a.length() == vector_b.length());
//
//  Vector<T> output;
//  for (Int i=0; i<(Int)vector_a.length(); ++i) {
//    output.push_back(vector_a[i]);
//    output.push_back(vector_b[i]);
//  }
//  return output;
//}
//
/** Decreases the sampling frequency of the input vector by keeping
 the first sample and then every `downsampling_factor`-th sample after the first. */
template<class T>
Vector<T> Downsample(
  const Vector<T>& vector,
  const int downsampling_factor) noexcept
{
  ASSERT(downsampling_factor >= 1);
  Vector<T> output;
  for (int i=0; i<vector.length(); i += downsampling_factor)
  {
    output.PushBack(vector[i]);
  }
  return output;
}

///**
// This is equivalent to Matlab's from:to. E.g. 3:5 outputs a vector [3,4,5].
// TODO: Implement fractional input.
// */
//template<class T>
//Vector<T> ColonOperator(const Int from, const Int to) noexcept {
//  if ((to-from) < 0) { return EmptyVector<T>(); }
//  const Int vector_length = (UInt) (to-from+1);
//  Vector<T> output(vector_length);
//  for (Int i=0; i<vector_length; ++i) {
//    output[i] = ((T) i) + ((T) from);
//  }
//  return output;
//}
//
///**
// This is equivalent to Matlab's from:step:to. E.g. 3:2:6 outputs a vector
// [3,4,6].
// */
//// TODO: implement negative step and fractional input.
//Vector<Real>
//ColonOperator(const Real from, const Real step, const Real to) noexcept;
//
//
///**
// Returns elements of vector `vector` from from_id to to_id
// (including extremes).
// */
//template<class T>
//Vector<T> Elements(const Vector<T>& vector,
//                                const Int from_id,
//                                const Int to_id) noexcept {
//  return Vector<T>(vector.begin() + ((Int)from_id),
//                        vector.begin() + ((Int)to_id)+1);
//}
//
//
//template<class T>
//Vector<T> GetSegment(const Vector<T>& vector,
//                                  const Int subset_id,
//                                  const Int subset_length,
//                                  bool zeropad_if_shorter = false) noexcept {
//  const Int size = vector.length();
//
//  const Int from_sample = subset_id * subset_length;
//  if (from_sample >= size) {
//    if (zeropad_if_shorter) {
//      return Zeros<T>(subset_length);
//    } else {
//      return Vector<T>(); // Return empty vector
//    }
//  }
//
//  const Int to_sample = Min(from_sample + subset_length - 1,
//                      size - 1);
//
//  const Int actual_length = to_sample - from_sample + 1;
//  if (zeropad_if_shorter && actual_length < subset_length) {
//    return ZeroPad(Elements(vector, from_sample, to_sample), subset_length);
//  } else {
//    return Elements(vector, from_sample, to_sample);
//  }
//}
//
///**
// Multiplies all the elements in the vector. Equivalent to Matlab's
// prod(vector).
// */
//template<class T>
//T Prod(const Vector<T>& vector) noexcept {
//  const Int num_elements = vector.length();
//  T output = (T) 1.0;
//  for (Int i=0; i<num_elements; ++i) { output *= vector[i]; }
//  return output;
//}
//
///** Dot product between two vectors. Equivalent to Matlab's dot(a,b) */
//template<class T>
//T Dot(
//    const Vector<T>& vector_a,
//    const Vector<T>& vector_b) noexcept {
//  const Int num_elements = (Int)vector_a.length();
//  ASSERT(num_elements == (Int)vector_b.length());
//
//  T output = (T) 0.0;
//  for (Int i=0; i<num_elements; ++i) {
//    output += vector_a[i]*vector_b[i];
//  }
//  return output;
//}
//
//Real Norm(const Vector<Real>& vector, Real l_norm = 2.0) noexcept;
//
template<typename T, int length>
void Print(const Vector<T,length>& vector) noexcept
{
  std::cout<<"\n------------\n";
  for (auto iter = vector.begin(); iter != vector.end(); ++iter)
  {
    std::cout<<*iter<<std::endl;
  }
  std::cout<<"------------\n";
}

///** Returns a real vector of `length` ones. */
//Vector<Real> Ones(Int length) noexcept;
//
//
//Vector<Real> Hann(const Int length) noexcept;
//
///** Returns a Hamming window of length `length' */
//Vector<Real> Hamming(const Int length) noexcept;
//
//Vector<Real> TukeyWin(const Int length,
//                                   const Real ratio) noexcept;
//
//

/** Equivalent to Matlab's linspace(min, max, num_elements); */
template<typename T>
Vector<Real,kDynamicLength> LinSpace(
  T min,
  T max,
  int num_elements) noexcept
{
  if (num_elements <= 1) { return UnaryVector(max); }
  T interval = max-min;
  T slot = interval / ((T) (num_elements-1));
  Vector<Real,kDynamicLength> output(num_elements);
  for (int i=0; i<num_elements; ++i)
  {
    output[i] = min + slot*((T) i);
  }
  return output;
}


template <typename T, int length>
T Sum(const Vector<T,length>& input) noexcept
{
  T output((T) 0.0);
  for (auto iter = input.begin(); iter != input.end(); iter++)
  {
    output += *iter;
  }
  return output;
}


template <typename T, int length>
/** Equivalent to Matlab's mean(input) */
Real Mean(const Vector<T,length>& input) noexcept
{
  return Sum(input) / ((T) input.length());
}

///**
// Returns the geometric mean of the input vector. Equivalent
// to Matlab's geomean(input)
// **/
//Real Geomean(const Vector<Real>& input) noexcept;
//
///**
// Weighted mean. Not implemented in Matlab (but should be). The weights are
// normalised inside the function. Hence Mean(input, ones(N)) gives the same
// result as Mean(input, ones(N)/N).
// */
//Real Mean(
//    const Vector<Real>& input,
//    const Vector<Real>& weigths) noexcept;
//
///**
// Returns the standard deviation of the `input` vector. Equivalent to Matlab's
// std(input). This includes the correction for having an unbiased estimator.
// */
//Real Std(const Vector<Real>& input) noexcept;
//
///** Var (unbiased estimator) */
//Real Var(const Vector<Real>& input) noexcept;
//
///** Weighted var (biased estimator) */
//Real Var(const Vector<Real>& input,
//                 const Vector<Real>& weights) noexcept;
//
///** Splits a string using a delimiter. */
//Vector<std::string> Split(
//    const std::string& string,
//    char delim) noexcept;
//
///** Converts roots to polynomial. Equivalent to Matlab's poly(roots) */
//Vector<Complex> Poly(const Vector<Complex> roots) noexcept;
//Vector<Complex> Poly(const Vector<Real> roots) noexcept;
//
///** Returns true if all elements are non negative */
//bool IsNonNegative(const Vector<Real>& input) noexcept;
//
//Matrix<Real> Cov(const Vector<Real>& x,
//                         const Vector<Real>& y) noexcept;
//
//Matrix<Real> Cov(const Vector<Vector<Real> >& input) noexcept;
//
//Real CovElement(const Vector<Real>& x,
//                        const Vector<Real>& y) noexcept;
//
///**
// Returns a vector containing the cumulative sum of
// the elements of X. Equivalent to Matlab's cumsum(input)
// */
//Vector<Real> CumSum(const Vector<Real>& input) noexcept;
//
///** Splits signal up into (overlapping) frames */
//Vector<Vector<Real> > Enframe(const Vector<Real>& input,
//                                                const Vector<Real>& window,
//                                                const Int frame_increment) noexcept;
//
//Vector<Real> OverlapAdd(const Vector<Vector<Real> >& frames,
//                                     const Vector<Real>& window,
//                                     const Int frame_increment) noexcept;
//
//Vector<Complex> ConvertToComplex(Vector<Real> input) noexcept;
//
//

/**
 Adds zero until the output vector has a length of total_length.
 If the length of input is smaller than total_length, than it returns the
 vector with the first total_length elements.
 */
template<typename T, int input_length, int output_length>
void ZeroPad(
  const Vector<T,input_length>& input,
  Vector<T,output_length>& output) noexcept
{
  auto input_iter = input.begin();
  auto output_iter = output.begin();
  while (output_iter != output.end())
  {
    if (input_iter < input.end())
    {
      *(output_iter++) = *(input_iter++);
    }
    else
    {
      *(output_iter++) = static_cast<T>(0.0);
    }
  }
}


} /**< namespace mcl */
