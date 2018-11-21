/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once


#include "elementaryop.h"
#include "basicop.h"
#include "mclintrinsics.h"
#include "vector.h"

namespace mcl {

template<typename T>
inline T Mod(
  const T x,
  const T y);
  
template<class T>
inline T Min(
  const T scalar_a,
  const T scalar_b) noexcept;
  
template<typename T>
bool IsNonNegative(
  const Vector<T>& vector) noexcept
{
  for (auto& element : vector)
  {
    if (element < 0.0)
    {
      return false;
    }
  }
  return true;
}
  
/** Equivalent to Matlab's length(input). */
template<class T, size_t length>
size_t Length(const Vector<T,length>& input) noexcept {
  return input.length();
}


template <class T, size_t length>
void SetToZero(Vector<T, length>& vector)
{
  for (auto& element : vector)
  {
    element = 0.0;
  }
}

/**
 Adds zero until the output vector has a length of total_length.
 If the length of input is smaller than total_length, than it returns the
 vector with the first total_length elements.
 */
template<typename T, size_t input_length, size_t output_length>
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

/** Returns a vector of zeros */
template <class T, size_t length>
Vector<T, length> Zeros() noexcept
{
  Vector<T, length> vector(length);
  SetToZero(vector);
  return std::move(vector);
}


template <class T>
Vector<T, kDynamicLength> Zeros(
  size_t length) noexcept
{
  Vector<T, kDynamicLength> vector(length);
  SetToZero(vector);
  return std::move(vector);
}

template <class T>
  Vector<T, kDynamicLength> EmptyVector() noexcept
{
  return Vector<T, kDynamicLength>(0);
}




  
template<typename T, size_t length>
void Multiply(
  const Vector<T>& input,
  const T gain,
  Vector<T>& output) noexcept
{
  MathIntrinsics<T,length>::Multiply(
    input,
    gain,
    output);
}

template<typename T, size_t length>
void MultiplyAdd(
  const Vector<T,length>& input_to_multiply,
  const T gain,
  const Vector<T,length>& input_to_add,
  Vector<T,length>& output) noexcept
{
  MathIntrinsics<T,length>::MultiplyAdd(
    input_to_multiply,
    gain,
    input_to_add,
    output);
}




/**
 Returns the point by point multiplication of the vector with the gain.
 Equivalent to Matlab's vector_a.*gain.
 */
template<typename T, size_t length>
Vector<T, length> Multiply(
  const Vector<T,length>& input,
  const T gain) noexcept
{
  Vector<T,length> output(input.length());
  MathIntrinsics<T,length>::Multiply(
    input,
    gain,
    output);
  return std::move(output);
}
  
  


/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<class T, size_t length>
void Add(
  const Vector<T, length>& vector,
  const T scalar,
  Vector<T, length>& output_vector) noexcept
{
  ASSERT(vector.length() == output_vector.length());
  for (size_t i=0; i<vector.length(); ++i)
  {
    output_vector[i] = vector[i]+scalar;
  }
}

template<class T, size_t length>
Vector<T, length> Add(
  const Vector<T, length>& vector,
  const T scalar) noexcept
{
  Vector<T, length> output(vector.length());
  Add(vector, scalar, output);
  return std::move(output);
}

/**
 Returns the subset of elements with indexes from_index and to_index.
 Equivalent to Matlab's vector(from_index:to_index). (Careful about the
 different indexes convention between here and Matlab.
 */
template<class T>
Vector<T> Subset(
    const Vector<T>& input,
    const size_t from_index,
    const size_t to_index) noexcept {
  ASSERT(from_index < input.length());
  ASSERT(to_index < input.length());
  ASSERT(from_index <= to_index);
  Vector<T> output(to_index-from_index+1);
  for (size_t i = from_index; i<=to_index; ++i) {
    output[i-from_index] = input[i];
  }
  return std::move(output);
}


/**
 Returns the concatenation of vector_a and vector_b. Equivalent to Matlab's
 [vector_a; vector_b].
 */
template<class T>
Vector<T> Concatenate(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  Vector<T> output(vector_a.length() + vector_b.length());
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
Vector<T> UnaryVector(
  const T& element) noexcept
{
  Vector<T> output(1);
  output[0] = element;
  return std::move(output);
}


/** Returns a vector with two elements. */
template<class T>
Vector<T> BinaryVector(
  const T& element_a,
  const T& element_b) noexcept
{
  Vector<T> output(2);
  output[0] = element_a;
  output[1] = element_b;
  return output;
}

/**
 Flips the vector. Equivalent to matlab's flipud or fliplr (which for vectors
 are equivalent).
 */
template<class T, size_t length>
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
template<class T, size_t length>
Vector<T,length> CircShift(
  const Vector<T,length>& vector,
  Int num_positions) noexcept
{
  Int N(static_cast<Int>(vector.length()));
  Vector<T> output(N);
  for (Int i=0; i<N; ++i)
  {
    Int index = static_cast<Int>(Mod(i - num_positions, N));
    output[i] = vector[index];
  }
  return output;
}

/** Equivalent to Matlab's conv(vector_a, vector_b). */
template<class T>
Vector<T> Conv(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  size_t N_a = vector_a.length();
  size_t N_b = vector_b.length();
  size_t out_length = N_a+N_b-1;

  Vector<T> moving_vector_temp = Concatenate(Zeros<T>(N_b-1), Flip(vector_a));
  Vector<T> moving_vector_a = Concatenate(moving_vector_temp, Zeros<T>(N_b-1));

  Vector<T> output = Zeros<T>(out_length);
  for (size_t n=0; n<out_length; ++n)
  {
    for (size_t m=0; m<N_b; ++m)
    {
      output[out_length-n-1] += moving_vector_a[n+m]*vector_b[m];
    }
  }
  return output;
}

/**
 Adds all the vectors and zero-pads short vectors if they have different
 lengths.
 */
template<class T>
Vector<T>
AddVectors(
  const Vector<Vector<T>>& vectors) noexcept
{
  // Get maximum length
  Vector<size_t> vector_lengths(vectors.length());
  for (size_t i=0; i<vectors.length(); ++i)
  {
    vector_lengths[i] = vectors[i].length();
  }
  size_t max_length(Max(vector_lengths));
  Vector<T> output = Zeros<T>(max_length);
  Vector<T> temp(max_length);
  for (auto& vector : vectors)
  {
    if (vector.length() == max_length)
    {
      output = Add(output, vector);
    }
    else
    {
      ZeroPad(vector, temp);
      output = Add(output, temp);
    }
  }
  return std::move(output);
}
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
/** Interleaves two vectors, with the first element of `vector_a` going
 first.*/
template<class T>
Vector<T> Interleave(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  ASSERT(vector_a.length() == vector_b.length());

  Vector<T> output;
  for (Int i=0; i<(Int)vector_a.length(); ++i)
  {
    output.push_back(vector_a[i]);
    output.push_back(vector_b[i]);
  }
  return output;
}
//
/** Decreases the sampling frequency of the input vector by keeping
 the first sample and then every `downsampling_factor`-th sample after the first. */
template<class T>
Vector<T> Downsample(
  const Vector<T>& vector,
  const Int downsampling_factor) noexcept
{
  ASSERT(downsampling_factor >= 1);
  Vector<T> output;
  for (size_t i=0; i<vector.length(); i += downsampling_factor)
  {
    output.PushBack(vector[i]);
  }
  return output;
}

/**
 This is equivalent to Matlab's from:to. E.g. 3:5 outputs a vector [3,4,5].
 TODO: Implement fractional input.
 */
template<class T>
Vector<T> ColonOperator(
  const T from,
  const T to) noexcept
{
  if ((to-from) < T(0.0))
  {
    return EmptyVector<T>();
  }
  else
  {
    const size_t vector_length = Floor(to-from+1);
    Vector<T> output(vector_length);
    for (size_t i=0; i<vector_length; ++i)
    {
      output[i] = ((T) i) + ((T) from);
    }
    return output;
  }
}

/**
 This is equivalent to Matlab's from:step:to. E.g. 3:2:6 outputs a vector
 [3,4,6].
 */
// TODO: implement negative step and fractional input.
template<class T>
Vector<T> ColonOperator(
  const T from,
  const T step,
  const T to) noexcept
{
  ASSERT(std::isgreater(step, 0));
  Vector<T> output;
  output.PushBack(from);
  size_t i = 1;
  while (std::islessequal(((T) i)*step+from, to)) {
    output.PushBack(((T) i++)*step+from);
  }
  return output;
}


/**
 Returns elements of vector `vector` from from_id to to_id
 (including extremes).
 */
template<class T>
Vector<T> Elements(
  const Vector<T>& vector,
  const size_t from_id,
  const size_t to_id) noexcept
{
  return Vector<T>(vector.begin() + from_id,
                   vector.begin() + to_id+1);
}



template<class T>
Vector<T> GetSegment(
  const Vector<T>& vector,
  const size_t subset_id,
  const size_t subset_length,
  bool zeropad_if_shorter = false) noexcept
{
  const size_t size = vector.length();

  const size_t from_sample = subset_id * subset_length;
  if (from_sample >= size)
  {
    if (zeropad_if_shorter)
    {
      return Zeros<T>(subset_length);
    }
    else
    {
      return Vector<T>(); // Return empty vector
    }
  }
  const size_t to_sample = Min(from_sample + subset_length - 1, size - 1);
  const size_t actual_length = to_sample - from_sample + 1;
  if (zeropad_if_shorter && actual_length < subset_length)
  {
    Vector<T> output(subset_length);
    ZeroPad(Elements(vector, from_sample, to_sample), output);
    return std::move(output);
  }
  else
  {
    return std::move(Elements(vector, from_sample, to_sample));
  }
}

/**
 Multiplies all the elements in the vector. Equivalent to Matlab's
 prod(vector).
 */
template<typename T, size_t length>
T Prod(
  const Vector<T, length>& vector) noexcept
{
  T output = (T) 1.0;
  for (auto& element : vector)
  {
    output *= element;
  }
  return output;
}

/** Dot product between two vectors. Equivalent to Matlab's dot(a,b) */
template<class T>
T Dot(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  const size_t num_elements = vector_a.length();
  ASSERT(vector_a.length() == vector_b.length());
  T output = (T) 0.0;
  for (size_t i=0; i<num_elements; ++i)
  {
    output += vector_a[i]*vector_b[i];
  }
  return output;
}

template<class T, size_t length>
T Norm(const Vector<T,length>& vector, T l_norm = 2.0) noexcept
{
const size_t num_elements = vector.length();
  T output = 0.0;
  for (size_t i=0; i<num_elements; ++i) {
    output += std::pow(std::fabs(vector[i]), l_norm);
  }
  return std::pow(output, 1.0/l_norm);
}

template<typename T, size_t length>
void Print(const Vector<T,length>& vector) noexcept
{
  std::cout<<"\n------------\n";
  for (auto iter = vector.begin(); iter != vector.end(); ++iter)
  {
    std::cout<<*iter<<std::endl;
  }
  std::cout<<"------------\n";
}

/** Returns a real vector of `length` ones. */
template<typename T>
Vector<T> Ones(size_t length) noexcept {
  return Vector<T>(length, static_cast<T>(1.0));
}



/** Equivalent to Matlab's linspace(min, max, num_elements); */
template<typename T>
Vector<T> LinSpace(
  T min,
  T max,
  size_t num_elements) noexcept
{
  if (num_elements <= 1) { return UnaryVector(max); }
  T interval = max-min;
  T slot = interval / ((T) (num_elements-1));
  Vector<T> output(num_elements);
  for (size_t i=0; i<num_elements; ++i)
  {
    output[i] = min + slot*((T) i);
  }
  return output;
}

template<typename T>
Vector<T> Hann(
  const size_t length) noexcept
{
  Vector<T> w = Ones<T>(length);
  for (size_t i=0; i<length; ++i) {
    w[i] = (1.0-cos(2.0*PI*((T)i)/((T)(length-1))))/2.0;
  }
  return w;
}

/** Returns a Hamming window of length `length' */
template<typename T>
Vector<T> Hamming(
  const size_t length) noexcept
{
  const T alpha = 0.54;
  const T beta = 1.0-alpha;
  Vector<T> w = Ones<T>(length);
  for (size_t i=0; i<length; ++i) {
    w[i] = alpha-beta*cos(2.0*PI*i/(length-1));
  }
  return w;
}

template<typename T>
Vector<T> TukeyWin(
  const size_t length,
  const T ratio) noexcept
{
  if (length == 1) { return UnaryVector<T>(1.0); }
  if (ratio <= 0) { return Ones<T>(length); }
  else if (ratio >= 1.0) { return Hann<T>(length); }
  else {
    Vector<T> t = LinSpace(0.0, 1.0, length);
    // Defines period of the taper as 1/2 period of a sine wave.
    T per = ratio/2.0;
    size_t tl = floor(per*(((T) length)-1.0))+1;
    size_t th = length-tl+1;
    // Window is defined in three sections: taper, constant, taper
    // w1 = ((1+cos(PI/per*(t(1:tl) - per)))/2);
    Vector<T> w = Ones<T>(length);
    for (size_t i=0; i<tl; ++i) {
      w[i] = (1.0+cos(PI/per*(t[i] - per)))/2.0;
    }
    for (size_t i=th-1; i<length; ++i) {
      w[i] = (1.0+cos(PI/per*(t[i] - 1.0 + per)))/2.0;
    }
    return w;
  }
}




template <typename T, size_t length>
T Sum(const Vector<T,length>& input) noexcept
{
  T output((T) 0.0);
  for (auto iter = input.begin(); iter != input.end(); iter++)
  {
    output += *iter;
  }
  return output;
}


template <typename T, size_t length>
/** Equivalent to Matlab's mean(input) */
T Mean(const Vector<T,length>& input) noexcept
{
  return Sum(input) / ((T) input.length());
}




/**
 Returns the geometric mean of the input vector. Equivalent
 to Matlab's geomean(input)
 **/
template<typename T, size_t length>
T Geomean(
  const Vector<T,length>& input) noexcept
{
// TODO: Throw error for negative entries
  return Pow(Prod(input), 1.0/((T)input.length()));
}

/**
 Weighted mean. Not implemented in Matlab (but should be). The weights are
 normalised inside the function. Hence Mean(input, ones(N)) gives the same
 result as Mean(input, ones(N)/N).
 */
template<typename T, size_t length>
T Mean(
  const Vector<T,length>& input,
  const Vector<T,length>& weights) noexcept
{
  ASSERT(input.length() == weights.length());
  ASSERT(IsNonNegative(weights));
  
  // Normalise the weigths
  Vector<T,length> normalised_weights = Multiply(weights, 1.0/Sum(weights));
  ASSERT(Sum(normalised_weights) == 1.0);
  return Sum(Multiply(input, normalised_weights));
}

/**
 Returns the standard deviation of the `input` vector. Equivalent to Matlab's
 std(input). This includes the correction for having an unbiased estimator.
 */
template<typename T, size_t length>
T Std(
  const Vector<T,length>& input) noexcept
{
  return sqrt(Var(input));
}

/** Var (unbiased estimator) */
template<typename T, size_t length>
T Var(
  const Vector<T,length>& input) noexcept
{
  T mean = Mean(input);
  T output(0.0);
  for (size_t i=0; i<input.length(); ++i)
  {
    output += pow(input[i] - mean,2.0);
  }
  return output/((T) (input.length()-1));
}

/** Weighted var (biased estimator) */
template<typename T, size_t length>
T Var(
  const Vector<T,length>& input,
  const Vector<T,length>& weights) noexcept
{
  ASSERT(IsNonNegative(weights));

  T weighted_mean = Mean(input, weights);
  Vector<T> temp = Pow(Add(input, -weighted_mean), 2.0);
  Vector<T> norm_weights = Multiply<T>(weights, 1.0/Sum(weights));

  return (Sum(Multiply(norm_weights, temp)));
}

/** Converts roots to polynomial. Equivalent to Matlab's poly(roots) */
template<typename T>
Vector<Complex<T>> Poly(const Vector<Complex<T>> roots) noexcept
{
  Int n(Length(roots));
  Vector<Complex<T>> output = Zeros<Complex<T>>(n+1);
  // c = [1 zeros(1,n,class(x))];
  output[0] = 1.0;
  // for j=1:n
  for (Int j=1; j<=n; ++j)
  {
    // c(2:(j+1)) = c(2:(j+1)) - e(j).*c(1:j);
    Vector<Complex<T>> temp(output);
    for (Int i=2; i<=(j+1); ++i)
    {
      output[i-1] = temp[i-1] - roots[j-1]*temp[i-2];
    }
  }
  return output;
}

template<typename T>
Vector<Complex<T>> Poly(
  const Vector<T> roots) noexcept
{
  return Poly(ComplexVector(roots));
}


//
///** Returns true if all elements are non negative */
//bool IsNonNegative(const Vector<T>& input) noexcept;
//
//Matrix<T> Cov(const Vector<T>& x,
//                         const Vector<T>& y) noexcept;
//
//Matrix<T> Cov(const Vector<Vector<T> >& input) noexcept;
//
//T CovElement(const Vector<T>& x,
//                        const Vector<T>& y) noexcept;
//
///**
// Returns a vector containing the cumulative sum of
// the elements of X. Equivalent to Matlab's cumsum(input)
// */
//Vector<T> CumSum(const Vector<T>& input) noexcept;
//
/** Splits signal up into (overlapping) frames */
template<typename T>
Vector<Vector<T> > Enframe(
  const Vector<T>& input,
  const Vector<T>& window,
  const Int frame_increment) noexcept
{
  Vector<Vector<T>> output;
  size_t i = 0;
  while ((i + window.length())<=input.length()) {
    size_t from_sample = i;
    size_t to_sample = i + window.length()-1;
    
    ASSERT(from_sample>=0 && from_sample<input.length());
    ASSERT(to_sample>=0 && to_sample<input.length());
    
    output.PushBack(Multiply(Elements(input, from_sample, to_sample), window));
    
    i = i + frame_increment;
  }
  return std::move(output);
}

template<typename T>
Vector<T> OverlapAdd(
  const Vector<Vector<T> >& frames,
  const Vector<T>& window,
  const size_t frame_increment) noexcept
{
  const size_t num_frames = frames.length();
  Vector<T> output(window.length()+(num_frames-1)*frame_increment);
  for (Int frame_i=0; frame_i<num_frames; ++frame_i)
  {
    if (frames[frame_i].length() != window.length())
    {
      ASSERT_WITH_MESSAGE(false, "Frame length different from window length");
    }
    for (Int k=0; k<(Int)window.length(); ++k)
    {
      output[frame_i*frame_increment+k] += window[k] * frames[frame_i][k];
    }
  }
  return output;
}

template<typename T, size_t length>
Vector<Complex<T>,length> ConvertToComplex(
  Vector<T,length> input) noexcept
{
  Vector<Complex<T>,length> output;
  for (Int i=0; i<(Int)input.length(); ++i)
  {
    output.PushBack(Complex(input[i], 0.0));
  }
  return output;
}


template<typename T>
Vector<T> CumSum(
  const Vector<T>& input) noexcept
{
  const size_t N = input.length();
  Vector<T> output(input.length());
  output[N-1] = Sum(input);
  for (size_t i=N-2; i>=0; --i)
  {
    output[i] = output[i+1]-input[i+1];
  }
  return output;
}

/** Splits a string using a delimiter. */
inline Vector<std::string> Split(
  const std::string &s,
  char delim) noexcept
{
  Vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim))
  {
    elems.PushBack(item);
  }
  return elems;
}
  

} /**< namespace mcl */
