/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "basicop.hpp"

#if defined(MCL_APPLE_ACCELERATE)
  #include <Accelerate/Accelerate.h>
#elif defined(MCL_AVX_ACCELERATE)
  #include <pmmintrin.h>
  #include <xmmintrin.h>
  #include <immintrin.h>
#endif

#ifdef MCL_NEON_ACCELERATE
  #include "arm_neon.hpp"
#endif

#ifdef MCL_ENVWINDOWS
  #define ALIGNED(n) __declspec(align(n))
#else
  #define ALIGNED(n) __attribute__ ((aligned (n)))
#endif


namespace mcl
{
  

template<typename T>
inline void Add(
  const Vector<Complex<T>>& input_a,
  const Vector<Complex<T>>& input_b,
  Vector<Complex<T>>& output) noexcept
{
  AddSerial(input_a, input_b, output);
}


inline void Add(
  const Vector<double>& input_a,
  const Vector<double>& input_b,
  Vector<double>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vaddD(
    &input_a[0], 1,
    &input_b[0], 1,
    &output[0], 1,
    output.length());
#else
  AddSerial(input_a, input_b, output);
#endif
}


inline void Add(
  const Vector<float>& input_a,
  const Vector<float>& input_b,
  Vector<float>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vadd(
    &input_a[0], 1,
    &input_b[0], 1,
    &output[0], 1,
    output.length());
#else
  AddSerial(input_a, input_b, output);
#endif
}


/**
 Returns the point by point multiplication of the vector with the gain.
 Equivalent to Matlab's vector_a.*gain.
 */

inline void Multiply(
  const Vector<double>& input,
  const double gain,
  Vector<double>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vmulD(&input[0], 1,
             &gain, 0,
             &output[0], 1,
             output.length());
#else
  MultiplySerial(input, gain, output);
#endif
}

/**
 Returns the point by point multiplication of the vector with the gain.
 Equivalent to Matlab's vector_a.*gain.
 */

inline void Multiply(
  const Vector<float>& input,
  const float gain,
  Vector<float>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vmul(
    &input[0], 1,
    &gain, 0,
    &output[0], 1,
    output.length());
#else
  MultiplySerial(input, gain, output);
#endif
}

template<typename T>
inline void Multiply(
  const Vector<Complex<T>>& input,
  const Complex<T> gain,
  Vector<Complex<T>>& output) noexcept
{
  MultiplySerial(input, gain, output);
}

inline void Multiply(
  const Vector<double>& input_a,
  const Vector<double>& input_b,
  Vector<double>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vmulD(
    &input_a[0], 1,
    &input_b[0], 1,
    &output[0], 1,
    output.length());
#else
  MultiplySerial(input_a, input_b, output);
#endif
}

inline void Multiply(
  const Vector<float>& input_a,
  const Vector<float>& input_b,
  Vector<float>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vmul(
    &input_a[0], 1,
    &input_b[0], 1,
    &output[0], 1,
    output.length());
#else
  MultiplySerial(input_a, input_b, output);
#endif
}

template<typename T>
inline void Multiply(
  const Vector<Complex<T>>& input_a,
  const Vector<Complex<T>>& input_b,
  Vector<Complex<T>>& output) noexcept
{
  MultiplySerial(input_a, input_b, output);
}

inline void MultiplyAdd(
  const Vector<double>& input_to_multiply,
  const double gain,
  const Vector<double>& input_to_add,
  Vector<double>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vmaD(
    &input_to_multiply[0], 1,
    &gain, 0,
    &input_to_add[0], 1,
    &output[0], 1,
    output.length());
#else
  MultiplyAddSerial(input_to_multiply, gain, input_to_add, output);
#endif
}


inline void MultiplyAdd(
  const Vector<float>& input_to_multiply,
  const float gain,
  const Vector<float>& input_to_add,
  Vector<float>& output) noexcept
{
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  vDSP_vma(
    &input_to_multiply[0], 1,
    &gain, 0,
    &input_to_add[0], 1,
    &output[0], 1,
    output.length());
#else
  MultiplyAddSerial(input_to_multiply, gain, input_to_add, output);
#endif
}






/** Calculates the convolution of the `input` vector with impulse response
  `kernel` vector (length P) and outputs it into `output` vector (length N).
  The length of the input vector should at least be N+P-1.
  You shouldn't think of this as Matlab's conv, but rather as the output of
  an FIR filter where input begins with P samples of the previous step
  (all zeroes if this is the first one) and N new samples (if it has more
  thank N new samples, the excess is ignored).
  The reason why this is implemented in this way is that this is how Apple's
  vDSP_conv is implemented.
  @param[in] input input vector
  @param[in] kernel filter kernel
  @param[out] output output vector
*/
template<typename T>
inline void ConvSerial(
  const Vector<T>& input,
  const Vector<T>& kernel,
  Vector<T>& output) noexcept
{
  const size_t N = output.length();
  const size_t P = kernel.length();
  ASSERT(N+P >= 1 && input.length() >= N+P-1);
  for (size_t n=0; n<N; ++n)
  {
    for (size_t p=0; p<P; ++p)
    {
      output[n] += input[n+p] * kernel[P-1-p];
    }
  }
}


#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA

inline void ConvApple(
  const Vector<double>& padded_input,
  const Vector<double>& coefficients,
  Vector<double>& output)
{
  vDSP_convD(
    &padded_input[0], 1,
    &coefficients[0]+coefficients.length()-1, -1,
    &output[0], 1,
    output.length(),
    coefficients.length());
}
#endif


#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
inline void ConvApple(
  const Vector<float>& input,
  const Vector<float>& kernel,
  Vector<float>& output) noexcept
{
  vDSP_conv(
    &input[0], 1,
    &coefficients[0]+coefficients.length()-1, -1,
    &output[0], 1,
    output.length(),
    coefficients.length());
}
#endif


#if defined(MCL_AVX_ACCELERATE) && MCL_AVX_ACCELERATE
inline void ConvAvx(
  const Vector<float>& input,
  const Vector<float>& kernel,
  Vector<float>& output) noexcept
{
#if defined(MCL_ENVWINDOWS) // defined(MCL_AVX_ACCELERATE)
  // Some Intel CPUs do not support AVX instructions.
  // Here we check whether they AVX supported or not, and in that case
  // we filter serially.
  if (! RuntimeArchInfo::GetInstance().IsAvxSupported())
  {
    ConvSerial(input_data, num_samples, output_data);
    return;
  }
#endif
  
  const size_t batch_size = 8;
  ALIGNED(16) __m256 input_frame;
  ALIGNED(16) __m256 coefficient;
  ALIGNED(16) __m256 accumulator;
  
  for(size_t n=0; (n+batch_size)<=num_samples; n+=batch_size)
  {
    accumulator = _mm256_setzero_ps();
    for(size_t k=0; k<length_; k++)
    {
      coefficient = _mm256_set1_ps((float) coefficients_[length_ - k - 1]);
      input_frame = _mm256_loadu_ps(extended_input_data + n + k);
      accumulator = _mm256_add_ps(_mm256_mul_ps(coefficient, input_frame), accumulator);
    }
    _mm256_storeu_ps(output_data_float+n, accumulator);
  }
  
  const size_t num_samples_completed = num_samples - (num_samples % batch_size);
  
  for (size_t n=num_samples_completed; n<num_samples; ++n)
  {
    for (size_t p=0; p<length_; ++p)
    {
      output_data[n] += coefficients_[length_-p-1] * extended_input_data[n+p];
    }
  }
}
#endif


#if defined(MCL_NEON_ACCELERATE) && MCL_NEON_ACCELERATE
inline void ConvNeon(
  const Vector<float>& input,
  const Vector<float>& kernel,
  Vector<float>& output) noexcept
{
  const Int batch_size = 4;
  float32x4_t input_frame;
  float32x4_t coefficient;
  float32x4_t accumulator;

  for(size_t n=0; (n+batch_size)<=num_samples; n+=batch_size)
  {
    accumulator = vdupq_n_f32(0.0);
    for(size_t k=0; k<length_; k++)
    {
      coefficient = vdupq_n_f32(coefficients_[length_ - k - 1]);
      input_frame = vld1q_f32(extended_input_data + n + k);
      accumulator = vmlaq_f32(accumulator, coefficient, input_frame);
    }
    vst1q_f32(output_data_float+n, accumulator);
  }
  
  const size_t num_samples_completed = num_samples - (num_samples % batch_size);
  
  for (size_t n=num_samples_completed; n<num_samples; ++n)
  {
    for (size_t p=0; p<length_; ++p)
    {
      output_data[n] += coefficients_[length_-p-1] * extended_input_data[n+p];
    }
  }
}
#endif


/** Calculates the convolution of the `input` vector with impulse response
  `kernel` vector (length P) and outputs it into `output` vector (length N).
  The length of the input vector should at least be N+P-1.
  You shouldn't think of this as Matlab's conv, but rather as the output of
  an FIR filter where input begins with P samples of the previous step
  (all zeroes if this is the first one) and N new samples (if it has more
  thank N new samples, the excess is ignored).
  The reason why this is implemented in this way is that this is how Apple's
  vDSP_conv is implemented.
  @param[in] input input vector
  @param[in] kernel filter kernel
  @param[out] output output vector
*/
template<typename T>
inline void Conv(
  const Vector<T>& input,
  const Vector<T>& kernel,
  Vector<T>& output) noexcept
{
  ASSERT(input.length() >= kernel.length()+output.length()-1);
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  MathIntrinsics<T>::ConvApple(input, kernel, output);
#elif defined(MCL_AVX_ACCELERATE) && MCL_AVX_ACCELERATE
  MathIntrinsics<T>::ConvAvx(input, kernel, output);
#elif defined(MCL_NEON_ACCELERATE) && MCL_NEON_ACCELERATE
  MathIntrinsics<T>::ConvNeon(input, kernel, output);
#else
  ConvSerial(input, kernel, output);
#endif
}




template<typename T>
inline Vector<T> Add(
  const Vector<T>& input_a,
  const Vector<T>& input_b) noexcept
{
  Vector<T> output(input_a.length());
  Add(input_a, input_b, output);
  return std::move(output);
}


/**
 Returns the point by point multiplication of the vector with the gain.
 Equivalent to Matlab's vector_a.*gain.
 */
template<typename T>
inline Vector<T> Multiply(
  const Vector<T>& input,
  const T gain) noexcept
{
  Vector<T> output(input.length());
  Multiply(input, gain, output);
  return std::move(output);
}


template<typename T>
inline Vector<T> Multiply(
  const Vector<T>& input_a,
  const Vector<T>& input_b) noexcept
{
  Vector<T> output(input_a.length());
  Multiply(input_a, input_b, output);
  return std::move(output);
}
  
  

 
/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<typename T>
inline void AddScalar(
  const Vector<T>& input,
  const T scalar,
  Vector<T>& output) noexcept
{
  auto input_iter = input.begin();
  auto output_iter = output.begin();
  while (output_iter != output.end())
  {
    *(output_iter++) = *(input_iter++) + scalar;
  }
}

/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<typename T>
inline Vector<T> AddScalar(
  const Vector<T>& vector,
  const T scalar) noexcept
{
  Vector<T> output(vector.length());
  AddScalar(vector, scalar, output);
  return std::move(output);
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


/**
 Returns the point by point subtraction of the two vectors.
 Equivalent to Matlab's vector_a-vector_b.
 */
template<class T>
inline Vector<T> Subtract(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  return Add(vector_a, Opposite(vector_b));
}

} /**< namespace mcl */
