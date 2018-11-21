/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vector.h"

namespace mcl
{


/** Fall-back multiply by a constant in case of no available optimisations. */
template<typename T, size_t length>
inline void MultiplySerial(
  const Vector<T,length>& input,
  const T gain,
  Vector<T,length>& output) noexcept
{
  ASSERT(input.length() == output.length());
  for (size_t i=0; i<output.length(); ++i)
  {
    output[i] = input[i]*gain;
  }
}

template<typename T, size_t length>
void MultiplyAddSerial(
  const Vector<T,length>& input_to_multiply,
  const T gain,
  const Vector<T,length>& input_to_add,
  Vector<T,length>& output) noexcept
{
  ASSERT(input_to_multiply.length() == input_to_add.length());
  ASSERT(input_to_add.length() == output.length());
  for (size_t i=0; i<input_to_multiply.length(); ++i)
  {
    output[i] = input_to_multiply[i]*gain + input_to_add[i];
  }
}



template<typename T, size_t length>
struct MathIntrinsics
{
};

template<size_t length>
struct MathIntrinsics<double,length>
{
  static void Multiply(
    const Vector<double>& input,
    const double gain,
    Vector<double>& output) noexcept {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
    vDSP_vmulD(input_data, 1,
               &gain, 0,
               output_data, 1,
               num_samples);
#else
    MultiplySerial(input, gain, output);
#endif
  }
  
  static void MultiplyAdd(
    const Vector<double,length>& input_to_multiply,
    const double gain,
    const Vector<double,length>& input_to_add,
    Vector<double,length>& output) noexcept
  {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
    vDSP_vmaD(
      input_data_mult, 1,
      &gain, 0,
      input_data_add, 1,
      output_data, 1,
      num_samples);
#else
    MultiplyAddSerial(input_to_multiply, gain, input_to_add, output);
#endif
  }
};


template<size_t length>
struct MathIntrinsics<float,length>
{
  static void Multiply(
    const Vector<float>& input,
    const float gain,
    Vector<float>& output) noexcept
  {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
    vDSP_vmul(
      input_data, 1,
      &gain, 0,
      output_data, 1,
      num_samples);
#else
    MultiplySerial(input, gain, output);
#endif
  }
  
  static void MultiplyAdd(
    const Vector<float,length>& input_to_multiply,
    const float gain,
    const Vector<float,length>& input_to_add,
    Vector<float,length>& output) noexcept
  {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
    vDSP_vma(
      input_data_mult, 1,
      &gain, 0,
      input_data_add, 1,
      output_data, 1,
      num_samples);
#else
    MultiplyAddSerial(input_to_multiply, gain, input_to_add, output);
#endif
  }
};


} /**< namespace mcl */
