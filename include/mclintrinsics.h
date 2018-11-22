/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vector.h"
#include "pointwiseop.h"

namespace mcl
{




template<typename T, size_t length>
struct MathIntrinsics
{
  static inline void Multiply(
    const Vector<T,length>& input,
    const T gain,
    Vector<T,length>& output) noexcept
  {
    MultiplySerial(input, gain, output);
  }
  
  static inline void Add(
    const Vector<T,length>& input_a,
    const Vector<T,length>& input_b,
    Vector<T,length>& output) noexcept
  {
    AddSerial(input_a, input_b, output);
  }
  
  static inline void Multiply(
    const Vector<T,length>& input_a,
    const Vector<T,length>& input_b,
    Vector<T,length>& output) noexcept
  {
    MultiplySerial(input_a, input_b, output);
  }
  
  static inline void MultiplyAdd(
    const Vector<T,length>& input_to_multiply,
    const T gain,
    const Vector<T,length>& input_to_add,
    Vector<T,length>& output) noexcept
  {
    MultiplyAddSerial(input_to_multiply, gain, input_to_add, output);
  }
};

template<size_t length>
struct MathIntrinsics<double,length>
{
  
  static inline void Add(
    const Vector<double,length>& input_a,
    const Vector<double,length>& input_b,
    Vector<double,length>& output) noexcept
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
  
  static inline void Multiply(
    const Vector<double,length>& input,
    const double gain,
    Vector<double,length>& output) noexcept
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
  
  static inline void Multiply(
    const Vector<double,length>& input_a,
    const Vector<double,length>& input_b,
    Vector<double,length>& output) noexcept
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
  
  static inline void MultiplyAdd(
    const Vector<double,length>& input_to_multiply,
    const double gain,
    const Vector<double,length>& input_to_add,
    Vector<double,length>& output) noexcept
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
};


template<size_t length>
struct MathIntrinsics<float,length>
{
  static inline void Add(
    const Vector<float,length>& input_a,
    const Vector<float,length>& input_b,
    Vector<float,length>& output) noexcept
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
  
  static inline void Multiply(
    const Vector<float,length>& input,
    const float gain,
    Vector<float,length>& output) noexcept
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
  
  static inline void Multiply(
    const Vector<float,length>& input_a,
    const Vector<float,length>& input_b,
    Vector<float,length>& output) noexcept
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
  
  static inline void MultiplyAdd(
    const Vector<float,length>& input_to_multiply,
    const float gain,
    const Vector<float,length>& input_to_add,
    Vector<float,length>& output) noexcept
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
};



} /**< namespace mcl */
