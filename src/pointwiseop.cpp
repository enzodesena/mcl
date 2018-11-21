/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "pointwiseop.h"
#include "mcltypes.h"
#include "vectorop.h"
#include <vector>

#ifdef MCL_APPLE_ACCELERATE
  #include <Accelerate/Accelerate.h>
#endif

namespace mcl {

  
void Multiply(const T* input_data_a,
              const T* input_data_b,
              Int num_samples,
              T* output_data) noexcept {
#ifdef MCL_APPLE_ACCELERATE
  #if MCL_DATA_TYPE_DOUBLE
  vDSP_vmulD(input_data_a, 1,
             input_data_b, 1,
             output_data, 1,
             num_samples);
  #else
  vDSP_vmul(input_data_a, 1,
            input_data_b, 1,
            output_data, 1,
            num_samples);
  #endif
#else
  for (Int i=0; i<num_samples; ++i) {
    output_data[i] = input_data_a[i]*input_data_b[i];
  }
#endif
}
  
void Add(const T* input_data_a,
                 const T* input_data_b,
                 Int num_samples,
                 T* output_data) noexcept {
#ifdef MCL_APPLE_ACCELERATE
  #if MCL_DATA_TYPE_DOUBLE
    vDSP_vaddD(input_data_a, 1,
               input_data_b, 1,
               output_data, 1,
               num_samples);
  #else
    vDSP_vadd(input_data_a, 1,
              input_data_b, 1,
              output_data, 1,
              num_samples);
  #endif
#else
  for (Int i=0; i<num_samples; ++i) {
    output_data[i] = input_data_a[i]+input_data_b[i];
  }
#endif
}
  

Vector<T> Inverse(const Vector<T>& vector) noexcept {
  Vector<T> output(vector.length());
  for (Int i=0; i<(Int)vector.length(); ++i) { output[i] = 1.0/vector[i]; }
  return output;
}



Vector<T> HalfWave(const Vector<T>& input) noexcept {
  Vector<T> output = Zeros<T>(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) {
    if (input[i] > 0.0)
      output[i] = input[i];
  }
  return output;
}



Vector<T> Abs(const Vector<Complex>& input) noexcept {
  Vector<T> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { 
    output[i] = std::abs(input[i]);
  }
  return output;
}



Vector<T> Cos(const Vector<T>& vector) noexcept {
  Int n(vector.length());
  Vector<T> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = cos(vector[i]); }
  return output;
}
  
Vector<T> Sin(const Vector<T>& vector) noexcept {
  Int n(vector.length());
  Vector<T> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = sin(vector[i]); }
  return output;
}

  
} // namespace mcl

