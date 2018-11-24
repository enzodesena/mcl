/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

//#include "firfilter.h"


namespace mcl {

  
  


void FilterAppleDsp(
  const Vector<double>& input_data,
  Vector<double>& output_data) noexcept {
  if (num_samples < length_ || (num_samples+length_-1) > MCL_MAX_VLA_LENGTH) {
    FilterSerial(input_data, num_samples, output_data);
    return;
  }
  
  MCL_STACK_ALLOCATE(mcl::T, padded_data, num_samples+length_-1); // TODO: handle stack overflow
  GetExtendedInput(input_data, num_samples, padded_data);
  
  vDSP_convD(padded_data, 1,
             coefficients_.data()+length_-1, -1,
             output_data, 1,
             num_samples, length_);
  
  // Reorganise state for the next run
  for (Int i=0; i<length_; ++i) {
    delay_line_[i] = input_data[num_samples-1-i];
  }
  
  counter_ = length_-1;
}
#endif



void FirFilter::Filter(const T* __restrict input_data,
                       const Int num_samples,
                       T* __restrict output_data) noexcept {
  if (updating_) { UpdateCoefficients(); }
  if (length_ == 1) {
    delay_line_[0] = input_data[num_samples-1];
    mcl::Multiply(input_data, num_samples, coefficients_[0], output_data);
    return;
  }
#if defined(MCL_APPLE_ACCELERATE)
  FilterAppleDsp(input_data, num_samples, output_data);
#elif defined(MCL_AVX_ACCELERATE) || defined(MCL_NEON_ACCELERATE)
  if (num_samples < length_ || (num_samples+length_-1) > MCL_MAX_VLA_LENGTH) {
    FilterSerial(input_data, num_samples, output_data);
    return;
  } else {
#if defined(MCL_ENVWINDOWS) // defined(MCL_AVX_ACCELERATE)
    // Some Intel CPUs do not support AVX instructions.
    // Here we check whether they AVX supported or not, and in that case
    // we filter serially.
    int cpu_info[4];
    __cpuid(cpu_info, 1);
    bool avx_supported = cpu_info[2] & (1 << 28) || false;
    if (avx_supported) {
      FilterSerial(input_data, num_samples, output_data);
      return;
    }
#endif
    
    MCL_STACK_ALLOCATE(float, extended_input_data, num_samples+length_-1); // TODO: handle stack overflow
    MCL_STACK_ALLOCATE(float, output_data_float, num_samples); // TODO: handle stack overflow
    GetExtendedInput<float>(input_data, num_samples, extended_input_data);
    
#ifdef MCL_AVX_ACCELERATE
    const Int batch_size = 8;
    ALIGNED(16) __m256 input_frame;
    ALIGNED(16) __m256 coefficient;
    ALIGNED(16) __m256 accumulator;
    
    for(Int n=0; (n+batch_size)<=num_samples; n+=batch_size) {
      accumulator = _mm256_setzero_ps();
      for(Int k=0; k<length_; k++) {
        coefficient = _mm256_set1_ps((float) coefficients_[length_ - k - 1]);
        input_frame = _mm256_loadu_ps(extended_input_data + n + k);
        accumulator = _mm256_add_ps(_mm256_mul_ps(coefficient, input_frame), accumulator);
        // accumulator = _mm256_fmadd_ps(coefficient, input_frame, accumulator);
      }
      _mm256_storeu_ps(output_data_float+n, accumulator);
    }
#else // MCL_NEON_ACCELERATE
    const Int batch_size = 4;
    float32x4_t input_frame;
    float32x4_t coefficient;
    float32x4_t accumulator;
    
    for(Int n=0; (n+batch_size)<=num_samples; n+=batch_size) {
      accumulator = vdupq_n_f32(0.0);
      for(Int k=0; k<length_; k++) {
        coefficient = vdupq_n_f32((float) coefficients_[length_ - k - 1]);
        input_frame = vld1q_f32(extended_input_data + n + k);
        accumulator = vmlaq_f32(accumulator, coefficient, input_frame);
      }
      vst1q_f32(output_data_float+n, accumulator);
    }
#endif
    
    const Int num_samples_left = num_samples % batch_size;
    const Int num_samples_completed = num_samples - num_samples_left;
    
    for (Int n=0; n<num_samples_completed; ++n) {
      output_data[n] = (T) output_data_float[n];
    }
    
    for (Int n=num_samples_completed; n<num_samples; ++n) {
      output_data[n] = 0.0;
      for (Int p=0; p<length_; ++p) {
        output_data[n] += coefficients_[length_-p-1] * extended_input_data[n+p];
      }
    }
    
    // Reorganise state for the next run
    for (Int i=0; i<length_; ++i) {
      delay_line_[i] = input_data[num_samples-1-i];
    }
    counter_ = length_-1;
  }
#else // defined(MCL_NO_ACCELERATE)
  FilterSerial(input_data, num_samples, output_data);
#endif
}
  

  
  

  
  

  
} // namespace mcl




