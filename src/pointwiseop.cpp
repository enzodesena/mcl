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

  
void Multiply(const Real* input_data_a,
              const Real* input_data_b,
              Int num_samples,
              Real* output_data) noexcept {
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
  
void Add(const Real* input_data_a,
                 const Real* input_data_b,
                 Int num_samples,
                 Real* output_data) noexcept {
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
  
std::vector<Complex> ComplexVector(const Vector<Real>& input) noexcept {
  std::vector<Complex> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { output[i] = Complex(input[i], 0.0); }
  return output;
}

Vector<Real> RealPart(const std::vector<Complex>& input) noexcept {
  Vector<Real> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { output[i] = input[i].real(); }
  return output;
}

Vector<Real> Imag(const std::vector<Complex>& input) noexcept {
  Vector<Real> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { output[i] = input[i].imag(); }
  return output;
}


Vector<Real> Inverse(const Vector<Real>& vector) noexcept {
  Vector<Real> output(vector.length());
  for (Int i=0; i<(Int)vector.length(); ++i) { output[i] = 1.0/vector[i]; }
  return output;
}


std::vector<Complex> Conj(const std::vector<Complex>& vector) noexcept {
  std::vector<Complex> output(vector.length());
  for (Int i=0; i<(Int)vector.length(); ++i) { 
    output[i] = Conj(vector[i]); 
  }
  return output;
}



Vector<Real> Pow(const Vector<Real>& vector, Real exponent) noexcept {
  Vector<Real> output(vector.length());
  for (Int i=0; i<(Int)vector.length(); ++i) { 
    output[i] = Pow(vector[i], exponent);
  }
  return output;
}


Vector<Real> HalfWave(const Vector<Real>& input) noexcept {
  Vector<Real> output = Zeros<Real>(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) {
    if (input[i] > 0.0)
      output[i] = input[i];
  }
  return output;
}


Vector<Real> Abs(const Vector<Real>& input) noexcept {
  Vector<Real> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { 
    output[i] = std::fabs(input[i]);
  }
  return output;
}


Vector<Real> Abs(const std::vector<Complex>& input) noexcept {
  Vector<Real> output(input.length());
  for (Int i=0; i<(Int)input.length(); ++i) { 
    output[i] = std::abs(input[i]);
  }
  return output;
}


Vector<Real> Log(const Vector<Real>& vector) noexcept {
  Int n(vector.length());
  Vector<Real> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = log(vector[i]); }
  return output;
}
  
Vector<Real> Log10(const Vector<Real>& vector) noexcept {
  Int n(vector.length());
  Vector<Real> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = log10(vector[i]); }
  return output;
}

Vector<Real> Cos(const Vector<Real>& vector) noexcept {
  Int n(vector.length());
  Vector<Real> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = cos(vector[i]); }
  return output;
}
  
Vector<Real> Sin(const Vector<Real>& vector) noexcept {
  Int n(vector.length());
  Vector<Real> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = sin(vector[i]); }
  return output;
}

std::vector<Int> ConvertToInt(const std::vector<UInt>& vector) noexcept {
  const Int length = vector.length();
  std::vector<Int> output(length);
  for (Int i=0; i<length; ++i) { output[i] = (Int) vector[i]; }
  return output;
}
  
} // namespace mcl

