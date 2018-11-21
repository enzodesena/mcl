/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "vectorop.h"
#include <fstream>
#include "mcltypes.h"
#include "pointwiseop.h"
#include "transformop.h"
#include <cmath>
#include "comparisonop.h"
#include <vector>
#include <cassert>

#ifdef MCL_APPLE_ACCELERATE
  #include <Accelerate/Accelerate.h>
#endif

namespace mcl {
  
void Multiply(const Real* input_data,
                      const Int num_samples,
                      const Real gain,
                      Real* output_data) noexcept {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  #ifdef MCL_DATA_TYPE_DOUBLE
  vDSP_vmulD(input_data, 1,
             &gain, 0,
             output_data, 1,
             num_samples);
  #else
  vDSP_vmul(input_data, 1,
            &gain, 0,
            output_data, 1,
            num_samples);
  #endif
#else
  for (Int i=0; i<num_samples; ++i) { output_data[i] = input_data[i]*gain; }
#endif
}
  
void MultiplyAdd(const Real* input_data_mult, const Real gain,
                 const Real* input_data_add, const Int num_samples,
                 Real* output_data) noexcept {
#if defined(MCL_APPLE_ACCELERATE_MMA) && MCL_APPLE_ACCELERATE_MMA
  #ifdef MCL_DATA_TYPE_DOUBLE
  vDSP_vmaD(input_data_mult, 1,
            &gain, 0,
            input_data_add, 1,
            output_data, 1,
            num_samples);
  #else
  vDSP_vma(input_data_mult, 1,
           &gain, 0,
           input_data_add, 1,
           output_data, 1,
           num_samples);
  #endif
#else
  for (Int i=0; i<num_samples; ++i) {
    output_data[i] = input_data_mult[i]*gain + input_data_add[i];
  }
#endif
}

  
  
Vector<std::string> Split(const std::string &s, char delim) noexcept {
  Vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
  
  
bool IsNonNegative(const Vector<Real>& input) noexcept {
  const Int num_elem = input.length();
  for (Int i=0; i<num_elem; ++i) {
    if (input[i] < 0.0) { return false; }
  }
  return true;
}
  
Matrix<Real> Cov(const Vector<Real>& x,
                 const Vector<Real>& y) noexcept {
  Vector<Vector<Real> > input(2);
  input[0] = x;
  input[1] = y;
  return Cov(input);
}
  
Matrix<Real> Cov(const Vector<Vector<Real> >& input) noexcept {
  const Int N = input.length();
  Matrix<Real> output(N, N);
  for (Int i=0; i<N; ++i) {
    for (Int j=0; j<N; ++j) {
      output.SetElement(i, j, CovElement(input[i], input[j]));
    }
  }
  return output;
}
  
Real CovElement(const Vector<Real>& x,
                const Vector<Real>& y) noexcept {
  const Int N = x.length();
  ASSERT(N == (Int)y.length());
  
  Real output = Sum(Multiply(Add(x, -Mean(x)), Add(y, -Mean(y))));
  // In case N>1 use the unbiased estimator of covariance.
  output = (N > 1) ? output/((Real) (N-1)) : output/((Real) (N));
  return output;
}
  
Vector<Real> CumSum(const Vector<Real>& input) noexcept {
  const Int N = input.length();
  Vector<Real> output(input.length());
  
  output[N-1] = Sum(input);
  for (Int i=N-2; i>=0; --i) { output[i] = output[i+1]-input[i+1]; }
  
  return output;
}
  


  
} // namespace mcl
