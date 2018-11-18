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

Vector<Real> Ones(Int length) noexcept {
  return Vector<Real>(length, (Real) 1.0);
}
  
Real Mean(const Vector<Real>& input,
          const Vector<Real>& weights) noexcept {
  ASSERT(input.length() == weights.length());
  ASSERT(IsNonNegative(weights));
  
  // Normalise the weigths
  Vector<Real> normalised_weights = Multiply<Real>(weights, 1.0/Sum(weights));
  ASSERT(IsEqual(Sum(normalised_weights), 1.0));
  return Sum(Multiply(input, normalised_weights));
}

Real Geomean(const Vector<Real>& input) noexcept {
  // TODO: Throw error for negative entries
  return Pow(Prod(input), 1.0/((Real)input.length()));
}

Real Std(const Vector<Real>& input) noexcept {
  return sqrt(Var(input));
}
  
Real Var(const Vector<Real>& input) noexcept {
  Real mean = Mean(input);
  Real output(0.0);
  for (Int i=0; i<(Int)input.length(); ++i) { output += pow(input[i] - mean,2.0); }
  return output/((Real) (input.length()-1));
}
  
Real Var(const Vector<Real>& input,
         const Vector<Real>& weights) noexcept {
  ASSERT(IsNonNegative(weights));
  
  Real weighted_mean = Mean(input, weights);
  Vector<Real> temp = Pow(Add(input, -weighted_mean), 2.0);
  Vector<Real> norm_weights = Multiply<Real>(weights, 1.0/Sum(weights));
  
  return (Sum(Multiply(norm_weights, temp)));
}
  
std::vector<std::string> Split(const std::string &s, char delim) noexcept {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
  
std::vector<Complex> Poly(const std::vector<Complex> roots) noexcept {
  Int n(Length(roots));
  std::vector<Complex> output = Zeros<Complex>(n+1);
  // c = [1 zeros(1,n,class(x))];
  output[0] = 1.0;
  
  // for j=1:n
  for (Int j=1; j<=n; ++j) {
    // c(2:(j+1)) = c(2:(j+1)) - e(j).*c(1:j);
    std::vector<Complex> temp(output);
    for (Int i=2; i<=(j+1); ++i) {
      output[i-1] = temp[i-1] - roots[j-1]*temp[i-2];
    }
  }
  
  return output;
}
  
std::vector<Complex> Poly(const Vector<Real> roots) noexcept {
  return Poly(ComplexVector(roots));
}
  
Vector<Real>
ColonOperator(const Real from, const Real step, const Real to) noexcept {
  ASSERT(std::isgreater(step, 0));
  Vector<Real> output;
  output.push_back(from);
  Int i = 1;
  while (std::islessequal(((Real) i)*step+from, to)) {
    output.push_back(((Real) i++)*step+from);
  }
  return output;
}

Vector<Real> TukeyWin(const Int length, const Real ratio) noexcept {
  if (length == 1) { return UnaryVector<Real>(1.0); }
  if (ratio <= 0) { return Ones(length); }
  else if (ratio >= 1.0) { return Hann(length); }
  else {
    Vector<Real> t = LinSpace(0.0, 1.0, length);
    // Defines period of the taper as 1/2 period of a sine wave.
    Real per = ratio/2.0;
    Int tl = (Int) floor(per*(((Real) length)-1.0))+1;
    Int th = length-tl+1;
    // Window is defined in three sections: taper, constant, taper
    // w1 = ((1+cos(PI/per*(t(1:tl) - per)))/2);
    Vector<Real> w = Ones(length);
    for (Int i=0; i<tl; ++i) {
      w[i] = (1.0+cos(PI/per*(t[i] - per)))/2.0;
    }
    for (Int i=th-1; i<length; ++i) {
      w[i] = (1.0+cos(PI/per*(t[i] - 1.0 + per)))/2.0;
    }
    return w;
  }
}
  
Vector<Real> Hann(const Int length) noexcept {
  Vector<Real> w = Ones(length);
  for (Int i=0; i<length; ++i) {
    w[i] = (1.0-cos(2.0*PI*((Real)i)/((Real)(length-1))))/2.0;
  }
  return w;
}
  
Real Norm(const Vector<Real>& vector, Real l_norm) noexcept {
  const Int num_elements = vector.length();
  Real output = 0.0;
  for (Int i=0; i<num_elements; ++i) {
    output += std::pow(std::fabs(vector[i]), l_norm);
  }
  return std::pow(output, 1.0/l_norm);
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
  std::vector<Vector<Real> > input(2);
  input[0] = x;
  input[1] = y;
  return Cov(input);
}
  
Matrix<Real> Cov(const std::vector<Vector<Real> >& input) noexcept {
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
  
Vector<Real> Hamming(const Int length) noexcept {
  const Real alpha = 0.54;
  const Real beta = 1.0-alpha;
  Vector<Real> w = Ones(length);
  for (Int i=0; i<length; ++i) {
    w[i] = alpha-beta*cos(2.0*PI*i/(length-1));
  }
  return w;
}

std::vector<Vector<Real> > Enframe(const Vector<Real>& input,
                                        const Vector<Real>& window,
                                        const Int frame_increment) noexcept {
  std::vector<Vector<Real> > output;
  Int i = 0;
  while ((i + window.length())<=input.length()) {
    Int from_sample = i;
    Int to_sample = i + window.length()-1;
    
    ASSERT(from_sample>=0 && from_sample<(Int)input.length());
    ASSERT(to_sample>=0 && to_sample<(Int)input.length());
    
    output.push_back(Multiply(Elements(input, from_sample, to_sample),
                              window));
    
    i = i + frame_increment;
  }
  return output;
}

Vector<Real> OverlapAdd(const std::vector<Vector<Real> >& frames,
                             const Vector<Real>& window,
                             const Int frame_increment) noexcept {
  const Int num_frames = frames.length();
  Vector<Real> output(window.length()+(num_frames-1)*frame_increment);
  for (Int frame_i=0; frame_i<num_frames; ++frame_i) {
    if (frames[frame_i].length() != window.length()) {
      ASSERT_WITH_MESSAGE(false, "Frame length different from window length");
    }
    for (Int k=0; k<(Int)window.length(); ++k) {
      output[frame_i*frame_increment+k] += window[k] * frames[frame_i][k];
    }
  }
  return output;
}
  
std::vector<Complex> ConvertToComplex(Vector<Real> input) noexcept {
  std::vector<Complex> output;
  for (Int i=0; i<(Int)input.length(); ++i) {
    output.push_back(Complex(input[i], 0.0));
  }
  return output;
}
  
} // namespace mcl
