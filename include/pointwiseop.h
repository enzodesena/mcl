/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#ifndef MCL_POINTWISE_H
#define MCL_POINTWISE_H

#include <cassert>
#include "mcltypes.h"
#include <vector>
#include <limits>

namespace mcl {

/**
 Returns the point by point addition of the two vectors.
 Equivalent to Matlab's vector_a+vector_b.
 */
template<class T> 
Vector<T> Add(const Vector<T>& vector_a,
                   const Vector<T>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  
  Vector<T> output((Int)vector_a.length());
  for (Int i=0; i<(Int)vector_a.length(); ++i) {
    output[i] = vector_a[i]+vector_b[i];
  }
  return output;
}
  
void Add(const Real* input_data_a,
                 const Real* input_data_b,
                 Int num_samples,
                 Real* output_data) noexcept;


template<>
inline Vector<Real> Add(const Vector<Real>& vector_a,
                             const Vector<Real>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  
  Vector<Real> output((Int)vector_a.length());
  Add(vector_a.data(), vector_b.data(), (Int)vector_a.length(),
      output.data());
  return output;
}
  
/** Returns the opposite vector.Equivalent to Matlab's -vector. */
template<class T> 
Vector<T> Opposite(const Vector<T>& vector) noexcept {
  // Checking we are not dealing with unsigned types.
  // The assert below responds false to complex. TODO: fix this
  //ASSERT(std::numeric_limits<T>::is_signed);
  
  Vector<T> output(vector.length());
  for (Int i=0; i<(Int)vector.length(); ++i) { output[i] = -vector[i]; }
  return output;
}
  
  
/** Returns the inverse vector.Equivalent to Matlab's 1./vector. */
Vector<Real> Inverse(const Vector<Real>& vector) noexcept;

/** 
 Returns the point by point subtraction of the two vectors.
 Equivalent to Matlab's vector_a-vector_b.
 */
template<class T> 
Vector<T> Subtract(const Vector<T>& vector_a,
                        const Vector<T>& vector_b) noexcept {
  return Add(vector_a, Opposite(vector_b));
}


/** 
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T> 
Vector<T> Multiply(const Vector<T>& vector_a,
                        const Vector<T>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  
  Vector<T> output((Int)vector_a.length());
  for (Int i=0; i<(Int)vector_a.length(); ++i) {
    output[i] = vector_a[i]*vector_b[i];
  }
  return output;
}

void Multiply(const Real* input_data_a, const Real* input_data_b,
              Int num_samples, Real* output_data) noexcept;
  
  
template<>
inline Vector<Real> Multiply(const Vector<Real>& vector_a,
                                  const Vector<Real>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  
  Vector<Real> output((Int)vector_a.length());
  Multiply(vector_a.data(), vector_b.data(), (Int)vector_a.length(),
           output.data());
  return output;
}
  
/** 
 Returns the point by point multiplication of the two vectors.
 Equivalent to Matlab's vector_a.*vector_b.
 */
template<class T>
Vector<T> Divide(const Vector<T>& vector_a,
                      const Vector<T>& vector_b) noexcept {
  ASSERT(vector_a.length() == vector_b.length());
  Vector<T> output((Int)vector_a.length());
  for (Int i=0; i<(Int)vector_a.length(); ++i) {
    output[i] = vector_a[i]/vector_b[i];
  }
  return output;
}


/** Equivalent to Matlab's exp(vector). */
template<class T>
Vector<T> Exp(const Vector<T>& vector) noexcept {
  Int n(vector.length());
  Vector<T> output(vector.length());
  for (Int i=0; i<n; ++i) { output[i] = exp(vector[i]); }
  return output;
}
  
  
  
/** 
 Returns the vector with conjugate entries.
 Equivalent to Matlab's conj(vector).
 */
std::vector<Complex> Conj(const std::vector<Complex>& vector) noexcept;

/** Transform real vector into complex vector with null imaginary part */
std::vector<Complex>
ComplexVector(const Vector<Real>& input) noexcept;

/** Equivalent to Matlab's real(input). */
Vector<Real> RealPart(const std::vector<Complex>& input) noexcept;

/** Equivalent to Matlab's imag(input). */
Vector<Real> Imag(const std::vector<Complex>& input) noexcept;


/** 
 Returns the point-wise poser to exponent.
 Equivalent to Matlab's vector.^exponent
 */
Vector<Real> Pow(const Vector<Real>& vector,
                              Real exponent) noexcept;


/** Equivalent to Matlab's abs(vector) */
Vector<Real> Abs(const Vector<Real>& input) noexcept;

/** Equivalent to Matlab's abs(vector) */
Vector<Real> Abs(const std::vector<Complex>& input) noexcept;

/** Equivalent to Matlab's vector.*(vector>0) */
Vector<Real> HalfWave(const Vector<Real>& vector) noexcept;

/** Equivalent to Matlab's cos(vector) */
Vector<Real> Cos(const Vector<Real>& vector) noexcept;

/** Equivalent to Matlab's sin(vector) */
Vector<Real> Sin(const Vector<Real>& vector) noexcept;
  
/** 
 Returns the natural logarithm of the elements of vector.
 Equivalent to Matlab's log(vector).
 */
Vector<Real> Log(const Vector<Real>& vector) noexcept;
  
/**
 Returns the 10-base logarithm of the elements of vector.
 Equivalent to Matlab's log10(vector).
 */
Vector<Real> Log10(const Vector<Real>& vector) noexcept;
  
std::vector<Int> ConvertToInt(const std::vector<UInt>& vector) noexcept;
  
bool PointWiseOpTest();
  
} /**< namespace mcl */

#endif
