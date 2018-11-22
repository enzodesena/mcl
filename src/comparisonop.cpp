/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "comparisonop.h"
#include "mcltypes.h"
#include "vectorop.h"
#include "pointwiseop.h"
#include "quaternion.h"
#include <vector>




namespace mcl {
  

bool IsEqual(Real num_a, Real num_b, Real precision) {
  if (isnan(num_a) || isnan(num_b)) { return false; }
  return ((Real) std::fabs(((double) num_a) - ((double) num_b))) < precision;
}
  


Vector<bool> IsNan(Vector<Real> input) {
  Vector<bool> output;
  for (Int i=0; i<(Int)input.length(); ++i) {
    output.push_back(IsNan(input[i]));
  }
  return output;
}


Vector<bool> IsInf(Vector<Real> input) {
  Vector<bool> output;
  for (Int i=0; i<(Int)input.length(); ++i) {
    output.push_back(IsInf(input[i]));
  }
  return output;
}

Vector<bool> Not(Vector<bool> input) {
  Vector<bool> output;
  for (Int i=0; i<(Int)input.length(); ++i) {
    output.push_back(!input[i]);
  }
  return output;
}

bool All(Vector<bool> input) {
  for (Int i=0; i<(Int)input.length(); ++i) {
    if (input[i] == false) {
      return false;
    }
  }
  return true;
}

bool None(Vector<bool> input) {
  return All(Not(input));
}

bool Any(Vector<bool> input) {
  for (Int i=0; i<(Int)input.length(); ++i) {
    if (input[i] == true) {
      return true;
    }
  }
  return false;
}

bool IsSmallerOrEqual(const Real num_a, const Real num_b, const Real precision) {
  return num_a <= (num_b + precision);
}

bool IsLargerOrEqual(const Real num_a, const Real num_b, const Real precision) {
  return num_a >= (num_b - precision);
}

bool AreAllSmallerOrEqual(const Vector<Real>& vector_a,
                          const Vector<Real>& vector_b) {
  if ((Int)vector_a.length() != (Int)vector_b.length())
    return false;
  
  for (Int i=0; i<(Int)(Int)vector_a.length(); ++i) {
    if (! IsSmallerOrEqual(vector_a[i], vector_b[i])) { return false; }
  }
  
  return true;
}

bool IsEqual(Complex num_a, Complex num_b, Real precision) {
  return (std::fabs(num_a.real() - num_b.real()) < precision) &
  (std::fabs(num_a.imag() - num_b.imag()) < precision);
}



bool IsEqual(const Vector<Int>& vector_a,
             const Vector<Int>& vector_b) {
  if ((Int)vector_a.length() != (Int)vector_b.length())
    return false;
  
  for (Int i=0; i<(Int)vector_a.length(); ++i) {
    if (vector_a[i] != vector_b[i]) { return false; }
  }
  return true;
}


  
bool IsEqual(const Real* input_data_a,
             const Real* input_data_b,
             const Int num_samples,
             Real precision) {
  for (Int i=0; i<num_samples; ++i) {
    if (! mcl::IsEqual(input_data_a[i], input_data_b[i], precision)) {
      return false;
    }
  }
  return true;
}
  

bool IsEqual(const Real* input_data_a,
             const Vector<Real> input_data_b,
             Real precision) {
  return IsEqual(input_data_a, input_data_b.data(), input_data_b.length(),
                 precision);
}

  
bool IsEqual(const Vector<Real> input_data_b,
             const Real* input_data_a,
             Real precision) {
  return IsEqual(input_data_a, input_data_b.data(), input_data_b.length(),
                 precision);
}
  
} // namespace mcl
