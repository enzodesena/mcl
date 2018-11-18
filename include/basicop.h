/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include "pointwiseop.h"
#include <vector>
#include <limits>

using std::vector;

namespace mcl {
  
/**
 Returns the index associated to the maximum value in the vector. The index
 counts starting from 0. If there are two maxima,
 the index of the first one is returned.
 */
template<class T>
Int MinIndex(const Vector<T>& input) noexcept {
  T min_value = std::numeric_limits<T>::max();
  Int min_index = 0;
  for (Int i=0; i<(Int)input.length(); ++i) {
    if (input[i] < min_value) {
      min_value = input[i];
      min_index = i;
    }
  }
  return min_index;
}
  
/** Returns the maximum value of the vector. */
template<class T>  
T Min(const Vector<T>& input) {
  return input[MinIndex(input)];
}


/** 
 Returns the index associated to the maximum value in the vector. The index
 counts starting from 0. If there are two maxima, 
 the index of the first one is returned.
 */
template<class T>
Int MaxIndex(const Vector<T>& input) noexcept {
  return MinIndex(Opposite(input));
}
  
template<>
Int MaxIndex<UInt>(const std::vector<UInt>& input) noexcept;

  
/** Returns the maximum value of the vector. */
template<class T>
T Max(const Vector<T>& input) noexcept {
  return input[MaxIndex(input)];
}




/** 
 Returns the indexes of the local peaks in the vector.
 Equivalent to Matlab's findpeaks.
 */
std::vector<UInt>
FindPeaksIndexes(const Vector<Real>& vector,
                         const Real min_peak_height = std::numeric_limits<Real>::min());

/** 
 Returns the values local peaks in the vector.
 Equivalent to Matlab's findpeaks.
 */
Vector<Real>
FindPeaks(const Vector<Real>& vector,
                  const Real min_peak_height = std::numeric_limits<Real>::min());

  
bool BasicOpTest();
  
} /**< namespace mcl */
