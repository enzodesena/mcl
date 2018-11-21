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

namespace mcl {
  
/**
 Returns the index associated to the maximum value in the vector. The index
 counts starting from 0. If there are two maxima,
 the index of the first one is returned.
 */
template<class T>
Int MinIndex(
  const Vector<T>& input) noexcept
{
  T min_value = std::numeric_limits<T>::max();
  Int min_index = 0;
  for (Int i=0; i<(Int)input.length(); ++i)
  {
    if (input[i] < min_value)
    {
      min_value = input[i];
      min_index = i;
    }
  }
  return min_index;
}
  
/** Returns the maximum value of the vector. */
template<class T>  
T Min(
  const Vector<T>& input)
{
  return input[MinIndex(input)];
}


/** 
 Returns the index associated to the maximum value in the vector. The index
 counts starting from 0. If there are two maxima, 
 the index of the first one is returned.
 */
template<class T>
Int MaxIndex(
  const Vector<T>& input) noexcept
{
  return MinIndex(Opposite(input));
}
  
template<>
Int MaxIndex<UInt>(const Vector<UInt>& input) noexcept {
  return MinIndex(Opposite(Convert<UInt,Int>(input)));
}
  
/** Returns the maximum value of the vector. */
template<class T>
T Max(
  const Vector<T>& input) noexcept
{
  return input[MaxIndex(input)];
}


/** 
 Returns the indexes of the local peaks in the vector.
 Equivalent to Matlab's findpeaks.
 */
template<typename T>
Vector<size_t> FindPeaksIndexes(
  const Vector<T>& vector,
  const T min_peak_height = std::numeric_limits<T>::min())
{
  // Allocate new vectors for the indexes of the local maxima
  Vector<size_t> indexes;
  for (size_t i=1; i<vector.length()-1; ++i)
  {
    if ((vector[i] > min_peak_height) &
        (vector[i] > vector[i-1]) &
        (vector[i] > vector[i+1]))
    {
      indexes.PushBack(i);
    }
  }
  return indexes;
}
  

/** 
 Returns the values local peaks in the vector.
 Equivalent to Matlab's findpeaks.
 */
template<typename T>
Vector<T> FindPeaks(
  const Vector<T>& vector,
  const T min_peak_height = std::numeric_limits<T>::min())
{
  Vector<size_t> indexes = FindPeaksIndexes(vector, min_peak_height);
  Vector<T> output(indexes.length());
  for (size_t i=0; i<indexes.length(); ++i)
  {
    output[i] = vector[indexes[i]];
  }
  return output;
}

  
  
bool BasicOpTest();
  
} /**< namespace mcl */
