/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vectorop.h"
#include "mcltypes.h"
#include "pointwiseop.h"
#include <limits>

namespace mcl
{

// Forward declarations
template<typename TOrigin, typename TDestination>
Vector<TDestination> Cast(
  const Vector<TOrigin>& vector) noexcept;
  
template<typename T>
inline Vector<T> Opposite(
  const Vector<T>& input) noexcept;
// End of forward declarations
  
/**
 Returns the index associated to the maximum value in the vector. The index
 counts starting from 0. If there are two maxima,
 the index of the first one is returned.
 */
template<class T>
inline Int MinIndex(
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
inline T Min(
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
inline size_t MaxIndex(
  const Vector<T>& input) noexcept
{
  return MinIndex(Opposite(input));
}

template<>
inline size_t MaxIndex<UInt>(
  const Vector<UInt>& input) noexcept
{
  return MinIndex(Opposite(Cast<UInt,Int>(input)));
}
  
template<>
inline size_t MaxIndex<size_t>(
  const Vector<size_t>& input) noexcept
{
  return MinIndex(Opposite(Cast<size_t,Int>(input)));
}
  
  
/** Returns the maximum value of the vector. */
template<class T>
inline T Max(
  const Vector<T>& input) noexcept
{
  return input[MaxIndex(input)];
}


/** 
 Returns the indexes of the local peaks in the vector.
 Equivalent to Matlab's findpeaks.
 */
template<typename T>
inline Vector<size_t> FindPeaksIndexes(
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
inline Vector<T> FindPeaks(
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
