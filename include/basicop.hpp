/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vectorop.hpp"
#include "pointwiseop.hpp"
#include <limits>

namespace mcl
{
  
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
  
  
template<class T>  
inline T Min(
  const Vector<T>& input)
{
  return input[MinIndex(input)];
}


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
  
  
template<class T>
inline T Max(
  const Vector<T>& input) noexcept
{
  return input[MaxIndex(input)];
}


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

  
  
} /**< namespace mcl */
