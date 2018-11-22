/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#pragma once

#ifndef VERY_SMALL
  #define VERY_SMALL (0.0001)
#endif

#if defined(__APPLE__)
#if (__GNUC__ >= 4)
#include <cmath>
#define isnan(x) std::isnan(x)
#else
#include <math.h>
#define isnan(x) __isnand((double)x)
#endif
#endif

#include "mcltypes.h"
#include "quaternion.h"
#include "vector.h"
#include <vector>

namespace mcl {


inline bool IsEqual(
  const double num_a,
  const double num_b)
{
  return num_a == num_b;
}

inline bool IsEqual(
  const float num_a,
  const float num_b)
{
  return num_a == num_b;
}

template<typename T, size_t length_a, size_t length_b>
inline bool IsEqual(
  const Vector<T,length_a>& vector_a,
  const Vector<T,length_b>& vector_b)
{
  if (vector_a.length() != vector_b.length())
  {
    return false;
  }
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end())
  {
    if (*(iter_a++) != *(iter_b++)) return false;
  }
  return true;
}

template<
  typename VectorT,
  size_t length_vector_a,
  size_t length_vector_b,
  typename PrecisionT>
inline bool ConditonCheckerWithPrecision(
  const Vector<VectorT,length_vector_a>& vector_a,
  const Vector<VectorT,length_vector_b>& vector_b,
  bool (*condition_checker)(VectorT, VectorT, PrecisionT),
  PrecisionT precision)
{
  if (vector_a.length() != vector_b.length()) return false;
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end())
  {
    if (! condition_checker(*(iter_a++), *(iter_b++), precision))
    {
      return false;
    }
  }
  return true;
}


  
inline bool IsApproximatelyEqual(
  const double num_a,
  const double num_b,
  const double precision = VERY_SMALL)
{
  return std::fabs(num_a - num_b) < precision;
}

inline bool IsApproximatelyEqual(
  const float num_a,
  const float num_b,
  const float precision = VERY_SMALL);

inline bool IsApproximatelyEqual(
  const float num_a,
  const float num_b,
  const float precision)
{
  return std::abs(num_a - num_b) < precision;
}

template<typename T>
inline bool IsApproximatelyEqual(
  const Complex<T> num_a,
  const Complex<T> num_b,
  const T precision = VERY_SMALL)
{
  return std::abs(num_a - num_b) < precision;
}


template<
  typename VectorT,
  size_t length_vector_a,
  size_t length_vector_b,
  typename PrecisionT>
inline bool IsApproximatelyEqual(
  const Vector<VectorT,length_vector_a>& vector_a,
  const Vector<VectorT,length_vector_b>& vector_b,
  const PrecisionT precision = VERY_SMALL)
{
  return ConditonCheckerWithPrecision(
    vector_a,
    vector_b,
    &IsApproximatelyEqual,
    precision);
}




bool IsSmallerOrEqual(const Real num_a, const Real num_b,
                      const Real precision = VERY_SMALL);

bool IsLargerOrEqual(const Real num_a, const Real num_b,
                     const Real precision = VERY_SMALL);

//bool AreAllSmallerOrEqual(const Vector<Real>& vector_a,
//                          const Vector<Real>& vector_b);

//template<class T>
//bool IsEqual(const Vector<T>& vector_a, const Vector<T>& vector_b,
//             Real precision = VERY_SMALL) noexcept {
//  if ((Int)vector_a.length() != (Int)vector_b.length())
//    return false;
//
//  for (Int i=0; i<(Int)(Int)vector_a.length(); ++i) {
//    if (! IsEqual(vector_a[i], vector_b[i], precision))
//      return false;
//  }
//  return true;
//}


template<typename T, size_t length_a, size_t length_b>
inline bool ConditonChecker(
  const Vector<T,length_a>& vector_a,
  const Vector<T,length_b>& vector_b,
  bool (*condition_checker)(T, T))
{
  if (vector_a.length() != vector_b.length()) return false;
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end())
  {
    if (! condition_checker(*(iter_a++), *(iter_b++)))
    {
      return false;
    }
  }
  return true;
}

//template<typename T, size_t length_a, size_t length_b>
//bool IsEqual(
//  const Vector<T,length_a>& vector_a,
//  const Vector<T,length_b>& vector_b)
//{
//  return ConditonChecker(vector_a, vector_b, [](auto a, auto b){ return a == b; });
//}

template<typename T>
inline bool IsEqual(
  const Quaternion<T>& q_a,
  const Quaternion<T>& q_b)
{
  return q_a.w() == q_b.w() && q_a.x() == q_b.x() &
    q_a.y() == q_b.y() && q_a.z() == q_b.z();
}


template<typename T>
inline bool IsEqual(
  const Point<T>& point_a,
  const Point<T>& point_b,
  const T precision = VERY_SMALL)
{
  return mcl::IsEqual(point_a.x(), point_b.x(), precision) &&
  mcl::IsEqual(point_a.y(), point_b.y(), precision) &&
  mcl::IsEqual(point_a.z(), point_b.z(), precision);
}

template<typename T>
inline bool IsEqual(
  const Vector<Point<T>>& points_a,
  const Vector<Point<T>>& points_b)
{
  const Int num_points = (Int)points_a.length();
  if (num_points != (Int)points_b.length()) { return false; }
  for (Int i=0; i<num_points; ++i) {
    if (! IsEqual(points_a[i], points_b[i])) { return false; }
  }
  return true;
}
  
/** Returns true if num is nan */
template<typename T>
inline bool IsNan(T num)
{
  return isnan(num);
}

/** Returns true if num is nan */
template<typename T>
inline Vector<bool> IsNan(Vector<T> input)
 {
  Vector<bool> output(input.length());
  for (Int i=0; i<input.length(); ++i)
  {
    output[i] = IsNan(input[i]);
  }
  return output;
}

///** Returns opposite bool as input */
//Vector<bool> Not(Vector<bool> input);
//
///** Returns true if all bools are true */
//bool All(Vector<bool> input);
//
///** Returns true if any one of the bools is true */
//bool Any(Vector<bool> input);
//
///** Opposite of All: returns true if none of the inputs are true */
//bool None(Vector<bool> input);

/** Returns true if num is +inf or -inf */
template<typename T>
inline bool IsInf(T num)
{
  return std::isinf(num);
}

///** Returns true if num is +inf or -inf */
//Vector<bool> IsInf(Vector<Real> input);

bool ComparisonOpTest();
  
  
  
} // namespace mcl

