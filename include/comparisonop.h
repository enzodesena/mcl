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


template<typename T>
bool IsApproximatelyEqual(
  T num_a,
  T num_b,
  T precision = VERY_SMALL)
{
  if (isnan(num_a) || isnan(num_b)) return false;
  return std::fabs(num_a - num_b) < precision;
}


//bool IsApproximatelyEqual(
//  double num_a,
//  double num_b,
//  double precision)
//{
//  if (isnan(num_a) || isnan(num_b)) return false;
//  return (std::fabs(num_a - num_b)) < precision;
//}

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


template<typename T, int length_a, int length_b>
bool IsEqual(
  const Vector<T,length_a>& vector_a,
  const Vector<T,length_b>& vector_b)
{
  if (vector_a.length() != vector_b.length()) return false;
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end())
  {
    if (*(iter_a++) != *(iter_b++)) return false;
  }
  return true;
}

template<typename T, int length_a, int length_b>
bool ConditonCheckerWithPrecision(
  const Vector<T,length_a>& vector_a,
  const Vector<T,length_b>& vector_b,
  bool (*condition_checker)(T, T, T),
  T precision)
{
  if (vector_a.length() != vector_b.length()) return false;
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end())
  {
    if (! condition_checker(*(iter_a++), *(iter_b++), precision)) return false;
  }
  return true;
}

template<typename T, int length_a, int length_b>
bool IsVectorApproximatelyEqual(
  const Vector<T,length_a>& vector_a,
  const Vector<T,length_b>& vector_b,
  const T precision)
{
  return ConditonCheckerWithPrecision(vector_a, vector_b, &IsApproximatelyEqual, precision);
}
//
//template<typename T, int length_a, int length_b>
//bool IsEqual(
//  const Vector<T,length_a>& vector_a,
//  const Vector<T,length_b>& vector_b)
//{
//  return IfAll(vector_a, vector_b, &IsEqual);
//}


//bool IsEqual(const std::vector<Int>& vector_a, const std::vector<Int>& vector_b);

bool IsEqual(const Quaternion& quaternion_a, const Quaternion& quaternion_b);



bool IsEqual(const Point& point_a, const Point& point_b,
             const Real precision = VERY_SMALL);

bool IsEqual(std::vector<Point> points_a, std::vector<Point> points_b);
  
/** Returns true if num is nan */
bool IsNan(Real num);

///** Returns true if num is nan */
//std::vector<bool> IsNan(Vector<Real> input);
//
///** Returns opposite bool as input */
//std::vector<bool> Not(std::vector<bool> input);
//
///** Returns true if all bools are true */
//bool All(std::vector<bool> input);
//
///** Returns true if any one of the bools is true */
//bool Any(std::vector<bool> input);
//
///** Opposite of All: returns true if none of the inputs are true */
//bool None(std::vector<bool> input);

/** Returns true if num is +inf or -inf */
bool IsInf(Real num);

///** Returns true if num is +inf or -inf */
//std::vector<bool> IsInf(Vector<Real> input);

bool ComparisonOpTest();
  
  
  
} // namespace mcl

