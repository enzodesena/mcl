/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#pragma once

#ifndef VERY_SMALL
  #define VERY_SMALL (0.0000000000001)
#endif

#ifndef VERY_SMALLF
  #define VERY_SMALLF (0.0000000000001f)
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
#include <functional>

namespace mcl {

// Forward declarations
template<typename T>
class Quaternion;
template<typename T>
struct Point;

template<typename T, typename U>
inline void ForEach(
  const Vector<T>& input_vector,
  std::function<U(T)> operation,
  Vector<U>& output_vector) noexcept;
// End of formard declarations

////////////////////////////////////////////////////////////////////////////////
// Comparisons on trivial types
////////////////////////////////////////////////////////////////////////////////

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
  const float precision = VERY_SMALLF)
{
  return std::abs(num_a - num_b) < precision;
}

inline bool IsApproximatelyEqual(
  const size_t num_a,
  const size_t num_b,
  const size_t /*precision*/)
{
  return num_a == num_b;
}

inline bool IsApproximatelyEqual(
  const Int num_a,
  const Int num_b,
  const Int /*precision*/) noexcept
{
  return num_a == num_b;
}

inline bool IsApproximatelyEqual(
  const UInt num_a,
  const UInt num_b,
  const UInt /*precision*/) noexcept
{
  return num_a == num_b;
}

template<typename T>
inline bool IsApproximatelyEqual(
  const Complex<T> num_a,
  const Complex<T> num_b,
  const T precision = VERY_SMALL) noexcept
{
  return std::abs(num_a - num_b) < precision;
}

template<typename T>
inline bool IsSmallerOrEqual(
  const T num_a,
  const T num_b) noexcept
{
  return std::islessequal(num_a, num_b);
}

template<typename T>
inline bool IsLargerOrEqual(
  const T num_a,
  const T num_b) noexcept
{
  return std::isgreaterequal(num_a, num_b);
}

//bool IsSmallerOrEqual(const Real num_a, const Real num_b, const Real precision) {
//  return num_a <= (num_b + precision);
//}
//
//bool IsLargerOrEqual(const Real num_a, const Real num_b, const Real precision) {
//  return num_a >= (num_b - precision);
//}


/** Returns true if num is +inf or -inf */
template<typename T>
inline bool IsInf(T num)
{
  return std::isinf(num);
}

/** Returns true if num is nan */
template<typename T>
inline bool IsNan(
  T num) noexcept
{
  return isnan(num);
}

////////////////////////////////////////////////////////////////////////////////
// Vector checks
////////////////////////////////////////////////////////////////////////////////

/** Checks a condition on all the elements of the vector. Returns true
 if all the conditions are true. False otherwise. */
template<typename T>
inline bool IsAnyConditionTrue(
  const Vector<T>& vector,
  std::function<bool(T)> condition_checker) noexcept
{
  for (const auto& element : vector)
  {
    if (condition_checker(element))
    {
      return true;
    }
  }
  return false;
}


/** Checks a condition on all the elements of the vector. Returns true
 if all the conditions are true. False otherwise. */
template<typename T>
inline bool AreAllConditionsTrue(
  const Vector<T>& vector,
  std::function<bool(T)> condition_checker) noexcept
{
  for (const auto& element : vector)
  {
    if (! condition_checker(element))
    {
      return false;
    }
  }
  return true;
}

template<typename T>
inline bool AreAllConditionsTrue(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b,
  std::function<bool(T,T)> condition_checker) noexcept
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


template<typename T>
inline bool IsEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  return AreAllConditionsTrue<T>(
    vector_a,
    vector_b,
    [] (T a, T b) { return a == b; });
}


template<typename T, typename TPrecision>
inline bool IsApproximatelyEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b,
  const TPrecision precision = VERY_SMALL) noexcept
{
  return AreAllConditionsTrue<T>(
    vector_a,
    vector_b,
    [precision] (T a, T b) { return IsApproximatelyEqual(a, b, precision); });
}


template<typename T>
inline Vector<bool> IsNan(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  ForEach(input, [] (T element) { return IsNan(element); }, output);
  return std::move(output);
}


/** Returns opposite bool as input */
inline Vector<bool> Not(
  const Vector<bool>& input) noexcept
{
  Vector<bool> output(input.length());
  ForEach<bool, bool>(input, [] (bool value) { return !value; }, output);
  return output;
}


/** Returns true if all bools are true */
inline bool AreAllTrue(
  const Vector<bool>& input) noexcept
{
  return AreAllConditionsTrue<bool>(
    input,
    [] (bool element) { return element; });
}



inline bool None(
  Vector<bool> input) noexcept
{
  return AreAllTrue(Not(input));
}



inline bool Any(
  const Vector<bool>& input) noexcept
{
  return IsAnyConditionTrue<bool>(input, [] (bool element) { return element; });
}


template<typename T>
inline Vector<bool> IsInf(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.length());
  ForEach(input, [] (T element) { return IsInf(element); }, output);
  return std::move(output);
}

template<typename T>
bool AreAllSmallerOrEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b)
{
  return AreAllConditionsTrue<T>(
    vector_a,
    vector_b,
    [] (T a, T b) { return IsSmallerOrEqual(a, b); });
}



template<typename T>
inline bool IsEqual(
  const Quaternion<T>& q_a,
  const Quaternion<T>& q_b) noexcept
{
  return q_a.w() == q_b.w() && q_a.x() == q_b.x() &
    q_a.y() == q_b.y() && q_a.z() == q_b.z();
}


template<typename T>
inline bool IsApproximatelyEqual(
  const Point<T>& point_a,
  const Point<T>& point_b,
  const T precision = VERY_SMALL) noexcept
{
  return
    mcl::IsApproximatelyEqual(point_a.x(), point_b.x(), precision) &&
    mcl::IsApproximatelyEqual(point_a.y(), point_b.y(), precision) &&
    mcl::IsApproximatelyEqual(point_a.z(), point_b.z(), precision);
}

template<typename T>
inline bool IsEqual(
  const Point<T>& point_a,
  const Point<T>& point_b)
{
  return
    point_a.x() == point_b.x() &&
    point_a.y() == point_b.y() &&
    point_a.z() == point_b.z();
}


template<typename T>
inline bool IsEqual(
  const Vector<Point<T>>& points_a,
  const Vector<Point<T>>& points_b) noexcept
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
inline Vector<bool> IsNan(Vector<T> input) noexcept
 {
  Vector<bool> output(input.length());
  for (size_t i=0; i<input.length(); ++i)
  {
    output[i] = IsNan(input[i]);
  }
  return output;
}



/** Returns true if any one of the bools is true */

inline bool Any(
  Vector<bool> input) noexcept
{
  return AreAllConditionsTrue<bool>(input, [] (bool element) { return !element; });
}

///** Opposite of All: returns true if none of the inputs are true */
//bool None(Vector<bool> input);


///** Returns true if num is +inf or -inf */
//Vector<bool> IsInf(Vector<Real> input);


template<typename T>
inline bool IsApproximatelyEqual(
  const Quaternion<T>& quat_a,
  const Quaternion<T>& quat_b,
  const T precision = VERY_SMALL)
{
  return
    IsApproximatelyEqual(quat_a.w(), quat_b.w()) &&
    IsApproximatelyEqual(quat_a.x(), quat_b.x()) &&
    IsApproximatelyEqual(quat_a.y(), quat_b.y()) &&
    IsApproximatelyEqual(quat_a.z(), quat_b.z());
}

bool ComparisonOpTest();
  
  
  
} // namespace mcl

