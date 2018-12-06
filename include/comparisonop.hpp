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

#include "vector.hpp"
#include <vector>
#include <functional>

namespace mcl
{
////////////////////////////////////////////////////////////////////////////////
// Comparisons on trivial types
////////////////////////////////////////////////////////////////////////////////

inline bool IsEqual(
  const double num_a,
  const double num_b) { return num_a == num_b; }

inline bool IsEqual(
  const float num_a,
  const float num_b) { return num_a == num_b; }

inline bool IsApproximatelyEqual(
  const double num_a,
  const double num_b,
  const double precision = VERY_SMALL) { return std::fabs(num_a - num_b) < precision; }

inline bool IsApproximatelyEqual(
  const float num_a,
  const float num_b,
  const float precision = VERY_SMALLF) { return std::abs(num_a - num_b) < precision; }

inline bool IsApproximatelyEqual(
  const size_t num_a,
  const size_t num_b,
  const size_t /*precision*/) { return num_a == num_b; }

inline bool IsApproximatelyEqual(
  const Int num_a,
  const Int num_b,
  const Int /*precision*/) noexcept { return num_a == num_b; }

inline bool IsApproximatelyEqual(
  const UInt num_a,
  const UInt num_b,
  const UInt /*precision*/) noexcept { return num_a == num_b; }

template<typename T>
bool IsApproximatelyEqual(
  const Complex<T> num_a,
  const Complex<T> num_b,
  const T precision = VERY_SMALL) noexcept { return std::abs(num_a - num_b) < precision; }

template<typename T>
bool IsSmallerOrEqual(
  const T num_a,
  const T num_b) noexcept { return std::islessequal(num_a, num_b); }

template<typename T>
bool IsLargerOrEqual(
  const T num_a,
  const T num_b) noexcept { return std::isgreaterequal(num_a, num_b); }

//bool IsSmallerOrEqual(const Real num_a, const Real num_b, const Real precision) {
//  return num_a <= (num_b + precision);
//}
//
//bool IsLargerOrEqual(const Real num_a, const Real num_b, const Real precision) {
//  return num_a >= (num_b - precision);
//}

/** Returns true if num is +inf or -inf */
template<typename T>
bool IsInf(
  T num) { return std::isinf(num); }

/** Returns true if num is nan */
template<typename T>
bool IsNan(
  T num) noexcept { return isnan(num); }

////////////////////////////////////////////////////////////////////////////////
// Vector checks
////////////////////////////////////////////////////////////////////////////////

/** Checks a condition on all the elements of the vector. Returns true
 if all the conditions are true. False otherwise. */
template<typename T>
bool IsAnyConditionTrue(
  const Vector<T>& vector,
  std::function<bool(
    T)> condition_checker) noexcept
{
  for (size_t i = 0; i < vector.size(); ++i) { if (condition_checker(vector[i])) { return true; } }
  return false;
}

/** Checks a condition on all the elements of the vector. Returns true
 if all the conditions are true. False otherwise. */
template<typename T>
bool AreAllConditionsTrue(
  const Vector<T>& vector,
  std::function<bool(
    T)> condition_checker) noexcept
{
  for (size_t i = 0; i < vector.size(); ++i) { if (! condition_checker(vector[i])) { return false; } }
  return true;
}

template<typename T>
bool AreAllConditionsTrue(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b,
  std::function<bool(
    T,
    T)> condition_checker) noexcept
{
  if (vector_a.size() != vector_b.size()) return false;
  auto iter_a(vector_a.begin());
  auto iter_b(vector_b.begin());
  while (iter_a != vector_a.end()) { if (! condition_checker(*(iter_a++), *(iter_b++))) { return false; } }
  return true;
}

template<typename T>
bool IsEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b) noexcept
{
  return AreAllConditionsTrue<T>
  (
    vector_a,
    vector_b,
    [](
    T a,
    T b)
    {
      return a == b;
    });
}

template<typename T, typename TPrecision>
bool IsApproximatelyEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b,
  const TPrecision precision = VERY_SMALL) noexcept
{
  return AreAllConditionsTrue<T>
  (
    vector_a,
    vector_b,
    [precision](
    T a,
    T b)
    {
      return IsApproximatelyEqual(a, b, precision);
    });
}

//
//template<typename T>
//inline Vector<bool> IsNan(
//  const Vector<T>& input) noexcept
//{
//  Vector<T> output(input.size());
//  ForEach(input, [] (T element) { return IsNan(element); }, output);
//  return std::move(output);
//}
//
//
///** Returns opposite bool as input */
//inline Vector<bool> Not(
//  const Vector<bool>& input) noexcept
//{
//  Vector<bool> output(input.size());
//  ForEach<bool, bool>(
//    input,
//    [] (const bool value) { return !value; },
//    output);
//  return output;
//}
//
//
///** Returns true if all bools are true */
//inline bool AreAllTrue(
//  const Vector<bool>& input) noexcept
//{
//  return AreAllConditionsTrue<bool>(
//    input,
//    [] (const bool element) { return element; });
//}
//
//
///** Opposite of All: returns true if none of the inputs are true */
//inline bool None(
//  Vector<bool> input) noexcept
//{
//  return AreAllTrue(Not(input));
//}
//
//
//
//inline bool Any(
//  const Vector<bool>& input) noexcept
//{
//  return IsAnyConditionTrue<bool>(
//    input,
//    [] (const bool element) { return element; });
//}
//
///** Returns true if num is +inf or -inf */
//template<typename T>
//inline Vector<bool> IsInf(
//  const Vector<T>& input) noexcept
//{
//  Vector<T> output(input.size());
//  ForEach(input, [] (T element) { return IsInf(element); }, output);
//  return std::move(output);
//}

template<typename T>
bool AreAllSmallerOrEqual(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b)
{
  return AreAllConditionsTrue<T>
  (
    vector_a,
    vector_b,
    [](
    T a,
    T b)
    {
      return IsSmallerOrEqual(a, b);
    });
}

/** Returns true if num is nan */
template<typename T>
Vector<bool> IsNan(
  Vector<T> input) noexcept
{
  Vector<bool> output(input.size());
  ForEach
  (input,
   [](
   bool e)
   {
     return IsNan(e);
   },
   output);
  return output;
}

/**
 Returns true if the imaginary part is exactly zero
 (tests equality with the type's 0).
 */
template<typename T>
bool IsReal(
  const Complex<T>& input) noexcept { return input.imag() == T(0.0); }

/**
 Returns true if the imaginary part is approximately zero. The precision used
 is VERY_SMALL in equality operations, hence use only for testing.
 */
template<typename T>
bool IsApproximatelyReal(
  const Complex<T>& input,
  const T precision = VERY_SMALL) noexcept { return IsApproximatelyEqual(input.imag(), T(0.0), precision); }

template<typename T>
bool IsApproximatelyReal(
  const Vector<Complex<T>>& vector,
  const T precision = VERY_SMALL) noexcept
{
  for (auto& element : vector) { if (! IsApproximatelyReal(element, precision)) { return false; } }
  return true;
}
} // namespace mcl
