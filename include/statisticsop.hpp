/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vectorop.hpp"
#include "mclintrinsics.hpp"
#include "matrixop.hpp"

namespace mcl
{
template<typename T>
/** Equivalent to Matlab's mean(input) */
T Mean(
  const Vector<T>& input) noexcept { return Sum(input) / ((T)input.size()); }

/**
 Returns the geometric mean of the input vector. Equivalent
 to Matlab's geomean(input)
 **/
template<typename T>
T Geomean(
  const Vector<T>& input) noexcept
{
  // TODO: Throw error for negative entries
  return Pow(Prod(input), 1.0 / ((T)input.size()));
}

/**
 Weighted mean. Not implemented in Matlab (but should be). The weights are
 normalised inside the function. Hence Mean(input, ones(N)) gives the same
 result as Mean(input, ones(N)/N).
 */
template<typename T>
T Mean(
  const Vector<T>& input,
  const Vector<T>& weights) noexcept
{
  ASSERT(input.size() == weights.size());
  ASSERT(IsNonNegative(weights));

  // Normalise the weigths
  Vector<T> normalised_weights = Multiply(weights, 1.0 / Sum(weights));
  ASSERT(Sum(normalised_weights) == 1.0);
  return Sum(Multiply(input, normalised_weights));
}

/**
 Returns the standard deviation of the `input` vector. Equivalent to Matlab's
 std(input). This includes the correction for having an unbiased estimator.
 */
template<typename T>
T Std(
  const Vector<T>& input) noexcept { return sqrt(Var(input)); }

/** Var (unbiased estimator) */
template<typename T>
T Var(
  const Vector<T>& input) noexcept
{
  T mean = Mean(input);
  T output(0.0);
  for (size_t i = 0; i < input.size(); ++i) { output += pow(input[i] - mean, 2.0); }
  return output / ((T)(input.size() - 1));
}

/** Weighted var (biased estimator) */
template<typename T>
T Var(
  const Vector<T>& input,
  const Vector<T>& weights) noexcept
{
  ASSERT(IsNonNegative(weights));

  T weighted_mean = Mean(input, weights);
  Vector<T> tt = AddScalar(input, -weighted_mean);
  Vector<T> temp = Pow(tt, 2.0);
  Vector<T> norm_weights = Multiply<T>(weights, 1.0 / Sum(weights));

  return (Sum(Multiply(norm_weights, temp)));
}

/** Returns the Pearson linear correlation between `vector_a` and `vector_b` */
template<typename T>
T Corr(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  T pearson_num_lin = Sum
  (mcl::Multiply
    (
      AddScalar(x, -Mean(x)),
      AddScalar(y, -Mean(y))));
  T pearson_den_lin = sqrt(Sum(Pow(AddScalar(x, -Mean(x)), 2.0))) *
    sqrt(Sum(Pow(AddScalar(y, -Mean(y)), 2.0)));
  return pearson_num_lin / pearson_den_lin;
}

/**
 Calculates the entropy of a discreate random variable with given `pdf'.
 It normalises the pdf if its sum is not 1.
 Note: this function is identical to Matlab's only for uint8 values.
 */
template<typename T>
T Entropy(
  const Vector<T> pdf,
  T base) noexcept
{
  Vector<T> normalised_pdf(Multiply(pdf, 1.0 / Sum(pdf)));
  return -Sum(Multiply(pdf, Log(pdf))) / log(base);
}

template<typename T>
Matrix<T> Cov(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  Vector<Vector<T>> input(2);
  input[0] = x;
  input[1] = y;
  return Cov(input);
}

template<typename T>
Matrix<T> Cov(
  const Vector<Vector<T>>& input) noexcept
{
  const size_t N = input.size();
  Matrix<T> output(N, N);
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j) { output.SetElement(i, j, CovElement(input[i], input[j])); }
  }
  return output;
}

template<typename T>
T CovElement(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  ASSERT(x.size() == y.size());
  const size_t N = x.size();

  T output = Sum(Multiply(AddScalar(x, -Mean(x)), AddScalar(y, -Mean(y))));
  // In case N>1 use the unbiased estimator of covariance.
  output = (N > 1) ? output / ((T)(N - 1)) : output / ((T)(N));
  return output;
}
} // namespace MCL
