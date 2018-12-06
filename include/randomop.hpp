/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

// Using TR1 and old c++ library

#include "vector.hpp"
#include <random>

namespace mcl
{
/**
 RandomGenerator class
 */
class RandomGenerator
{
public:
  RandomGenerator() noexcept
    : distribution_norm_(std::normal_distribution<double>(0.0, 1.0))
    , distribution_uniform_(std::uniform_real_distribution<double>(0.0, 1.0))
  {
    generator_.seed(1u);
  }


  RandomGenerator(
    unsigned int seed) noexcept
    : distribution_norm_(std::normal_distribution<double>(0.0, 1.0))
    , distribution_uniform_(std::uniform_real_distribution<double>(0.0, 1.0))
  {
    generator_.seed(seed);
  }


  /** 
   Returns a vector containing pseudorandom values drawn from the standard
   normal distribution. Equivalent to Matlab's randn(size,1);
   */
  Vector<Real> Randn(
    const Int size) noexcept
  {
    Vector<Real> output(size);
    for (Int i = 0; i < size; i++)
    {
      output[i] = distribution_norm_(generator_);
    }
    return output;
  }


  /**
   Returns a vector containing pseudorandom values drawn from the uniform
    distribution. Equivalent to Matlab's rand(size,1);
   */
  Vector<Real> Rand(
    const Int size) noexcept
  {
    Vector<Real> output(size);
    for (Int i = 0; i < size; i++)
    {
      output[i] = distribution_uniform_(generator_);
    }
    return output;
  }


  /** 
   Returns a single pseudorandom value drawn from the uniform
   distribution. Equivalent to Matlab's rand(1,1);
   */
  Real Rand() noexcept
  {
    return Rand(1)[0];
  }


  /**
   Returns a single pseudorandom integer value, uniformly distributed 
   between the given extrema (which are included).
   */
  Int RandInt(
    const Int minimum,
    const Int maximum) noexcept
  {
    // number of possible outcomes
    // e.g. max = 2, min = 0 => num_outcomes = 3
    const Int num_outcomes = maximum - minimum + 1;
    Int output = (Int)floor(
        (double)(distribution_uniform_(generator_) * ((double)num_outcomes))) +
      minimum;

    ASSERT(output >= minimum);
    ASSERT(output <= maximum);
    return output;
  }


  /** Set seed of the for random generator */
  void SetSeed(
    unsigned int seed) noexcept
  {
    generator_.seed(seed);
  }


private:
  std::default_random_engine generator_;
  std::normal_distribution<double> distribution_norm_;
  std::uniform_real_distribution<double> distribution_uniform_;
};
} // namespace mcl
