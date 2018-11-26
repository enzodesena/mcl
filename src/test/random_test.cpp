/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "randomop.hpp"
#include "comparisonop.hpp"
#include "mcltypes.hpp"
#include <vector>
#include <cassert>
#include "vectorop.hpp"

namespace mcl {

bool RandomGeneratorTest() {
  
  RandomGenerator rand_gen;

  Vector<Real> rand_vector_a = rand_gen.Randn(5);
  ASSERT(rand_vector_a.length() == 5);
  // Check whether seeds are different:
  Vector<Real> rand_vector_b = rand_gen.Randn(5);
  ASSERT(!IsEqual(rand_vector_a, rand_vector_b));
  
  ASSERT(Abs(Mean(rand_gen.Randn(100000))) < 0.1);
  ASSERT(Abs(Std(rand_gen.Randn(100000))-1.0) < 0.1);
  
  // Testing uniform distribution
  ASSERT(! IsEqual(rand_gen.Rand(1), rand_gen.Rand(1)));
  ASSERT(! IsEqual(rand_gen.Rand(1), rand_gen.Rand(1)));
  
  // Test that output is between 0 and 1
  ASSERT(Abs(rand_gen.Rand(1)[0]-0.5) < 0.5);
  ASSERT(Abs(rand_gen.Rand(1)[0]-0.5) < 0.5);
  ASSERT(Abs(rand_gen.Rand(1)[0]-0.5) < 0.5);
  ASSERT(Abs(rand_gen.Rand(1)[0]-0.5) < 0.5);
  ASSERT(Abs(rand_gen.Rand(1)[0]-0.5) < 0.5);
  
  ASSERT(rand_gen.Rand(5).length() == 5);
  ASSERT(Abs(Mean(rand_gen.Rand(100000))-0.5)<0.05);
  
  // Test integer generator
  const Int num_samples = 1000000;
  Vector<Int> rand_int_vector(num_samples);
  const Int min_value = -2;
  const Int max_value = 3;
  Vector<Int> num_occurrances(max_value-min_value+1, 0);
  for (Int i=0; i<num_samples; ++i) {
    Int output = rand_gen.RandInt(min_value, max_value);
    rand_int_vector[i] = output;
    num_occurrances[output-min_value]++;
  }
  ASSERT(mcl::Min(rand_int_vector) >= min_value);
  ASSERT(mcl::Max(rand_int_vector) <= max_value);
  ASSERT(mcl::IsApproximatelyEqual(((Real) mcl::Min(num_occurrances)) / ((Real) num_samples),
                      ((Real) mcl::Max(num_occurrances)) / ((Real) num_samples),
                      1.0E-2));
  
  // Test integer generator
  const Int num_samples_b = 1000000;
  Vector<Int> rand_int_vector_b(num_samples_b);
  const Int min_value_b = 0;
  const Int max_value_b = 1;
  Vector<Int> num_occurrances_b(max_value_b-min_value_b+1, 0);
  for (size_t i=0; i<num_samples_b; ++i) {
    Int output = rand_gen.RandInt(min_value_b, max_value_b);
    rand_int_vector_b[i] = output;
    num_occurrances_b[output-min_value_b]++;
  }
  ASSERT(mcl::Min(rand_int_vector_b) >= min_value_b);
  ASSERT(mcl::Max(rand_int_vector_b) <= max_value_b);
  ASSERT(mcl::IsApproximatelyEqual(((Real) mcl::Min(num_occurrances_b)) / ((Real) num_samples_b),
                      ((Real) mcl::Max(num_occurrances_b)) / ((Real) num_samples_b),
                      1.0E-2));
  
  return true;
}

} // namespace mcl
