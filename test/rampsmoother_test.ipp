/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "rampsmoother.hpp"

namespace mcl
{
inline bool RampSmootherTest()
{
  RampSmoother smoother(1.0);
  assert(!smoother.IsUpdating());
  assert(IsEqual(smoother.GetNextValue(), 1.0));
  assert(IsEqual(smoother.GetNextValue(), 1.0));
  assert(IsEqual(smoother.GetNextValue(), 1.0));
  assert(!smoother.IsUpdating());

  smoother.SetTargetValue(3.0, 2.0); // ramp time
  assert(smoother.IsUpdating());
  assert(IsEqual(smoother.GetNextValue(), 2.0));
  assert(smoother.IsUpdating());
  assert(IsEqual(smoother.GetNextValue(), 3.0));
  assert(!smoother.IsUpdating());
  assert(IsEqual(smoother.GetNextValue(), 3.0));

  smoother.SetTargetValue(1.0, 0.0); // Force to value
  assert(!smoother.IsUpdating());
  assert(IsEqual(smoother.GetNextValue(), 1.0));
  assert(!smoother.IsUpdating());

  Vector<Real> input_samples = {1.0, 1.0, 1.0, 1.0};
  smoother.SetTargetValue(2.0, 2.0); // ramp tim
  Vector<Real> output_samples(4);
  Vector<Real> output_samples_cmp = {1.5, 2.0, 2.0, 2.0};
  smoother.GetNextValuesMultiply(input_samples, output_samples);
  ASSERT(IsEqual(output_samples, output_samples_cmp));

  smoother.SetTargetValue(15000.0, 6.0); // ramp time
  RampSmoother smoother_copy(smoother);
  AVOID_UNUSED_WARNING(smoother_copy);
  ASSERT(IsEqual(smoother.GetNextValue(), smoother_copy.GetNextValue(1)));
  smoother.GetNextValue();
  smoother.GetNextValue();
  ASSERT(IsEqual(smoother.GetNextValue(), smoother_copy.GetNextValue(3)));
  
  return true;
}
} // namespace
