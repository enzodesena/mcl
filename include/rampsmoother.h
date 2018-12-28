/*
 MCL
 Copyright (c) 2018, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

namespace mcl
{

template<typename T>
class RampSmoother
{
public:
  /**
   @param[in] initial_value The initial assigned value. */
  RampSmoother(
    const T initial_value) noexcept;

  T GetNextValue() noexcept;

  T GetNextValue(
    const size_t num_jumps) noexcept;


  /** Takes an array of values (`input_data`), multiplies them by the
   next values coming out of the smoother, and writes the result into
   an output array (`output_data`).
   @param[in] input The input vector.
   @param[out] output The output vector (data will be overwritten. */
  void GetNextValuesMultiply(
    const mcl::Vector<T>& input,
    mcl::Vector<T>& output) noexcept;

  /** Takes an array of values (`input_data`), multiplies them by the
   next values coming out of the smoother, and adds the result to
   an input-output array (`input_output_data`).
   @param[in] input_data The input data to multiply by.
   @param[in,out] input_output_data The data onto which we will add the result
   of the multiplication. */
  void GetNextValuesMultiplyAdd(
    const T* input_data,
    const Int num_samples,
    T* input_output_data) noexcept;

  T target_value() const noexcept;

  /**
    @param[in] target_value the target value
    @param[in] ramp_samples the number of calls to `GetNextValue()` it will
    take to reach the target value. The number of calls can be fractional.
    */
  void SetTargetValue(
    const T target_value,
    const double ramp_samples) noexcept;

  bool IsUpdating() const noexcept;

private:
  T current_value_;
  T target_value_;
  T step_;
  Int countdown_;
};

} // namespace

#include "rampsmoother.ipp"
