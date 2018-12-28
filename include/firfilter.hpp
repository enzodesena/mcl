/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "digitalfilter.hpp"
#include "mclintrinsics.hpp"
#include "rampsmoother.h"

namespace mcl
{
/** FIR Filter */
template<typename T>
class FirFilter
{
public:
  /** Constructs a default FIR filter, i.e. identical filter */
  FirFilter() noexcept
    : FirFilter(UnaryVector<T>(1.0))
  {
  }


  /** Constructs an FIR filter with impulse response B. */
  FirFilter(
    const Vector<T>& B,
    size_t max_expected_input_length = 0) noexcept
    : state_(B.size(), T(0.0))
    , numerator_coeffs_(B)
    , counter_(B.size() - 1)
    , length_(B.size())
    , numerator_coeffs_smoothers_(B.size(), T(0.0))
    , temp_conv_data_(Vector(max_expected_input_length + length_ - 1, T(0.0)))
    , temp_add_data_(Vector(max_expected_input_length, T(0.0)))
  {
    for (size_t i=0; i<B.size(); ++i)
    {
      numerator_coeffs_smoothers_[i].SetTargetValue(B[i], 0);
    }
  }

  T Filter(
    const T input_sample) noexcept
  {
    if (numerator_coeffs_smoothers_[0].IsUpdating())
    {
      UpdateCoefficients();
    }
    return FilterSample(input_sample);
  }

  /**
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input.
   As a second example, if B=[0,1], A=[1], you will have
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  T FilterSample(
    const T input_sample) noexcept
  {
    if (length_ == 1)
    {
      state_[0] = input_sample;
      return input_sample * numerator_coeffs_[0];
    }
    state_[counter_] = input_sample;
    T result = T(0.0);
    size_t index = counter_;
    for (size_t i = 0; i < length_; ++i)
    {
      result += numerator_coeffs_[i] * state_[index++];
      if (index >= length_)
      {
        index = 0;
      }
    }
    counter_ = (counter_) ? counter_-1 : length_ - 1;
    return result;
  }


  /** 
   Updates the filter coefficients. You can set how long it takes to 
   update the coefficients (using linear interpolation between old and new
   impulse response). If an update is requested while another is already in
   progress, the new interpolation will pick up from where the old one left
   off to avoid audible artifacts.
   @param[in] update_length How many calls to Filter it takes to update the
   coefficients. A value of 0 means that the update is instantaneous. A call
   to Filter(const T input) counts as one, just like
   Filter(const Vector<T>& input).
   */
  void SetNumeratorCoeffs(
    const Vector<T>& impulse_response,
    const size_t update_length = 0) noexcept
  {
    ASSERT(impulse_response.size() == length_);
    
    for (size_t i=0; i<numerator_coeffs_.size(); ++i)
    {
      numerator_coeffs_smoothers_[i].SetTargetValue(impulse_response[i], update_length);
    }
  }


  /** Resets the state of the filter */
  void SetStateToZero() noexcept
  {
    SetToZero(state_);
  }


  /** Returns the impulse response of the filter */
  Vector<T> impulse_response() const noexcept
  {
    return numerator_coeffs_;
  }


  void FilterSerial(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    for (size_t i=0; i<input.size(); ++i) {
      output[i] = FilterSample(input[i]);
    }
  }


  void FilterAdd(
    const Vector<T>& input_to_filter,
    const Vector<T>& input_to_add,
    Vector<T>& output) noexcept
  {
    if (temp_add_data_.size() < input_to_filter.size())
    {
      temp_add_data_ = Vector<T>(input_to_filter.size());
    }
    Vector<T> temp = MakeReference(temp_add_data_, 0, input_to_filter.size());
    Filter(input_to_filter, temp);
    Add(temp, input_to_add, output);
  }
  
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    ASSERT(input.size() == output.size());
    
    if (numerator_coeffs_smoothers_[0].IsUpdating())
    {
      UpdateCoefficients();
    }
    
    if (input.size() < length_)
    {
      FilterSerial(input, output);
    }
    else
    {
      if (temp_conv_data_.size() < (input.size() + length_ - 1))
      {
        temp_conv_data_ = Vector<T>(input.size() + length_ - 1);
      }
      GetExtendedInput(input, temp_conv_data_);
      Conv(temp_conv_data_, numerator_coeffs_, output);
      // Reorganise state for the next run
      for (size_t i = 0; i < length_; ++i)
      {
        state_[i] = input[input.size() - 1 - i];
      }
      counter_ = length_ - 1;
    }
  }
private:
  /* This is the current vector of coefficients. When the filter is updating
   this will in general be different from target_coefficients_. */
  Vector<T> state_;
  Vector<T> numerator_coeffs_;
  size_t counter_;
  size_t length_;
  Vector<RampSmoother<T>> numerator_coeffs_smoothers_;
  Vector<T> temp_conv_data_;
  Vector<T> temp_add_data_;

  void GetExtendedInput(
    const Vector<T>& input,
    Vector<T>& extended_input) noexcept
  {
    // Stage 1
    for (size_t i = 0; i < counter_; ++i)
    {
      extended_input[i] = state_[counter_ - i - 1];
    }

    // Stage 2
    // Starting from counter_ in padded_data
    // Ending in counter_+(length_-counter_-1)
    for (size_t i = counter_; i < (length_ - 1); ++i)
    {
      extended_input[i] = state_[length_ - 1 - (i - counter_)];
    }

    // Stage 3
    // Append input signal
    for (size_t i = (length_ - 1); i < (length_ - 1 + input.size()); ++i)
    {
      extended_input[i] = input[i - (length_ - 1)];
    }
  }


  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients() noexcept
  {
    for (size_t i=0; i<numerator_coeffs_.size(); ++i)
    {
      numerator_coeffs_[i] = numerator_coeffs_smoothers_[i].GetNextValue();
    }
  }

};
} // namespace mcl
