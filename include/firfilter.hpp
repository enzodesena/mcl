/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "digitalfilter.hpp"
#include "mclintrinsics.hpp"

namespace mcl
{
/** FIR Filter */
template<typename T>
class FirFilter : private DigitalFilterInterface<T>
{
private:
  /* This is the current vector of coefficients. When the filter is updating
   this will in general be different from target_coefficients_. */
  Vector<T> delay_line_;
  Vector<T> current_coefficients_;
  size_t counter_;
  size_t length_;
  Vector<T> target_coefficients_;
  Vector<T> starting_coefficients_;
  size_t update_index_;
  size_t update_length_;
  bool updating_;
  Vector<T> temp_conv_data_;
  Vector<T> temp_add_data_;

  void GetExtendedInput(
    const Vector<T>& input,
    Vector<T>& extended_input) noexcept
  {
    // Stage 1
    for (size_t i = 0; i < counter_; ++i)
    {
      extended_input[i] = delay_line_[counter_ - i - 1];
    }

    // Stage 2
    // Starting from counter_ in padded_data
    // Ending in counter_+(length_-counter_-1)
    for (size_t i = counter_; i < (length_ - 1); ++i)
    {
      extended_input[i] = delay_line_[length_ - 1 - (i - counter_)];
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
    ASSERT(update_index_>=0 && update_index_<=update_length_);
    ASSERT(target_coefficients_.size() == starting_coefficients_.size());
    ASSERT(target_coefficients_.size() == current_coefficients_.size());
    
    if (update_length_ == 0)
    {
      current_coefficients_ = target_coefficients_;
    }
    else
    {
      T weight_new = ((T)update_index_ + T(1)) / ((T)update_length_ + T(1));
      T weight_old = T(1) - weight_new;
      std::function<T(T,T)> operation =
        [weight_old, weight_new] (T start_value, T target_value)
        {
          return start_value*weight_old + target_value*weight_new;
        };
      ForEach(starting_coefficients_, target_coefficients_, operation, current_coefficients_);
    }
    

    if (update_index_ == update_length_)
    {
      updating_ = false;
    }
    else
    {
      update_index_++;
    }
  }


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
    : delay_line_(B.size(), T(0.0))
    , current_coefficients_(B)
    , counter_(B.size() - 1)
    , length_(B.size())
    , target_coefficients_(B)
    , starting_coefficients_(B)
    , update_index_(0)
    , update_length_(0)
    , updating_(false)
    , temp_conv_data_(Vector(max_expected_input_length + length_ - 1, T(0.0)))
    , temp_add_data_(Vector(max_expected_input_length, T(0.0)))
  {
  }


  /**
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input.
   As a second example, if B=[0,1], A=[1], you will have
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  T Filter(
    const T input_sample) noexcept
  {
    if (updating_)
    {
      UpdateCoefficients();
    }
    if (length_ == 1)
    {
      delay_line_[0] = input_sample;
      return input_sample * current_coefficients_[0];
    }
    delay_line_[counter_] = input_sample;
    T result = T(0.0);
    size_t index = counter_;
    for (size_t i = 0; i < length_; ++i)
    {
      result += current_coefficients_[i] * delay_line_[index++];
      if (index >= length_)
      {
        index = 0;
      }
    }
    if (counter_-- == 0)
    {
      counter_ = length_ - 1;
    }
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
  void SetImpulseResponse(
    const Vector<T>& impulse_response,
    const size_t update_length = 0) noexcept
  {
    if (mcl::IsEqual(impulse_response, target_coefficients_))
    {
      return;
    }

    if (impulse_response.size() != length_)
    {
      // If the impulse response changes length, then reinitialise.
      length_ = impulse_response.size();
      delay_line_ = Vector<T>(length_, T(0.0));
      counter_ = length_ - 1;
      target_coefficients_ = impulse_response;
      starting_coefficients_ = impulse_response;
      current_coefficients_ = impulse_response;
      updating_ = false;
      update_length_ = 0;
      update_index_ = 0;
    }
    else
    {
      updating_ = true;
      update_length_ = update_length;
      update_index_ = 0;
      target_coefficients_ = impulse_response;
      starting_coefficients_ = current_coefficients_;
    }
    ASSERT(target_coefficients_.size() == starting_coefficients_.size());
  }


  /** Resets the state of the filter */
  void Reset() noexcept override
  {
    SetToZero(delay_line_);
  }


  /** Returns the impulse response of the filter */
  Vector<T> impulse_response() const noexcept
  {
    return target_coefficients_;
  }


  void FilterSerial(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    std::function<T(T)> operation = [this](T value) -> T { return Filter(value); };
    ForEach(input.begin(), input.end(), operation, output.begin());
  }


  void FilterAdd(
    const Vector<T>& input_to_filter,
    const Vector<T>& input_to_add,
    Vector<T>& output) noexcept override
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
    Vector<T>& output) noexcept override
  {
    ASSERT(input.size() == output.size());
    
    if (updating_)
    {
      UpdateCoefficients();
    }
    
    if (length_ == 1)
    {
      Multiply(input, current_coefficients_, output);
      delay_line_[0] = input[input.size() - 1];
    }
    else if (input.size() < length_)
    {
      FilterSerial(input, output);
    }
    else
    {
      if (temp_conv_data_.size() < input.size() + length_ - 1)
      {
        temp_conv_data_ = Vector<T>(input.size() + length_ - 1);
      }
      GetExtendedInput(input, temp_conv_data_);
      Conv(temp_conv_data_, current_coefficients_, output);
      // Reorganise state for the next run
      for (size_t i = 0; i < length_; ++i)
      {
        delay_line_[i] = input[input.size() - 1 - i];
      }
      counter_ = length_ - 1;
    }
  }
};
} // namespace mcl
