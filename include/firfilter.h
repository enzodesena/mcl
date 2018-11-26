/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once
#include <vector>

#include "mcltypes.h"
#include "mclintrinsics.h"
#include "digitalfilter.h"
#include "vectorop.h"

namespace mcl {
/** FIR Filter */
template<typename T>
class FirFilter {
private:

  /* This is the current vector of coefficients. When the filter is updating
   this will in general be different from impulse_response_. */
  Vector<T> delay_line_;
  Vector<T> coefficients_;
  size_t counter_;
  size_t length_;
  Vector<T> impulse_response_;
  Vector<T> impulse_response_old_;
  size_t update_index_;
  size_t update_length_;
  bool updating_;
  
  
  void GetExtendedInput(
    const Vector<T>& input,
    Vector<T> extended_input) noexcept
  {
    
    // Stage 1
    for (size_t i=0; i<counter_; ++i)
    {
      extended_input[i] = delay_line_[counter_-i-1];
    }
    
    // Stage 2
    // Starting from counter_ in padded_data
    // Ending in counter_+(length_-counter_-1)
    for (size_t i=counter_; i<(length_-1); ++i)
    {
      extended_input[i] = delay_line_[length_-1-(i-counter_)];
    }
    
    // Stage 3
    // Append input signal
    for (size_t i=(length_-1); i<(length_-1+input.length()); ++i)
    {
      extended_input[i] = input[i-(length_-1)];
    }
  }
  
  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients() noexcept
  {
    ASSERT(update_index_>=0 && update_index_<=update_length_);
    ASSERT(impulse_response_.length() == impulse_response_old_.length());
    ASSERT(impulse_response_.length() == coefficients_.length());
    T weight_new = ((T)update_index_+T(1))/((T)update_length_+T(1));
    T weight_old = T(1)-weight_new;
    Multiply(impulse_response_, weight_new, coefficients_);
    MultiplyAdd(impulse_response_old_, weight_old,
                coefficients_,
                coefficients_);
    // The above is a lock-free equivalent version to
    // coefficients_ = mcl::Add(mcl::Multiply(impulse_response_, weight_new),
    //                          mcl::Multiply(impulse_response_old_, weight_old));
    
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
    : delay_line_(1)
    , coefficients_(mcl::UnaryVector<T>(1.0))
    , impulse_response_(mcl::UnaryVector<T>(1.0))
    , impulse_response_old_(mcl::UnaryVector<T>(1.0))
    , update_index_(0)
    , update_length_(0)
    , updating_(false)
    , counter_(0)
    , length_(1)
  {
    SetToZero(delay_line_);
  }
  
  /** Constructs an FIR filter with impulse response B. */
  FirFilter(
    const Vector<T>& B) noexcept
    : delay_line_(B.length())
    , coefficients_(B)
    , impulse_response_(B)
    , impulse_response_old_(B)
    , update_index_(0)
    , update_length_(0)
    , updating_(false)
    , counter_(B.length()-1)
    , length_(B.length())
  {
    SetToZero(delay_line_);
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
      return input_sample*coefficients_[0];
    }
    delay_line_[counter_] = input_sample;
    T result = T(0.0);
    size_t index = counter_;
    for (size_t i=0; i<length_; ++i)
    {
      result += coefficients_[i] * delay_line_[index++];
      if (index >= length_)
      {
        index = 0;
      }
    }
    
    if (--counter_ < 0)
    {
      counter_ = length_-1;
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
   to Filter(const T input) counts one, just like
   Filter(const Vector<T>& input).
   */
  void SetImpulseResponse(
    const Vector<T>& impulse_response,
    const Int update_length = 0) noexcept
  {
    if (mcl::IsEqual(impulse_response, impulse_response_))
    {
      return;
    }
  
    if ((Int)impulse_response.length() != length_)
    {
      // If the impulse response changes length, then reset everything.
      length_ = impulse_response.length();
      delay_line_.assign(length_, 0.0);
      counter_ = length_-1;
      impulse_response_ = impulse_response;
      impulse_response_old_ = impulse_response;
      coefficients_ = impulse_response;
      updating_ = false;
      update_length_ = 0;
      update_index_ = 0;
    }
    else
    {
      updating_ = true;
      update_length_ = update_length;
      update_index_ = 0;
      impulse_response_ = impulse_response;
      
      if (! updating_)
      { // If there is no update being carried out
        impulse_response_old_ = impulse_response_;
      }
      else
      {
        impulse_response_old_ = coefficients_;
      }
    }
    ASSERT(impulse_response_.length() == impulse_response_old_.length());
  }
  
  /** Resets the state of the filter */
  void Reset() noexcept
  {
    SetToZero(delay_line_);
  }
  
  /** Returns the impulse response of the filter */
  Vector<T> impulse_response() const noexcept
  {
    return impulse_response_;
  }
  
  
  void FilterSerial(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    auto input_iter(input.begin());
    auto output_iter(output.begin());
    while (output_iter != output.end())
    {
      *(output_iter++) = Filter(*(input_iter++));
    }
  }
  
  
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    ASSERT(input.length() == output.length());
    const size_t num_samples = input.length();
    if (updating_)
    {
      UpdateCoefficients();
    }
    if (length_ == 1)
    {
      delay_line_[0] = input[input.length()-1];
      Multiply(input, coefficients_, output);
      return;
    }
    if (num_samples < length_ || (num_samples+length_-1) > MCL_MAX_VLA_LENGTH)
    {
      FilterSerial(input, output);
      return;
    }
    else
    {
      Vector<T> extended_input(num_samples+length_-1, 0.0);
      GetExtendedInput(input, extended_input);
      Conv(extended_input, coefficients_, output);
      // Reorganise state for the next run
      for (size_t i=0; i<length_; ++i)
      {
        delay_line_[i] = input[num_samples-1-i];
      }
      counter_ = length_-1;
    }
  }
};


  



/** Tests */
bool FirFilterTest();
void FirSpeedTests();
  
} // namespace mcl

