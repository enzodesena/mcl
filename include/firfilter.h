/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once
#include <vector>

#include "mcltypes.h"
#include "digitalfilter.h"
#include "vectorop.h"

namespace mcl {
/** FIR Filter */
template<typename T>
class FirFilter : public DigitalFilter {
public:
  /** Constructs a default FIR filter, i.e. identical filter */
  FirFilter() noexcept
  {
  }
  
  /** Constructs an FIR filter with impulse response B. */
  FirFilter(
    const Vector<T>& B) noexcept
  {
  }
  
  /** 
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input. 
   As a second example, if B=[0,1], A=[1], you will have 
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  virtual T Filter(
    const T input_sample) noexcept;
  
  using DigitalFilter::Filter;
  
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
  
  }
  
  /** Resets the state of the filter */
  void Reset() noexcept
  {
  
  }
  
  /** Returns the impulse response of the filter */
  Vector<T> impulse_response() noexcept
  {
  
  }
  
  /** Constructs a filter for which output==gain*input always. */
  static FirFilter GainFilter(const T gain) noexcept;
  
  
  virtual ~FirFilter() {}
  
private:
#ifdef MCL_APPLE_ACCELERATE
  T FilterAppleDsp(T input_sample) noexcept;
  void FilterAppleDsp(const T* __restrict input_data, const Int num_samples,
                      T* __restrict output_data) noexcept;
#endif
  
  template<class T>
  void GetExtendedInput(const T* __restrict input_data, const Int num_samples,
                        T* __restrict extended_input_data) {
    
    // Stage 1
    for (Int i=0; i<counter_; ++i) {
      extended_input_data[i] = delay_line_[counter_-i-1];
    }
    
    // Stage 2
    // Starting from counter_ in padded_data
    // Ending in counter_+(length_-counter_-1)
    for (Int i=counter_; i<(length_-1); ++i) {
      extended_input_data[i] = delay_line_[length_-1-(i-counter_)];
    }
    
    // Stage 3
    // Append input signal
    for (Int i=(length_-1); i<(length_-1+num_samples); ++i) {
      extended_input_data[i] = input_data[i-(length_-1)];
    }
  }
  
  T FilterStraight(T input_sample) noexcept;
  
  Vector<T> FilterSequential(const Vector<T>& input) noexcept;
  
  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients() noexcept;
  
  Vector<T> impulse_response_;
  Vector<T> impulse_response_old_;
  Int update_index_;
  Int update_length_;
  
  bool updating_;
  
  /* This is the current vector of coefficients. When the filter is updating
   this will in general be different from impulse_response_. */
  Vector<T> coefficients_;
  Vector<T> delay_line_;
  Int counter_;
  Int length_;
};
  
  
class GainFilter : public DigitalFilter {
public:
  GainFilter(const T gain) : gain_(gain) {}
  
  virtual T Filter(const T input) noexcept {
    return input*gain_;
  }
  
  virtual void Filter(const T* input_data, const Int num_samples,
                      T* output_data) noexcept {
    Multiply(input_data, num_samples, gain_, output_data);
  }
  
  virtual void Reset() {}
private:
  T gain_;
};

class IdenticalFilter : public DigitalFilter {
public:
  IdenticalFilter() {}
  
  virtual T Filter(const T input) noexcept {
    return input;
  }
  
  virtual void Filter(const T* input_data, const Int num_samples,
                      T* output_data) noexcept {
    for (Int i=0; i<num_samples; ++i) {
      output_data[i] = input_data[i];
    }
  }
  
  virtual void Reset() {}
};


/** Tests */
bool FirTest();
void FirSpeedTests();
  
} // namespace mcl

