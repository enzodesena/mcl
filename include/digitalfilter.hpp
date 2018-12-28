/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mclintrinsics.hpp"
#include "rampsmoother.h"

namespace mcl
{
/** FIR Filter */
template<typename T>
class DigitalFilter
{
public:
  /** Constructs a default FIR filter, i.e. identical filter */
  DigitalFilter() noexcept;
  
  DigitalFilter(
    const Vector<T>& numerator_coeffs,
    const Vector<T>& denominator_coeffs = mcl::UnaryVector<T>(T(1.0)),
    size_t max_expected_input_length = 0) noexcept;

  T Filter(
    const T input_sample) noexcept;

  /**
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input.
   As a second example, if B=[0,1], A=[1], you will have
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  T FilterSample(
    const T input_sample) noexcept;

  Vector<T> GetNumeratorCoeffs() noexcept;
  
  Vector<T> GetDenominatorCoeffs() noexcept;

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
    const double update_length = 0) noexcept;

  void SetDenominatorCoeffs(
    const Vector<T>& denominator_coeffs,
    const double update_length = 0) noexcept;

  /** Resets the state of the filter */
  void SetStateToZero() noexcept;

  void FilterAdd(
    const Vector<T>& input_to_filter,
    const Vector<T>& input_to_add,
    Vector<T>& output) noexcept;
  
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept;
  
  Vector<T> Filter(
    const Vector<T>& input) noexcept;
private:
  /* This is the current vector of coefficients. When the filter is updating
   this will in general be different from target_coefficients_. */
  Vector<T> state_;
  Vector<T> numerator_coeffs_;
  Vector<T> denominator_coeffs_;
  size_t counter_;
  size_t length_;
  Vector<RampSmoother<T>> numerator_smoothers_;
  Vector<T> temp_conv_data_;
  Vector<T> temp_add_data_;

  void GetExtendedInput(
    const Vector<T>& input,
    Vector<T>& extended_input) noexcept;

  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients() noexcept;
  

  void FilterSerial(
    const Vector<T>& input,
    Vector<T>& output) noexcept;
};
} // namespace mcl


#include "digitalfilter.ipp"
