/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mclintrinsics.hpp"
#include "rampsmoother.hpp"

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
    Copies the numerator and denominator coefficients from another filter.
    If the filter order is the same, this also leaves the state intact.
    If desired, the coefficients can be update smoothly.
    @param[in] other_filter the other filter
    @param[in] update_length how many calls to `Filter` it takes to update
    the coefficients
  */
  void CopyCoefficientsFrom(
    const DigitalFilter<T>& other_filter,
    const double update_length = 0) noexcept;
  
  /** 
   Updates the filter coefficients. You can set how long it takes to 
   update the coefficients (using linear interpolation between old and new
   coefficients). If an update is requested while another is already in
   progress, the new interpolation will pick up from where the old one left
   off to avoid audible artifacts.
   @param[in] update_length How many calls to Filter it takes to update the
   coefficients. A value of 0 means that the update is instantaneous. A call
   to Filter(const T input) counts as one, just like
   Filter(const Vector<T>& input).
   */
  void SetNumeratorCoeffs(
    const Vector<T>& numerator_coeffs,
    const double update_length = 0) noexcept;

  /**
   Updates the filter coefficients. You can set how long it takes to
   update the coefficients (using linear interpolation between old and new
   coefficients). If an update is requested while another is already in
   progress, the new interpolation will pick up from where the old one left
   off to avoid audible artifacts.
   @param[in] update_length How many calls to Filter it takes to update the
   coefficients. A value of 0 means that the update is instantaneous. A call
   to Filter(const T input) counts as one, just like
   Filter(const Vector<T>& input).
   */
  void SetDenominatorCoeffs(
    const Vector<T>& denominator_coeffs,
    const double update_length = 0) noexcept;

  const Vector<T>& GetNumeratorCoeffs() const noexcept;
  
  const Vector<T>& GetDenominatorCoeffs() const noexcept;

  /** Resets the state of the filter */
  void ResetState() noexcept;

  void FilterAdd(
    const Vector<T>& input_to_filter,
    const Vector<T>& input_to_add,
    Vector<T>& output) noexcept;
  
  template<typename InputIterator, typename OutputIterator>
  void Filter(
    InputIterator input_begin,
    InputIterator input_end,
    OutputIterator output_begin) noexcept;
  
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
  Vector<T> temp_conv_input_;
  Vector<T> temp_conv_output_;
  Vector<T> temp_add_data_;

  /**
   @param[in] input_begin
   @param[in] input_end
   @param[out] extended_input
  */
  template<typename InputIterator>
  void GetExtendedInput(
    InputIterator input_begin,
    InputIterator input_end,
    Vector<T>& extended_input) noexcept;

  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients() noexcept;
  
  /**
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input.
   As a second example, if B=[0,1], A=[1], you will have
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  T FilterSample(
    const T input_sample) noexcept;
  
  template<typename InputIterator, typename OutputIterator>
  void FilterSerial(
    InputIterator input_begin,
    InputIterator input_end,
    OutputIterator output_begin) noexcept;
};


template<typename T>
DigitalFilter<T> GainFilter(T gain)
{
  return DigitalFilter<T>(UnaryVector<T>(gain));
}

} // namespace mcl

#include "digitalfilter.ipp"
