/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


namespace mcl
{

template<typename T>
DigitalFilter<T>::DigitalFilter() noexcept
  : DigitalFilter(UnaryVector<T>(1.0))
{
}


template<typename T>
DigitalFilter<T>::DigitalFilter(
  const Vector<T>& numerator_coeffs,
  const Vector<T>& denominator_coeffs,
  size_t max_expected_input_length) noexcept
  : state_(numerator_coeffs.size(), T(0.0))
  , numerator_coeffs_(numerator_coeffs)
  , denominator_coeffs_(denominator_coeffs)
  , counter_(numerator_coeffs.size() - 1)
  , length_(numerator_coeffs.size())
  , numerator_smoothers_(numerator_coeffs.size(), T(0.0))
  , denominator_smoothers_(denominator_coeffs.size(), T(0.0))
  , temp_conv_data_(Vector(max_expected_input_length + length_ - 1, T(0.0)))
  , temp_add_data_(Vector(max_expected_input_length, T(0.0)))
{
  ASSERT(denominator_coeffs.size() == 1 ||
    denominator_coeffs.size() == numerator_coeffs.size());
  ASSERT(IsApproximatelyEqual(denominator_coeffs[0], 1.0, std::numeric_limits<T>::epsilon()));
  for (size_t i=0; i<numerator_coeffs.size(); ++i)
  {
    numerator_smoothers_[i].SetTargetValue(numerator_coeffs[i], 0);
  }
  for (size_t i=0; i<denominator_coeffs.size(); ++i)
  {
    denominator_smoothers_[i].SetTargetValue(denominator_coeffs[i], 0);
  }
}

template<typename T>
T DigitalFilter<T>::Filter(
  const T input_sample) noexcept
{
  if (numerator_smoothers_[0].IsUpdating())
  {
    UpdateCoefficients();
  }
  return FilterSample(input_sample);
}

template<typename T>
T DigitalFilter<T>::FilterSample(
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


template<typename T>
void DigitalFilter<T>::SetNumeratorCoeffs(
  const Vector<T>& numerator_coeffs,
  const double update_length) noexcept
{
  ASSERT(numerator_coeffs.size() == length_);
  
  for (size_t i=0; i<length_; ++i)
  {
    numerator_smoothers_[i].SetTargetValue(numerator_coeffs[i], update_length);
  }
}


template<typename T>
void DigitalFilter<T>::SetDenominatorCoeffs(
  const Vector<T>& denominator_coeffs,
  const double update_length) noexcept
{
  ASSERT(denominator_coeffs.size() == denominator_coeffs_.size());
  ASSERT(IsApproximatelyEqual(denominator_coeffs[0], 1.0, std::numeric_limits<T>::epsilon()));
  
  for (size_t i=0; i<denominator_coeffs_.size(); ++i)
  {
    denominator_smoothers_[i].SetTargetValue(denominator_coeffs[i], update_length);
  }
}

template<typename T>
void DigitalFilter<T>::SetStateToZero() noexcept
{
  SetToZero(state_);
}


template<typename T>
Vector<T> DigitalFilter<T>::impulse_response() const noexcept
{
  return numerator_coeffs_;
}


template<typename T>
void DigitalFilter<T>::FilterSerial(
  const Vector<T>& input,
  Vector<T>& output) noexcept
{
  for (size_t i=0; i<input.size(); ++i) {
    output[i] = FilterSample(input[i]);
  }
}


template<typename T>
void DigitalFilter<T>::FilterAdd(
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


template<typename T>
void DigitalFilter<T>::Filter(
  const Vector<T>& input,
  Vector<T>& output) noexcept
{
  ASSERT(input.size() == output.size());
  
  if (numerator_smoothers_[0].IsUpdating())
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


template<typename T>
void DigitalFilter<T>::GetExtendedInput(
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

template<typename T>
void DigitalFilter<T>::UpdateCoefficients() noexcept
{
  for (size_t i=0; i<numerator_coeffs_.size(); ++i)
  {
    numerator_coeffs_[i] = numerator_smoothers_[i].GetNextValue();
  }
  for (size_t i=0; i<denominator_coeffs_.size(); ++i)
  {
    denominator_coeffs_[i] = denominator_smoothers_[i].GetNextValue();
  }
}
} // namespace mcl
