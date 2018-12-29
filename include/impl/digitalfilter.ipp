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
  if (denominator_coeffs_.size() == 1)
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
  else
  {
    // Transposed direct form II (this appears to be slightly slower)
//    T output = B_[0]*input + state_[0];
//    for (size_t i = 0; i <= size-2; ++i)
//    {
//      state_[i] = state_[i+1] + input*B_[i+1] - output*A_[i+1];
//    }

    // Direct form II
    T v = input_sample; // The temporary value in the recursive branch.
    T output(static_cast<T>(0.0));
    // The index i in both loops refers to the branch in the classic plot of a
    // direct form II, with the highest branch (the one multiplied by b(0) only)
    // being i=0.
    for (size_t i = 1; i < length_; ++i)
    {
      v += state_[i - 1] * (-denominator_coeffs_[i]);
      output += state_[i - 1] * numerator_coeffs_[i];
    }
    for (size_t i = (length_ - 1); i >= 1; --i)
    {
      state_[i] = state_[i - 1];
    }
    state_[0] = v;
    output += v * numerator_coeffs_[0];
    return output;
  }
}


template<typename T>
void DigitalFilter<T>::SetNumeratorCoeffs(
  const Vector<T>& numerator_coeffs,
  const double update_length) noexcept
{
  ASSERT(denominator_coeffs_.size() == 1 ||
    numerator_coeffs.size() == denominator_coeffs_.size());
  
  if (numerator_coeffs.size() != numerator_coeffs_.size())
  {
    // If the impulse response changes length, then reset everything.
    length_ = numerator_coeffs.size();
    state_ = mcl::Vector(length_, 0.0);
    counter_ = length_-1;
    numerator_coeffs_ = numerator_coeffs;
    numerator_smoothers_ = Vector<RampSmoother<T>>
    (
      numerator_coeffs.size(),
      RampSmoother(T(0.0)));
    // TODO: for N current length and M new length, make smooth transition for
    // the first N coefficients
  }
  
  if (update_length < 1.0)
  {
    std::copy
    (
      numerator_coeffs.begin(),
      numerator_coeffs.end(),
      numerator_coeffs_.begin());
  }
  else
  {
    for (size_t i=0; i<length_; ++i)
    {
      numerator_smoothers_[i].SetTargetValue(numerator_coeffs[i], update_length);
    }
  }
}


template<typename T>
void DigitalFilter<T>::SetDenominatorCoeffs(
  const Vector<T>& denominator_coeffs,
  const double /* update_length */ ) noexcept
{
  // TODO: implement slow update
  ASSERT(denominator_coeffs.size() == denominator_coeffs_.size());
  ASSERT(IsApproximatelyEqual(denominator_coeffs[0], 1.0, std::numeric_limits<T>::epsilon()));
  
  std::copy
  (
    denominator_coeffs.begin(),
    denominator_coeffs.end(),
    denominator_coeffs_.begin());
}


template<typename T>
void DigitalFilter<T>::CopyParametersFrom(
  const DigitalFilter<T>& other_filter,
  const double update_length) noexcept
{
  this->SetNumeratorCoeffs(other_filter.numerator_coeffs_);
  this->SetDenominatorCoeffs(other_filter.denominator_coeffs_);
}


template<typename T>
void DigitalFilter<T>::ResetState() noexcept
{
  SetToZero(state_);
}


template<typename T>
Vector<T> DigitalFilter<T>::GetNumeratorCoeffs() noexcept
{
  return numerator_coeffs_;
}
  
  
template<typename T>
Vector<T> DigitalFilter<T>::GetDenominatorCoeffs() noexcept
{
  return denominator_coeffs_;
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
  
  if (denominator_coeffs_.size() > 1 || input.size() < length_)
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
Vector<T> DigitalFilter<T>::Filter(
  const Vector<T>& input) noexcept
{
  Vector<T> output(input.size());
  Filter(input, output);
  return std::move(output);
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
}


/**
 Get wall filters of type wall_type and for FS given by sampling_frequency
 */
template<typename T>
DigitalFilter<T> WallFilter(
  WallType wall_type,
  T sampling_frequency)
{
  // TODO: implement for frequencies other than 44100
  if (! IsEqual(sampling_frequency, 44100))
  {
    Logger::GetInstance().LogError
    (
      "Attempting to use a wall filter "
      "designed for 44100 Hz sampling frequency with a sampling frequency "
      "of %f Hz. The filter response will be inaccurate.",
      sampling_frequency);
  }

  Vector<T> B;
  Vector<T> A;

  switch (wall_type)
  {
  case kRigid:
  {
    B[0] = 1.0;
    A[0] = 1.0;
    break;
  }
  case kCarpetPile:
  {
    // B_carpet_pile=[0.562666833756030  -1.032627191365576   0.469961155406544];
    // A_carpet_pile=[1.000000000000000  -1.896102349247713   0.896352947528892];
    Vector<T> B_(3);
    Vector<T> A_(3);
    B_[0] = 0.562666833756030;
    B_[1] = -1.03262719136557;
    B_[2] = 0.469961155406544;
    A_[0] = 1.000000000000000;
    A_[1] = -1.896102349247713;
    A_[2] = 0.896352947528892;
    B = B_;
    A = A_;
    break;
  }
  case kCarpetCotton:
  {
    // B_carpet_cotton = [0.687580695329600  -1.920746652319969   1.789915765926473  -0.556749690855965];
    // A_carpet_cotton = [1.000000000000000  -2.761840732459190   2.536820778736938  -0.774942833868750];
    Vector<T> B_(4);
    Vector<T> A_(4);
    B_[0] = 0.687580695329600;
    B_[1] = -1.920746652319969;
    B_[2] = 1.789915765926473;
    B_[3] = -0.556749690855965;
    A_[0] = 1.000000000000000;
    A_[1] = -2.761840732459190;
    A_[2] = 2.536820778736938;
    A_[3] = -0.774942833868750;
    B = B_;
    A = A_;
    break;
  }
  case kWallBricks:
  {
    // B_wall=[0.978495798553620  -1.817487798457697   0.839209660516074];
    // A_wall=[1.000000000000000  -1.858806492488240   0.859035906864860];
    Vector<T> B_(3);
    Vector<T> A_(3);
    B_[0] = 0.978495798553620;
    B_[1] = -1.817487798457697;
    B_[2] = 0.839209660516074;
    A_[0] = 1.000000000000000;
    A_[1] = -1.858806492488240;
    A_[2] = 0.859035906864860;
    B = B_;
    A = A_;
    break;
  }
  case kCeilingTile:
  {
    // B_ceiling=[0.168413736374283  -0.243270224986791   0.074863520490536];
    // A_ceiling=[1.000000000000000  -1.845049094190385   0.845565720138466];
    Vector<T> B_(3);
    Vector<T> A_(3);
    B_[0] = 0.168413736374283;
    B_[1] = -0.243270224986791;
    B_[2] = 0.074863520490536;
    A_[0] = 1.000000000000000;
    A_[1] = -1.845049094190385;
    A_[2] = 0.845565720138466;
    B = B_;
    A = A_;
    break;
  }
  default:
  {
    Logger::GetInstance().LogError(
      "Unrecognised type of wall filter. Reverting to a completely absorptive one.");
    Vector<T> B_(1);
    Vector<T> A_(1);
    B_[0] = 0.0;
    A_[0] = 1.0;
    B = B_;
    A = A_;
    break;
  }
  }
  return DigitalFilter<T>(B, A);
}


/** Returns a pinkifier filter */
template<typename T>
DigitalFilter<T> PinkifierFilter()
{
  Vector<T> poles(5);
  poles[0] = 0.9986823;
  poles[1] = 0.9914651;
  poles[2] = 0.9580812;
  poles[3] = 0.8090598;
  poles[4] = 0.2896591;
  Vector<T> zeros(5);
  zeros[0] = 0.9963594;
  zeros[1] = 0.9808756;
  zeros[2] = 0.9097290;
  zeros[3] = 0.6128445;
  zeros[4] = -0.0324723;

  Vector<Complex<T>> num = Poly(zeros);
  Vector<Complex<T>> den = Poly(poles);

  return DigitalFilter<T>(RealPart<T>(num), RealPart<T>(den));
}

} // namespace mcl
