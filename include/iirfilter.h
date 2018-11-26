/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "digitalfilter.h"
#include "vectorop.h"
#include "constants.h"

namespace mcl {

  
/** IIR Filter */
template<typename T>
class IirFilter : public DigitalFilterInterface<T> {
public:
  /** Constructs a default filter, i.e. identical filter*/
  IirFilter()
    : B_(mcl::UnaryVector<T>(1.0))
    , A_(mcl::UnaryVector<T>(1.0))
  {
    Int size = B_.length();
    state_ = Vector<T>(size);
    for (Int i=0; i<size; ++i) { state_[i] = 0.0; }
  }
  
  /** 
   Constructs an IIR filter (II-type direct implementation). B and A are numerator
   and denominator of the filter, respectively.
   */
  IirFilter(
    const Vector<T>& B,
    const Vector<T>& A)
    : B_(B)
    , A_(A)
  {
    // TODO: implement also for B.size != A.size
    ASSERT(B.length() == A.length());
    ASSERT(B.length() >= 1);
    
    A0_ = A[0];
    if (! IsApproximatelyEqual(A[0], 1.0, std::numeric_limits<T>::epsilon())) {
      B_ = Multiply(B, (T) 1.0 / A[0]);
      A_ = Multiply(A, (T) 1.0 / A[0]);
    }
    
    
    Int size = B.length();
    state_ = Vector<T>(size);
    for (Int i=0; i<size; ++i) { state_[i] = 0.0; }
  }
  
  /** 
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input. 
   As a second example, if B=[0,1], A=[1, 0], you will have
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  T FilterSample(
    const T input) noexcept
  {
    // Speed up return for simple gain filters
    if (B_.length() == 1) { return input*B_[0]; }
    size_t size = B_.length();
    ASSERT(size >= 1);
    T v = input; // The temporary value in the recursive branch.
    T output(0.0);
    // The index i in both loops refers to the branch in the classic plot of a
    // direct form II, with the highest branch (the one multiplied by b(0) only)
    // being i=0.
    for (size_t i=1; i<size; ++i)
    {
      v += state_[i-1]*(-A_[i]);
      output += state_[i-1]*B_[i];
    }
    for (size_t i=(size-1); i>=1; --i)
    {
      state_[i] = state_[i-1];
    }
    state_[0] = v;
    output += v*B_[0];
    return output;
  }
  
  T Filter(
    const T input) noexcept
  {
    return FilterSample(input);
  }
  
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    if (B_.length() == 1)
    {
      Multiply(input, B_[0], output);
    }
    else
    {
      ForEach<T,T>(
        input,
        [this] (T element) { return this->FilterSample(element); },
        output);
    }
  }
  
  Vector<T> Filter(
    const Vector<T>& input) noexcept
  {
    Vector<T> output(input.length());
    Filter(input, output);
    return std::move(output);
  }
  
  /** Returns the order of the filter. */
  Int order() const noexcept
  {
    return Max(B_.length(), A_.length())-1;
  }
  
  /** 
   Updates the filter coefficients. May cause articafts if the coefficients are
   updated too rapidly.
   */
  void SetCoefficients(
    const Vector<T>& B,
    const Vector<T>& A) noexcept
  {
    // TODO: implement case where length changes.
    ASSERT(B_.length() == B.length());
    ASSERT(A_.length() == A.length());
    B_ = B;
    A_ = A;
    A0_ = A[0];
  }
  
  /** Sets the coefficients as identical to those of another filter. */
  void SetCoefficients(
    const IirFilter& other_filter) noexcept
  {
    const Int filter_order = order();
    assert(filter_order == other_filter.order());
    
    for (Int i=0; i<=filter_order; ++i) {
      SetNumeratorCoefficient(i, other_filter.GetNumeratorCoefficient(i));
      SetDenominatorCoefficient(i, other_filter.GetDenominatorCoefficient(i));
    }
  }
  
  T GetNumeratorCoefficient(
    const size_t coeff_id) const noexcept
  {
    ASSERT(coeff_id>=0 && coeff_id<=order());
    return B_[coeff_id]*A0_;
  }
  
  T GetDenominatorCoefficient(
    const size_t coeff_id) const noexcept
  {
    return A_.at(coeff_id);
  }
  
  void SetNumeratorCoefficient(
    const size_t coeff_id,
    const T value) noexcept
  {
    ASSERT(coeff_id >= 0 && coeff_id<=order());
    B_[coeff_id] = value/A0_;
  }
  
  
  inline void SetDenominatorCoefficient(
    const size_t coeff_id,
    const T value) noexcept
  {
    ASSERT(coeff_id >= 0 &&coeff_id<(Int)A_.length());
    A_[coeff_id] = value;
    if (coeff_id == 0) { A0_ = value; }
  }
  
  /** Returns the forward coefficients */
  Vector<T> B() const
  {
    // Return the non-normalised version
    return Multiply(B_, A0_);
  }
  
  /** Returns the backward coefficients */
  Vector<T> A() const
  {
    // Return the non-normalised version
    return Multiply(A_, A0_);
  }
  
  void Reset() noexcept
  {
    SetToZero(state_);
  }
  
private:
  Vector<T> B_;
  Vector<T> A_;
  
  // By storing A[0] before normalisation we can output B() and A() before
  // normalisation while keeping the internal representation normalised
  T A0_;
  
  Vector<T> state_;
};
  
/** Filter bank abstract class */
template<typename T>
class IirFilterBank : public FilterBank<T> {
private:
  Vector<IirFilter<T>> filters_;
  
public:
  IirFilterBank(
    const Vector<IirFilter<T>>& filters) noexcept
    : filters_(filters)
  {
  }
  
  virtual Int num_filters() noexcept
  {
    return filters_.length();
  }
  
  /** Returns the output of the filter bank for an input equal to `input`. */
  virtual Vector<T> Filter(
    const T input) noexcept
  {
    const size_t N = filters_.length();
    Vector<T> outputs(N);
    for (size_t i=0; i<N; ++i)
    {
      outputs[i] = filters_[i].Filter(input);
    }
    return outputs;
  }
  
  /** Returns the output of the filter bank for a given input. */
  virtual Vector<Vector<T> >
  Filter(const Vector<T>& input)
  {
    const size_t N = filters_.length();
    Vector<Vector<T> > outputs(N);
    for (size_t i=0; i<N; ++i)
    {
      outputs[i] = filters_[i].Filter(input);
    }
    return outputs;
  }
  
  /** Resets the state of the filter */
  void Reset()
  {
    for (size_t i=0; i<filters_.length(); ++i)
    {
      filters_[i].Reset();
    }
  }
};
  
//
//  /** Implements a first-order IIR low-pass filter with a given decay constant. */
//  class RampSmoothingFilter : public DigitalFilter {
//  public:
//    
//    /**
//     @param[in] ramp_samples number of samples after which the value is
//     to 1/e away from target value. */
//    RampSmoothingFilter(const T ramp_samples) noexcept {
//      ASSERT_WITH_MESSAGE(std::isgreaterequal(ramp_samples, 0),
//                          "Decay constant cannot be negative ");
//      
//      
//    }
//    
//    virtual T Filter(const T input) noexcept {
//      return filter_.Filter(input);
//    }
//    
//    using DigitalFilter::Filter;
//    
//    virtual void Reset() noexcept { filter_.Reset(); }
//    
//    
//  private:
//  };
//  
///** Implements a first-order IIR low-pass filter with a given decay constant. */
//class LowPassSmoothingFilter : public DigitalFilter {
//public:
//  
//  /**
//   @param[in] ramp_samples number of samples after which the value is
//   to 1/e away from target value. */
//  LowPassSmoothingFilter(const T ramp_samples) noexcept {
//    ASSERT_WITH_MESSAGE(std::isgreaterequal(ramp_samples, 0),
//                        "Decay constant cannot be negative ");
//    
//    T a1 = exp(-1.0/ramp_samples);
//    T b0 = 1.0 - a1;
//    filter_ = IirFilter(mcl::BinaryVector(b0, 0.0),
//                        mcl::BinaryVector(1.0, -a1));
//  }
//  
//  virtual T Filter(const T input) noexcept {
//    return filter_.Filter(input);
//  }
//  
//  using DigitalFilter::Filter;
//  
//  virtual void Reset() noexcept { filter_.Reset(); }
//  
//  
//private:
//  IirFilter filter_;
//};
  

/** Constructs a filter for which output==gain*input always. */
template<typename T>
inline IirFilter<T> GainFilter(
  const T gain) noexcept
{
  Vector<T> B(1);
  Vector<T> A(1);
  B[0] = 1.0*gain;
  A[0] = 1.0;
  return IirFilter<T>(B,A);
}


/** Constructs a filter for which output==input always. */
template<typename T>
inline IirFilter<T> IdenticalFilter() noexcept
{
  return GainFilter(1.0);
}

/**
 Get wall filters of type wall_type and for FS given by sampling_frequency
 */
template<typename T>
inline IirFilter<T> WallFilter(
  WallType wall_type, T sampling_frequency)
{
  // TODO: implement for frequencies other than 44100
  if (! IsEqual(sampling_frequency, 44100))
  {
    mcl::Logger::GetInstance().LogError("Attempting to use a wall filter "
        "designed for 44100 Hz sampling frequency with a sampling frequency "
        "of %f Hz. The filter response will be inaccurate.", sampling_frequency);
  }
  
  Vector<T> B;
  Vector<T> A;
  
  switch (wall_type)
  {
    case kRigid:
    {
      B.PushBack(1.0);
      A.PushBack(1.0);
      break;
    }
    case kCarpetPile:
    {
      // B_carpet_pile=[0.562666833756030  -1.032627191365576   0.469961155406544];
      // A_carpet_pile=[1.000000000000000  -1.896102349247713   0.896352947528892];
      B.PushBack(0.562666833756030);
      B.PushBack(-1.032627191365576);
      B.PushBack(0.469961155406544);
      A.PushBack(1.000000000000000);
      A.PushBack(-1.896102349247713);
      A.PushBack(0.896352947528892);
      break;
    }
    case kCarpetCotton:
    {
      // B_carpet_cotton = [0.687580695329600  -1.920746652319969   1.789915765926473  -0.556749690855965];
      // A_carpet_cotton = [1.000000000000000  -2.761840732459190   2.536820778736938  -0.774942833868750];
      B.PushBack(0.687580695329600);
      B.PushBack(-1.920746652319969);
      B.PushBack(1.789915765926473);
      B.PushBack(-0.556749690855965);
      A.PushBack(1.000000000000000);
      A.PushBack(-2.761840732459190);
      A.PushBack(2.536820778736938);
      A.PushBack(-0.774942833868750);
      break;
    }
    case kWallBricks:
    {
      // B_wall=[0.978495798553620  -1.817487798457697   0.839209660516074];
      // A_wall=[1.000000000000000  -1.858806492488240   0.859035906864860];
      B.PushBack(0.978495798553620);
      B.PushBack(-1.817487798457697);
      B.PushBack(0.839209660516074);
      A.PushBack(1.000000000000000);
      A.PushBack(-1.858806492488240);
      A.PushBack(0.859035906864860);
      break;
    }
    case kCeilingTile:
    {
      // B_ceiling=[0.168413736374283  -0.243270224986791   0.074863520490536];
      // A_ceiling=[1.000000000000000  -1.845049094190385   0.845565720138466];
      B.PushBack(0.168413736374283);
      B.PushBack(-0.243270224986791);
      B.PushBack(0.074863520490536);
      A.PushBack(1.000000000000000);
      A.PushBack(-1.845049094190385);
      A.PushBack(0.845565720138466);
      break;
    }
    default:
    {
      mcl::Logger::GetInstance().LogError("Unrecognised type of wall filter. Reverting to a completely absorptive one.");
      B.PushBack(0.0);
      A.PushBack(1.0);
      break;
    }
  }
  return IirFilter<T>(B,A);
}

/** Returns a pinkifier filter */
template<typename T>
inline IirFilter<T> PinkifierFilter()
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
  
  return IirFilter<T>(RealPart<T>(num),RealPart<T>(den));
}

bool IirFilterTest();
  
} // namespace mcl
