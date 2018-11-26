/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vector.h"
#include "vectorop.h"
#include "mcltypes.h"

namespace mcl
{


//class DigitalFilter
//{
//public:
//  template<typename DigitalFilterType>
//  DigitalFilter(DigitalFilterType x)
//    : self_(std::make_unique<DigitalFilterModel<DigitalFilterType>>(std::move(x)))
//  {
//  }
//  
//  DigitalFilter(
//    const DigitalFilter& x)
//    : self_(x.self_->copy_())
//  {
//  }
//  
//  DigitalFilter(
//    DigitalFilter&& x) noexcept = default;
//  
//  DigitalFilter& operator=(
//    const DigitalFilter& x) noexcept
//  {
//    DigitalFilter tmp(x);
//    *this = std::move(tmp); // Using move assignment operator
//    return *this;
//  }
//  
//  /** Move assignment operator */
//  DigitalFilter& operator=(
//    DigitalFilter&& x) noexcept = default;
//  
//  
//  template<typename T, size_t input_length, size_t output_length>
//  void Filter(
//    const Vector<T,input_length>& input,
//    Vector<T,output_length>& output) noexcept
//  {
//    self_->Filter_(input, output);
//  }
//  
//  void Reset() noexcept
//  {
//    self_->Reset_();
//  }
//private:
//  struct DigitalFilterConcept {
//    virtual ~DigitalFilterConcept() = default;
//    template<typename T, size_t input_length, size_t output_length>
//    virtual void Filter_(
//      const Vector<T,input_length>& input,
//      Vector<T,output_length>& output) = 0;
//    virtual void Reset_() = 0;
//    virtual std::unique_ptr<DigitalFilterConcept> copy_() = 0;
//  };
//  
//  template<typename DigitalFilterType>
//  struct DigitalFilterModel final : DigitalFilterConcept
//  {
//    DigitalFilterModel(DigitalFilterType x)
//      : data_(std::move(x))
//    {
//    }
//    
//    std::unique_ptr<DigitalFilterConcept> copy_() override
//    {
//      return std::make_unique<DigitalFilterModel>(*this);
//    }
//    
//    template<typename T, size_t input_length, size_t output_length>
//    void Filter_(
//      const Vector<T,input_length>& input,
//      Vector<T,output_length>& output) noexcept override
//    {
//      data_.Filter(out, position);
//    }
//    
//    void Reset_() override noexcept
//    {
//      data_.Reset();
//    }
//    
//    DigitalFilterType data_;
//  };
//  
//  std::unique_ptr<DigitalFilterConcept> self_; // Concept is drawable object
//};

  
  
/** Filter bank abstract class */
template<typename T>
class FilterBank {
public:
  /** Returns the output of the filter bank for an input equal to `input`. Hello world! */
  virtual Vector<T> Filter(
    const T input) = 0;
  
  /** Returns the output of the filter bank for a given input. */
  virtual Vector<Vector<T>> Filter(
    const Vector<T>& input) = 0;
  
  /** Resets the state of the filter */
  virtual void Reset() = 0;
  
  virtual Int num_filters() = 0;
};

//template<typename T>
//class GainFilter : public DigitalFilter<T> {
//private:
//  T gain_;
//public:
//  GainFilter(
//    const T gain)
//    : gain_(gain)
//  {
//  }
//  
//  void Filter(
//    const Vector<T>& input,
//    Vector<T>& output) noexcept
//  {
//    ASSERT(input.length() == output.length());
//    Multiply(input, input.length(), gain_, output);
//  }
//  
//  void Reset()
//  {
//  }
//};
//
//template<typename T>
//class IdenticalFilter : public DigitalFilter<T> {
//public:
//  void Filter(
//    const Vector<T>& input,
//    Vector<T>& output) noexcept
//  {
//    output = input;
//  }
//  
//  void Reset()
//  {
//  }
//};
//  
  
} // namespace mcl
