/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "vector.hpp"

namespace mcl
{

template<typename T>
class DigitalFilterInterface
{
  public:
  virtual void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept = 0;

  virtual void Reset() noexcept = 0;
};


template<typename T>
class DigitalFilter
{
public:
  template<typename DigitalFilterT>
  DigitalFilter(DigitalFilterT x)
    : self_(std::make_unique<DigitalFilterModel<DigitalFilterT>>(std::move(x)))
  {
  }

  DigitalFilter(
    const DigitalFilter& x)
    : self_(x.self_->copy_())
  {
  }

  DigitalFilter(
    DigitalFilter&& x) noexcept = default;

  DigitalFilter& operator=(
    const DigitalFilter& x) noexcept
  {
    DigitalFilter tmp(x);
    *this = std::move(tmp); // Using move assignment operator
    return *this;
  }

  /** Move assignment operator */
  DigitalFilter& operator=(
    DigitalFilter&& x) noexcept = default;


  void Filter(
    const Vector<double>& input,
    Vector<double>& output) noexcept
  {
    self_->Filter_(input, output);
  }
  
  void Filter(
    const Vector<float>& input,
    Vector<float>& output) noexcept
  {
    self_->Filter_(input, output);
  }

  void Reset() noexcept
  {
    self_->Reset_();
  }
private:
  struct DigitalFilterConcept {
    virtual ~DigitalFilterConcept() = default;
    virtual void Filter_(
      const Vector<T>& input,
      Vector<T>& output) = 0;
    virtual void Reset_() = 0;
    virtual std::unique_ptr<DigitalFilterConcept> copy_() = 0;
  };

  template<typename DigitalFilterT>
  struct DigitalFilterModel final : DigitalFilterConcept
  {
    DigitalFilterT data_;
    
    DigitalFilterModel(DigitalFilterT x)
      : data_(std::move(x))
    {
    }

    std::unique_ptr<DigitalFilterConcept> copy_() override
    {
      return std::make_unique<DigitalFilterModel>(*this);
    }

    void Filter_(
      const Vector<T>& input,
      Vector<T>& output) noexcept override
    {
      data_.Filter(input, output);
    }

    void Reset_() noexcept override
    {
      data_.Reset();
    }

  };

  std::unique_ptr<DigitalFilterConcept> self_; // Concept is filterable object
};

  
  
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

template<typename T>
class GainFilter {
private:
  T gain_;
public:
  GainFilter(
    const T gain)
    : gain_(gain)
  {
  }
  
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    ASSERT(input.size() == output.size());
    Multiply(input, gain_, output);
  }
  
  void Reset()
  {
  }
};

  
} // namespace mcl
