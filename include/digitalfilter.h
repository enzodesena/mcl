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
  
/** Digital filter abstract class */
template<typename T>
class DigitalFilter
{
public:
  /** Returns the output of the filter for an input equal to `input` . */
  virtual T Filter(
    const T input) noexcept
  {
    return Filter(UnaryVector<T>(input))[0];
  }
  
  template<size_t length>
  Vector<T,length> Filter(
    const Vector<T,length>& input) noexcept
  {
    Vector<T,length> output(input.length());
    Filter(input, output);
    return std::move(output);
  }
  
  virtual void FilterSerial(
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
  
  /** Returns the output of the filter for an input signal equal to `input`. */
  virtual void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept = 0;
  
  /** Resets the state of the filter */
  virtual void Reset() = 0;
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
class GainFilter : public DigitalFilter<T> {
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
    ASSERT(input.length() == output.length());
    Multiply(input, input.length(), gain_, output);
  }
  
  void Reset()
  {
  }
};

template<typename T>
class IdenticalFilter : public DigitalFilter<T> {
public:
  void Filter(
    const Vector<T>& input,
    Vector<T>& output) noexcept
  {
    output = input;
  }
  
  void Reset()
  {
  }
};
  
  
} // namespace mcl
