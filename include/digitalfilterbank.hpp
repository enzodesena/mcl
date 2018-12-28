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
/** Filter bank abstract class */
template<typename T>
class DigitalFilterBank
{
public:
  /** Returns the output of the filter bank for an input equal to `input`. Hello world! */
  virtual Vector<T> Filter(
    T input) = 0;

  /** Returns the output of the filter bank for a given input. */
  virtual Vector<Vector<T>> Filter(
    const Vector<T>& input) = 0;

  /** Resets the state of the filter */
  virtual void Reset() = 0;

  virtual Int num_filters() = 0;
};


} // namespace mcl
