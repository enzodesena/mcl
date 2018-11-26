/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include <array>
#include <vector>
#include "mcltypes.h"

namespace mcl
{

/** Test function for the functions in this file */
bool VectorOpTest();



template<typename T>
class Vector : private std::vector<T>
{
public:
  using Iterator = typename std::vector<T>::iterator;
  using ConstIterator = typename std::vector<T>::const_iterator;
  
  Vector() noexcept : std::vector<T>()
  {
  }
  
  Vector(
    size_t initial_length,
    T value = T()) noexcept
    : std::vector<T>(initial_length, value)
  {
  }

  Vector(
    ConstIterator iter_begin,
    ConstIterator iter_end) noexcept
    : std::vector<T>(iter_begin, iter_end)
  {
  }
  
  Iterator begin() noexcept
  {
    return std::vector<T>::begin();
  }

  Iterator end() noexcept
  {
    return std::vector<T>::end();
  }

  ConstIterator begin() const noexcept
  {
    return std::vector<T>::begin();
  }

  ConstIterator end() const noexcept
  {
    return std::vector<T>::end();
  }

  Vector(
    std::initializer_list<T> list)
    : std::vector<T>(list)
  {
  }
  
  inline size_t length() const noexcept
  {
    return this->size();
  }
  
  inline T& operator[](
    const size_t index) noexcept
  {
    ASSERT(index>=0 && index < length());
    return std::vector<T>::data()[index];
  }
  
  inline const T& operator[](
    const size_t index) const noexcept
  {
    ASSERT(index>=0 && index < length());
    return std::vector<T>::data()[index];
  }
  
  inline void PushBack(const T& element) noexcept
  {
    std::vector<T>::push_back(element);
  }
};

//
//template<typename T>
//class VectorRef<T> {
//private:
//  Vector<T>& other_vector_;
//  size_t start_;
//  size_t length_;
//public:
//  using Iterator = typename std::vector<T>::iterator;
//  using ConstIterator = typename std::vector<T>::const_iterator;
//
//  inline size_t length() const noexcept
//  {
//    // The default for `num_elements_` is std::numeric_limits<size_t>::max(),
//    // which means this reference vector has the same length as the
//    // `other_vector`.
//    // If, on the other hand, num_elements_ is a smaller number, it means
//    // that the length of this reference vector is smaller than the
//    // length of the `other_vector`.
//    // If the length of the other vector is shortened (not possible in the
//    // current implementation), an assert will happen when using the []
//    // operator with an index ending up outside the vector.
//    return std::min(other_vector_.length(), length_);
//  }
//
//  inline Iterator begin() noexcept
//  {
//    return other_vector_.begin() + start_;
//  }
//
//  inline ConstIterator begin() const noexcept
//  {
//    return other_vector_.begin() + start_;
//  }
//
//  inline Iterator end() noexcept
//  {
//    return other_vector_.end() + start_ - (other_vector_.length() - length());
//  }
//
//  inline ConstIterator end() const noexcept
//  {
//    return other_vector_.end() + start_ - (other_vector_.length() - length());
//  }
//
//  Vector(
//    Vector<T>& other_vector,
//    size_t start = 0,
//    size_t length = std::numeric_limits<size_t>::max()) noexcept
//    : other_vector_(other_vector)
//    , start_(start)
//    , length_(length)
//  {
//    if (length_ < std::numeric_limits<size_t>::max()) // TODO: remove this if not debugging
//    {
//      ASSERT(length_ <= other_vector.length());
//    }
//    ASSERT(start >= 0 && start < other_vector.length());
//    ASSERT((start+this->length()) <= other_vector.length());
//  }
//
//  inline T& operator[](
//    const size_t index) noexcept
//  {
//    ASSERT(index>=0 && index < length());
//    const size_t other_index(index + start_);
//    ASSERT(other_index>=0 && other_index < other_vector_.length());
//    return other_vector_[other_index];
//  }
//
//  inline const T& operator[](
//    const size_t index) const noexcept
//  {
//    ASSERT(index>=0 && index < length());
//    const size_t other_index(index + start_);
//    ASSERT(other_index>=0 && other_index < other_vector_.length());
//    return other_vector_[other_index];
//  }
//};



} /**< namespace mcl */
