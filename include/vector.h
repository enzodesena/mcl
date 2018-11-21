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

constexpr size_t kDynamicLength = std::numeric_limits<size_t>::max();
constexpr size_t kReferenced = std::numeric_limits<size_t>::max()-1;

template<typename T = double, size_t vector_length = kDynamicLength>
class Vector : private std::array<T, vector_length>
{
public:
  using std::array<T,vector_length>::operator[];
  using Iterator = typename std::array<T, vector_length>::iterator;
  using std::array<T,vector_length>::begin;
  using std::array<T,vector_length>::end;
  using std::array<T,vector_length>::data;
  
  Vector() noexcept
  {
  }
  
  Vector(
    size_t length,
    T value = T()) noexcept
  {
    ASSERT(length == vector_length);
    for (auto iter = begin(); iter != end(); ++iter)
    {
      *iter = value;
    }
  }
  
  inline size_t length() const noexcept
  {
    return static_cast<int>(vector_length);
  }
  
  inline T& operator[](const size_t index) noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
  
  inline const T& operator[](const size_t index) const noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
};


template<typename T>
class Vector<T, kDynamicLength> : private std::vector<T>
{
public:
  using Iterator = typename std::vector<T>::iterator;
  using ConstIterator = typename std::vector<T>::const_iterator;
  using std::vector<T>::begin;
  using std::vector<T>::end;
//  using std::vector<T>::data;
  
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

  inline size_t length() const noexcept
  {
    return static_cast<int>(this->size());
  }
  
  inline T& operator[](const size_t index) noexcept
  {
    ASSERT(index>=0 && index < length());
    return std::vector<T>::data()[index];
  }
  
  inline const T& operator[](const size_t index) const noexcept
  {
    ASSERT(index>=0 && index < length());
    return std::vector<T>::data()[index];
  }
  
  inline void PushBack(const T& element) noexcept
  {
    std::vector<T>::push_back(element);
  }
};


template<typename T>
class Vector<T,kReferenced> {
public:
  Vector(Vector<T,kReferenced>& other_vector) : other_vector_(other_vector)
  {
  }
  
  inline size_t length()
  {
    return other_vector_.length();
  }
  
private:
  Vector<T,kReferenced>& other_vector_;
};



} /**< namespace mcl */
