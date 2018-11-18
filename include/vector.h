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

constexpr int kDynamicLength = -1;
constexpr int kReferenced = -2;

template<typename T = double, int vector_length = kDynamicLength>
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
    int length,
    T value = T()) noexcept
  {
    ASSERT(length == vector_length);
  }
  
  inline int length() const noexcept
  {
    return static_cast<int>(vector_length);
  }
  
  inline T& operator[](const int index) noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
  
  inline const T& operator[](const int index) const noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
};


template<typename T>
class Vector<T, kDynamicLength> : private std::vector<T>
{
public:
  Vector() noexcept : std::vector<T>()
  {
  }
  
  Vector(
    int initial_length,
    T value = T()) noexcept
    : std::vector<T>(initial_length, value)
  {
  }
  
  
  using std::vector<T>::operator[];
  using Iterator = typename std::vector<T>::iterator;
  using std::vector<T>::begin;
  using std::vector<T>::end;
  using std::vector<T>::data;

  inline int length() const noexcept
  {
    return static_cast<int>(this->size());
  }
  
  inline T& operator[](const int index) noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
  
  inline const T& operator[](const int index) const noexcept
  {
    ASSERT(index>=0 && index < length());
    return data()[index];
  }
  
  inline void PushBack(const T& element) noexcept
  {
    std::vector<T>::push_back(element);
  }
  
//  inline void SetElements(
//    const Iterator&
//    const Iterator& new_elements_begin,
//    const Iterator& new_elements_end) noexcept
//  {
//    auto iter =
//    for (auto iter = begin()+index;
//         iter < end() && new_elements_begin < new_elements_end;
//         iter++)
//    {
//      *iter = *(new_elements_begin++);
//    }
//  }
};


//template<typename T>
//class Vector<T,kReferenced> {
//public:
//  Vector(Vector<T,kReferenced>)
//  inline int length() {
//    return this->size();
//  }
//  
//  inline T At(int index) {
//    ASSERT(index>=0 & index<this->size());
//    return this[index];
//  }
//};



} /**< namespace mcl */
