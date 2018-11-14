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
  using std::array<T, vector_length>::operator[];
  using Iterator = typename std::array<T, vector_length>::iterator;
  using std::array<T, vector_length>::begin;
  using std::array<T, vector_length>::end;
  using std::array<T, vector_length>::data;
  
  Vector () noexcept {}
  
  Vector(int length) noexcept { ASSERT(length == vector_length); }
  
  inline int length()
  {
    return static_cast<int>(vector_length);
  }
  
  inline T& At(int index) noexcept
  {
    ASSERT(index > 0 & index < length());
    return std::array<T, vector_length>::at(index);
  }
  
  inline const T& At(int index) const noexcept
  {
    ASSERT(index > 0 & index < length());
    return std::array<T, vector_length>::at(index);
  }
  
};


template<typename T>
class Vector<T, kDynamicLength> : private std::vector<T>
{
public:
  using std::vector<T>::operator[];
  using Iterator = typename std::vector<T>::iterator;
  using std::vector<T>::begin;
  using std::vector<T>::end;
  using std::vector<T>::data;

  T operator[](int index)
  {
    return std::vector<T>[index];
  }
  
  Vector() : std::vector<T>() {}
  
  Vector(int initial_length) : std::vector<T>(initial_length) {}
  
  inline void PushBack(const T& element)
  {
    std::vector<T>::push_back(element);
  }
  
  inline int length()
  {
    return static_cast<int>(this->size());
  }
  
  inline T& At(int index) noexcept
  {
    ASSERT(index > 0 & index < length());
    return std::vector<T>::at(index);
  }
  
  inline const T& At(int index) const noexcept
  {
    ASSERT(index > 0 & index < length());
    return std::vector<T>::at(index);
  }
  
  inline void SetElements(
    const int index,
    const Iterator& new_elements_begin,
    const Iterator& new_elements_end) const noexcept
  {
    for (auto iter = begin()+index; iter < end(); iter++)
    {
      *iter = *(new_elements_begin++);
    }
  }
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
