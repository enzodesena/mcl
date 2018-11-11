/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once
#include <array>
#include <vector>


constexpr int kDynamicLength = -1;
constexpr int kReferenced = -2;


template<typename T = double, int vector_length = kDynamicLength>
class Vector : std::array<T, vector_length> {
public:
  using std::array<T, vector_length>::operator[];
  
  Vector() {
    puts(__PRETTY_FUNCTION__);
  }
  
  Vector(int length) { ASSERT(length == vector_length); }
  
  inline int length() {
    return vector_length;
  }
  
  inline T At(int index) {
    ASSERT(index>=0 & index<this->size());
    return this[index];
  }
};


template<typename T>
class Vector<T,kDynamicLength> : std::vector<T> {
public:
  using std::vector<T>::operator[];

  Vector(int initial_length) : std::vector<T>(initial_length) {
    puts(__PRETTY_FUNCTION__);
  }
  
  using std::vector<T>::data;
  
  inline int length() {
    return this->size();
  }
  
  inline T At(int index) {
    ASSERT(index>=0 & index<this->size());
    return this[index];
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
