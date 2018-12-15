/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include <vector>
#include "mcltypes.hpp"

namespace mcl
{
template<typename T>
class Vector
{
public:

  
  template<bool is_const = false>
  class ConstNoConstIterator : public std::iterator<std::random_access_iterator_tag, T>
  {
  public:
    typedef ConstNoConstIterator SelfType;
    typedef typename std::conditional<is_const,const T*,T*>::type Pointer;
    typedef typename std::conditional<is_const,const T&,T&>::type Reference;
    using DifferenceType = typename std::iterator<std::random_access_iterator_tag, T>::difference_type;

    ConstNoConstIterator() : ptr_(nullptr) {}
    ConstNoConstIterator(T* rhs) : ptr_(rhs) {}
    inline SelfType& operator+=(DifferenceType rhs) noexcept { ptr_ += rhs; return *this; }
    inline SelfType& operator-=(DifferenceType rhs) noexcept { ptr_ -= rhs; return *this; }
    inline Reference operator*() const noexcept { return *ptr_; }
    inline Reference operator[](DifferenceType rhs) const noexcept { return ptr_[rhs]; }
    inline Pointer operator->() const noexcept { return ptr_; }

    inline SelfType& operator++() noexcept { ++ptr_; return *this; }
    inline SelfType& operator--() noexcept { --ptr_; return *this; }
    inline SelfType operator++(int /* unused */) noexcept { SelfType temp(*this); ++ptr_; return temp; }
    inline SelfType operator--(int /* unused */) noexcept { SelfType temp(*this); --ptr_; return temp; }
    inline SelfType operator+(const SelfType& rhs) noexcept { return SelfType(ptr_+rhs.ptr); }
    inline DifferenceType operator-(const SelfType& rhs) const noexcept { return ptr_-rhs.ptr_; }
    inline SelfType operator+(DifferenceType rhs) const noexcept { return SelfType(ptr_+rhs); }
    inline SelfType operator-(DifferenceType rhs) const noexcept { return SelfType(ptr_-rhs); }
    friend inline SelfType operator+(DifferenceType lhs, const SelfType& rhs) noexcept { return SelfType(lhs+rhs.ptr_); }
    friend inline SelfType operator-(DifferenceType lhs, const SelfType& rhs) noexcept { return SelfType(lhs-rhs.ptr_); }

    inline bool operator==(const SelfType& rhs) const noexcept { return ptr_ == rhs.ptr_; }
    inline bool operator!=(const SelfType& rhs) const noexcept { return ptr_ != rhs.ptr_; }
    inline bool operator>(const SelfType& rhs) const noexcept { return ptr_ > rhs.ptr_; }
    inline bool operator<(const SelfType& rhs) const noexcept { return ptr_ < rhs.ptr_; }
    inline bool operator>=(const SelfType& rhs) const noexcept { return ptr_ >= rhs.ptr_; }
    inline bool operator<=(const SelfType& rhs) const noexcept { return ptr_ <= rhs.ptr_; }
  private:
    T* ptr_;
  };

  typedef ConstNoConstIterator<true> ConstIterator;
  typedef ConstNoConstIterator<false> Iterator;

  Vector() noexcept
    : size_(0)
    , data_(nullptr)
    , owns_data_(true)
  {
  }


  Vector(
    size_t size,
    T value = T()) noexcept
    : size_(size)
    , data_(Allocate(size_))
    , owns_data_(true)
  {
    for (size_t i = 0; i<size_; ++i)
    {
      data_[i] = value;
    }
  }


  Vector(
    std::initializer_list<T> list)
    : size_(list.size())
    , data_(Allocate(size_))
    , owns_data_(true)
  {
    auto iter = begin();
    for (auto& element : list)
    {
      *iter++ = element;
    }
  }


  Vector(
    ConstIterator other_begin,
    ConstIterator other_end) noexcept
    : size_(other_end-other_begin)
    , data_(Allocate(size_))
    , owns_data_(true)
  {
    auto iter = begin();
    while(other_begin != other_end)
    {
      *iter++ = *other_begin++;
    }
  }


  Vector(
    Iterator other_begin,
    Iterator other_end) noexcept
    : size_(other_end-other_begin)
    , data_(Allocate(size_))
    , owns_data_(true)
  {
    auto iter = begin();
    while(other_begin != other_end)
    {
      *iter++ = *other_begin++;
    }
  }


  Vector(
    const Vector& other)
    : size_(other.size_)
    , data_((other.owns_data_) ? Allocate(size_) : other.data_)
    , owns_data_(other.owns_data_)
  {
    if (owns_data_)
    {
      auto iter = begin();
      auto other_iter = other.begin();
      while(other_iter != other.end())
      {
        *iter++ = *other_iter++;
      }
    }
  }


  /** Copy assignment operator. If you are trying to assign the object onto
   itself, this operator has no effect. Also, there is no effect if you try
   to assign a buffer that is referencing itself. For
   instance, if A is a buffer that owns the data, and B is a buffer that
   refers to A's data, then the assignment A = B has no effect. */
  Vector& operator=(
    const Vector& other)
  {
    if (owns_data_ && other.data_ == data_)
    {
      return *this;
    }
    
    if (owns_data_)
    {
      Deallocate();
    }
    
    size_ = other.size_;
    if (other.owns_data_)
    {
      data_ = Allocate(size_);
      auto iter = begin();
      auto other_iter = other.begin();
      while(other_iter != other.end())
      {
        *iter++ = *other_iter++;
      }
    }
    else
    {
      data_ = other.data_;
    }
    
    owns_data_ = other.owns_data_;
    
    return *this;
  }


  ~Vector()
  {
    if (OwnsData())
    {
      Deallocate();
    }
    data_ = nullptr;
  }


  bool OwnsData() const noexcept
  {
    return owns_data_;
  }


  Iterator begin() noexcept
  {
    return Iterator(data_);
  }


  ConstIterator begin() const noexcept
  {
    return ConstIterator(data_);
  }


  Iterator end() noexcept
  {
    return Iterator(data_ + size_);
  }


  ConstIterator end() const noexcept
  {
    return ConstIterator(data_ + size_);
  }


  size_t size() const noexcept
  {
    return size_;
  }


  T& operator[](
    const size_t index) noexcept
  {
    ASSERT(index>=0 && index < size());
    return data_[index];
  }


  const T& operator[](
    const size_t index) const noexcept
  {
    ASSERT(index>=0 && index < size());
    return data_[index];
  }

  friend Vector MakeReference(
    Vector<T>& vector,
    size_t first_element_index = 0,
    size_t size = std::numeric_limits<size_t>::max()) noexcept
  {
    size = (size < std::numeric_limits<size_t>::max()) ? size : vector.size();
    return std::move(Vector(vector, first_element_index, size));
  }
  
private:
  size_t size_;
  T* data_;
  bool owns_data_;
  
  T* Allocate(size_t size)
  {
    return new T[size];
  }
  
  void Deallocate()
  {
    delete[] data_;
  }
  
  Vector(
    Vector<T>& referenced_vector,
    size_t first_element_index,
    size_t size) noexcept
    : size_(size)
    , data_(referenced_vector.data_ + first_element_index)
    , owns_data_(false)
  {
    ASSERT(first_element_index+size <= referenced_vector.size());
  }
};


template<typename TOrigin, typename TDestination>
Vector<TDestination> Cast(
  const Vector<TOrigin>& vector) noexcept
{
  const size_t length = vector.size();
  Vector<TDestination> output(length);
  for (size_t i = 0; i < length; ++i)
  {
    output[i] = static_cast<TDestination>(vector[i]);
  }
  return output;
}


template<typename T>
Vector<Complex<T>> CastToComplex(
  Vector<T> input) noexcept
{
  Vector<Complex<T>> output(input.size());
  for (Int i = 0; i < (Int)input.size(); ++i)
  {
    output[i] = Complex<T>(input[i], 0.0);
  }
  return output;
}
} /**< namespace mcl */
