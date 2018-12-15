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
  class GenIterator : public std::iterator<std::random_access_iterator_tag, T>
  {
  public:
    typedef typename std::conditional<is_const,const T*,T*>::type Pointer;
    typedef typename std::conditional<is_const,const T&,T&>::type Reference;
    using DifferenceType = typename std::iterator<std::random_access_iterator_tag, T>::difference_type;

    inline bool operator==(const GenIterator& rhs) const noexcept { return ptr_ == rhs.ptr_; }
    inline bool operator!=(const GenIterator& rhs) const noexcept { return ptr_ != rhs.ptr_; }
    inline bool operator>(const GenIterator& rhs) const noexcept { return ptr_ > rhs.ptr_; }
    inline bool operator<(const GenIterator& rhs) const noexcept { return ptr_ < rhs.ptr_; }
    inline bool operator>=(const GenIterator& rhs) const noexcept { return ptr_ >= rhs.ptr_; }
    inline bool operator<=(const GenIterator& rhs) const noexcept { return ptr_ <= rhs.ptr_; }
    
    bool IsWithinBoundsOf(const Vector<T>& vector) const noexcept
    {
      return ptr_ >= vector.data_
        && ptr_ < vector.data_+vector.size_;
    }
    
    GenIterator() : ptr_(nullptr) {}

#ifndef NDEBUG
    GenIterator(T* ptr, const Vector<T>& vector) noexcept : ptr_(ptr) , vector_ref_(vector) {}
    #define VECTOR_REF_ARG , vector_ref_
    #define VECTOR_RHS_REF_ARG , rhs.vector_ref_
#else
    GenIterator(T* ptr) : ptr_(ptr) {}
    #define VECTOR_REF_ARG
    #define VECTOR_RHS_REF_ARG
#endif

    inline Reference operator*() const noexcept { ASSERT(IsWithinBoundsOf(vector_ref_)); return *ptr_; }
    inline Reference operator[](DifferenceType rhs) const noexcept { ASSERT(IsWithinBoundsOf(vector_ref_)); return ptr_[rhs]; }
    inline Pointer operator->() const noexcept { ASSERT(IsWithinBoundsOf(vector_ref_)); return ptr_; }

    inline GenIterator& operator+=(DifferenceType rhs) noexcept { ptr_ += rhs; return *this; }
    inline GenIterator& operator-=(DifferenceType rhs) noexcept { ptr_ -= rhs; return *this; }
    inline GenIterator& operator++() noexcept { ++ptr_; return *this; }
    inline GenIterator& operator--() noexcept { --ptr_;  return *this; }
    inline DifferenceType operator-(const GenIterator& rhs) const noexcept { return ptr_-rhs.ptr_; }
    
    inline GenIterator operator++(int /* unused */) noexcept { GenIterator temp(*this); ++ptr_; return temp; }
    inline GenIterator operator--(int /* unused */) noexcept { GenIterator temp(*this); --ptr_; return temp; }
    inline GenIterator operator+(const GenIterator& rhs) noexcept { return GenIterator(ptr_+rhs.ptr VECTOR_REF_ARG); }
    inline GenIterator operator+(DifferenceType rhs) const noexcept { return GenIterator(ptr_+rhs VECTOR_REF_ARG); }
    inline GenIterator operator-(DifferenceType rhs) const noexcept { return GenIterator(ptr_-rhs VECTOR_REF_ARG); }
    friend inline GenIterator operator+(DifferenceType lhs, const GenIterator& rhs) noexcept { return GenIterator(lhs+rhs.ptr_ VECTOR_RHS_REF_ARG); }
    friend inline GenIterator operator-(DifferenceType lhs, const GenIterator& rhs) noexcept { return GenIterator(lhs-rhs.ptr_ VECTOR_RHS_REF_ARG); }

  private:
    T* ptr_;
#ifndef NDEBUG
    const Vector<T>& vector_ref_;
#endif
  };

  typedef GenIterator<false> Iterator;
  typedef GenIterator<true> ConstIterator;

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
      Deallocate(data_);
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
      Deallocate(data_);
    }
    data_ = nullptr;
  }


  bool OwnsData() const noexcept
  {
    return owns_data_;
  }

#ifndef NDEBUG
    #define THIS_VECTOR_REF_ARG , *this
#else
    #define THIS_VECTOR_REF_ARG
#endif

  Iterator begin() noexcept
  {
    return Iterator(data_ THIS_VECTOR_REF_ARG);
  }


  ConstIterator begin() const noexcept
  {
    return ConstIterator(data_ THIS_VECTOR_REF_ARG);
  }


  Iterator end() noexcept
  {
    return Iterator(data_ + size_ THIS_VECTOR_REF_ARG);
  }


  ConstIterator end() const noexcept
  {
    return ConstIterator(data_ + size_ THIS_VECTOR_REF_ARG);
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
  
  static T* Allocate(size_t size)
  {
    return new T[size];
  }
  
  static void Deallocate(T* data)
  {
    delete[] data;
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
