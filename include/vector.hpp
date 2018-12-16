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
private:
  std::vector<T> data_;
  T* data_ptr_;
  size_t size_;
  bool owns_data_;
  
  
  Vector(
    Vector<T>& other,
    size_t first_element_index,
    size_t size) noexcept
    : data_()
    , data_ptr_(other.data_.data() + first_element_index)
    , size_(size)
    , owns_data_(false)
  {
    ASSERT(first_element_index+size <= other.size());
  }
  
  
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
    
    GenIterator() : ptr_(nullptr) {}

    inline GenIterator& operator+=(DifferenceType rhs) noexcept { ptr_ += rhs; return *this; }
    inline GenIterator& operator-=(DifferenceType rhs) noexcept { ptr_ -= rhs; return *this; }
    inline GenIterator& operator++() noexcept { ++ptr_; return *this; }
    inline GenIterator& operator--() noexcept { --ptr_;  return *this; }
    inline DifferenceType operator-(const GenIterator& rhs) const noexcept { return ptr_-rhs.ptr_; }
    inline GenIterator operator++(int /* unused */) noexcept { GenIterator temp(*this); ++ptr_; return temp; }
    inline GenIterator operator--(int /* unused */) noexcept { GenIterator temp(*this); --ptr_; return temp; }

    bool IsWithinBoundsOf(const Vector<T>& vector) const noexcept
    {
      return ptr_ >= vector.data_ptr_ && ptr_ < vector.data_ptr_+vector.size_;
    }
    
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
  
  
public:
  typedef GenIterator<false> Iterator;
  typedef GenIterator<true> ConstIterator;

  Vector() noexcept
    : data_()
    , data_ptr_(nullptr)
    , size_(0)
    , owns_data_(true)
  {
  }


  Vector(
    size_t size,
    T value = T()) noexcept
    : data_(size, value)
    , data_ptr_(data_.data())
    , size_(size)
    , owns_data_(true)
  {
  }


  Vector(
    std::initializer_list<T> list)
    : data_(list)
    , data_ptr_(data_.data())
    , size_(list.size())
    , owns_data_(true)
  {
  }


  Vector(
    ConstIterator other_begin,
    ConstIterator other_end) noexcept
    : data_(other_begin, other_end)
    , data_ptr_(data_.data())
    , size_(data_.size())
    , owns_data_(true)
  {
  }


  Vector(
    Iterator other_begin,
    Iterator other_end) noexcept
    : data_(other_begin, other_end)
    , data_ptr_(data_.data())
    , size_(data_.size())
    , owns_data_(true)
  {
  }


  Vector(
    const Vector& other)
    : data_(other.data_)
    , data_ptr_((other.owns_data_) ? data_.data() : other.data_ptr_)
    , size_(other.size_)
    , owns_data_(other.owns_data_)
  {
  }


  /** Copy assignment operator. If you are trying to assign the object onto
   itself, this operator has no effect. Also, there is no effect if you try
   to assign a buffer that is referencing itself. For
   instance, if A is a buffer that owns the data, and B is a buffer that
   refers to A's data, then the assignment A = B has no effect. */
  Vector& operator=(
    const Vector& other)
  {
    if (owns_data_ && other.data_.data() == data_.data())
    {
      return *this;
    }

    data_ = other.data_;
    data_ptr_ = (other.owns_data_) ? data_.data() : other.data_ptr_;
    size_ = other.size_;
    
    return *this;
  }


  ~Vector()
  {
    data_ptr_ = nullptr;
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
    return Iterator(data_ptr_ THIS_VECTOR_REF_ARG);
  }


  ConstIterator begin() const noexcept
  {
    return ConstIterator(data_ptr_ THIS_VECTOR_REF_ARG);
  }


  Iterator end() noexcept
  {
    return Iterator(data_ptr_ + size_ THIS_VECTOR_REF_ARG);
  }


  ConstIterator end() const noexcept
  {
    return ConstIterator(data_ptr_ + size_ THIS_VECTOR_REF_ARG);
  }


  size_t size() const noexcept
  {
    return size_;
  }


  T& operator[](
    const size_t index) noexcept
  {
    ASSERT(index>=0 && index < size());
    return data_ptr_[index];
  }


  const T& operator[](
    const size_t index) const noexcept
  {
    ASSERT(index>=0 && index < size());
    return data_ptr_[index];
  }

  friend Vector MakeReference(
    Vector<T>& vector,
    size_t first_element_index = 0,
    size_t size = std::numeric_limits<size_t>::max()) noexcept
  {
    size = (size < std::numeric_limits<size_t>::max()) ? size : vector.size();
    return std::move(Vector(vector, first_element_index, size));
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
