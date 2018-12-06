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
  using Iterator = typename std::vector<T>::iterator;
  using ConstIterator = typename std::vector<T>::const_iterator;

  Vector() noexcept
    : data_()
    , data_ptr_(nullptr)
    , size_(0)
    , begin_()
    , end_()
    , const_begin_()
    , const_end_()
    , owns_data_(true)
  {
  }

  Vector(
    size_t size,
    T value = T()) noexcept
    : data_(std::vector<T>(size, value))
    , data_ptr_(data_.data())
    , size_(data_.size())
    , begin_(data_.begin())
    , end_(data_.end())
    , const_begin_(data_.begin())
    , const_end_(data_.end())
    , owns_data_(true)
  {
  }

  Vector(
    std::initializer_list<T> list)
    : data_(std::vector<T>(list))
    , data_ptr_(data_.data())
    , size_(data_.size())
    , begin_(data_.begin())
    , end_(data_.end())
    , const_begin_(data_.begin())
    , const_end_(data_.end())
    , owns_data_(true)
  {
  }

  Vector(
    ConstIterator iter_begin,
    ConstIterator iter_end) noexcept
    : data_(std::vector<T>(iter_begin, iter_end))
    , data_ptr_(data_.data())
    , size_(data_.size())
    , begin_(data_.begin())
    , end_(data_.end())
    , const_begin_(data_.begin())
    , const_end_(data_.end())
    , owns_data_(true)
  {
  }

  Vector(
    Vector<T>& referenced_vector,
    size_t first_element_index,
    size_t size) noexcept
    : data_()
    , data_ptr_(referenced_vector.data_.data() + first_element_index)
    , size_(size)
    , begin_(referenced_vector.begin() + (data_ptr_ - referenced_vector.data_.data()))
    , end_(begin_ + size_)
    , const_begin_(referenced_vector.begin() + (data_ptr_ - referenced_vector.data_.data()))
    , const_end_(const_begin_ + size_)
    , owns_data_(false)
  {
    ASSERT(first_element_index+size <= referenced_vector.size());
  }

  Vector(
    const Vector& other)
    : data_(other.data_)
    , data_ptr_((other.owns_data_) ? data_.data() : other.data_ptr_)
    , size_(other.size_)
    , begin_((other.owns_data_) ? data_.begin() : other.begin_)
    , end_((other.owns_data_) ? data_.end() : other.end_)
    , const_begin_((other.owns_data_) ? data_.begin() : other.const_begin_)
    , const_end_((other.owns_data_) ? data_.end() : other.const_end_)
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
    if (this != &other)
    {
      if (owns_data_ && other.data_.data() == data_.data()) { return *this; }

      data_ = other.data_;
      data_ptr_ = (other.owns_data_) ? data_.data() : other.data_ptr_;
      size_ = other.size_;
      begin_ = (other.owns_data_) ? data_.begin() : other.begin_;
      const_begin_ = (other.owns_data_) ? data_.begin() : other.const_begin_;
      end_ = (other.owns_data_) ? data_.end() : other.end_;
      const_end_ = (other.owns_data_) ? data_.end() : other.const_end_;
      owns_data_ = other.owns_data_;
    }
    return *this;
  }

  ~Vector() { data_ptr_ = nullptr; }

  bool OwnsData() const noexcept { return owns_data_; }

  Iterator begin() noexcept { return begin_; }

  ConstIterator begin() const noexcept { return begin_; }

  Iterator end() noexcept { return end_; }

  ConstIterator end() const noexcept { return end_; }

  size_t size() const noexcept { return size_; }

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

private:
  std::vector<T> data_;

  T* data_ptr_;
  size_t size_;

  Iterator begin_;
  Iterator end_;
  ConstIterator const_begin_;
  ConstIterator const_end_;

  bool owns_data_;
};

template<typename TOrigin, typename TDestination>
Vector<TDestination> Cast(
  const Vector<TOrigin>& vector) noexcept
{
  const size_t length = vector.size();
  Vector<TDestination> output(length);
  for (size_t i = 0; i < length; ++i) { output[i] = static_cast<TDestination>(vector[i]); }
  return output;
}

template<typename T>
Vector<Complex<T>> CastToComplex(
  Vector<T> input) noexcept
{
  Vector<Complex<T>> output(input.size());
  for (Int i = 0; i < (Int)input.size(); ++i) { output[i] = Complex<T>(input[i], 0.0); }
  return output;
}
} /**< namespace mcl */
