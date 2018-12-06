/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

// This file contains definitions of matrix operations and classes

#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include "basicop.hpp"

#if MCL_LOAD_EIGEN
  #include <Eigen/Dense>
#endif

namespace mcl
{
/** Matrix class */
template<class T>
class Matrix
{
public:
  /** Default constructor with empty (0x0) matrix */
  Matrix() noexcept
    : num_rows_(0)
    , num_columns_(0)
  {
  }

  /**
   Constructs a matrix with default entries
   */
  Matrix(
    size_t num_rows,
    size_t num_columns) noexcept
    : data_(Vector<Vector<T>>(num_rows))
    , num_rows_(num_rows)
    , num_columns_(num_columns) { for (size_t i = 0; i < num_rows; ++i) { data_[i] = Vector<T>(num_columns); } }

  /**
   Constructs a matrix from a vector of vectors (outer vector represents rows)
   */
  Matrix(
    const Vector<Vector<T>> vectors) noexcept
  {
    num_rows_ = vectors.size();
    if (num_rows_ > 0)
    {
      num_columns_ = vectors[0].size();
      for (size_t i = 1; i < num_rows_; ++i)
      {
        // Check that all rows have the same number of columns
        if ((Int)vectors[i].size() != num_columns_)
        {
          ASSERT_WITH_MESSAGE(false, "One or more rows do not have the same number of columns");
        }
      }
      data_ = vectors;
    }
    else { num_columns_ = 0; }
  }

  /** Sets element in given row and column */
  void SetElement(
    const size_t index_row,
    const size_t index_column,
    const T& element) noexcept
  {
    ASSERT_WITH_MESSAGE(index_row < num_rows_ && index_column < num_columns_,
      "Out-of-bounds index in setting matrix element");
    data_[index_row][index_column] = element;
  }

  /** Sets an entire column */
  void SetColumn(
    size_t index_column,
    Vector<T> column) noexcept
  {
    ASSERT(column.size() == num_rows_);
    ASSERT(index_column < num_columns_);
    for (size_t i = 0; i < num_rows_; ++i) { SetElement(i, index_column, column[i]); }
  }

  /** Sets an entire row */
  void SetRow(
    size_t index_row,
    Vector<T> row) noexcept
  {
    ASSERT(row.size() == num_columns_);
    ASSERT(index_row < num_rows_);
    data_[index_row] = row;
  }

  /** Accesses an element in given row and column */
  T GetElement(
    const size_t index_row,
    const size_t index_column) const noexcept
  {
    ASSERT(index_row < num_rows_);
    ASSERT(index_column < num_columns_);
    return data_[index_row][index_column];
  }

  /** Accesses an entire row */
  Vector<T> GetRow(
    size_t index_row) const noexcept
  {
    ASSERT(index_row < num_rows_);
    return data_[index_row];
  }

  /** Accesses an entire column */
  Vector<T> GetColumn(
    size_t index_column) const noexcept
  {
    ASSERT(index_column < num_columns_);
    Vector<T> output(num_rows_);
    for (size_t i = 0; i < num_rows_; ++i) { output[i] = data_[i][index_column]; }
    return output;
  }

  /** Returns the serialised matrix. Equivalent to Matlab's matrix(:) */
  Vector<T> Serial() const noexcept
  {
    Vector<T> serial(num_columns() * num_rows());

    size_t k = 0;
    for (size_t j = 0; j < num_columns(); ++j)
    {
      for (size_t i = 0; i < num_rows(); ++i) { serial[k++] = GetElement(i, j); }
    }
    return serial;
  }

  /** Returns the number of rows */
  size_t num_rows() const noexcept { return num_rows_; }

  /** Returns the number of columns */
  size_t num_columns() const noexcept { return num_columns_; }

  /** Writes the matrix to a file. The optional parameter `precision` sets
   the number of decimal positions in the output file*/
  void Save(
    std::string file_name,
    Int precision = 5)
  {
    std::ofstream output_file;
    output_file.open(file_name.c_str());
    output_file << std::fixed;
    output_file << std::setprecision((int)precision);
    for (size_t i = 0; i < num_rows_; ++i)
    {
      for (size_t j = 0; j < num_columns_; ++j) { output_file << data_.at(i).at(j) << " "; }
      output_file << std::endl;
    }
    output_file.close();
  }

  /** Returns the raw data */
  Vector<Vector<T>> data() noexcept { return data_; }

  /**
   Reads a matrix. Elements have to be separated by tabs and there
   must be no ending empty line (e.g. 5 lines == 5 rows).
   */
  static Matrix Load(
    std::string file_name)
  {
    std::string line;
    std::ifstream in_file(file_name.c_str());
    ASSERT(in_file.is_open());

    // First: lets count the number of rows
    size_t number_of_rows = 0;
    size_t number_of_columns = 0;
    while (std::getline(in_file, line))
    {
      Vector<std::string> elements = Split(line, '\t');
      if (number_of_columns == 0) { number_of_columns = (Int)elements.size(); }
      else { ASSERT(number_of_columns == (Int)elements.size()); }

      ++number_of_rows;
    }

    // Reset pointer
    in_file.clear();
    in_file.seekg(0, std::ios::beg);

    // Create new matrix
    // TODO: recognize complex matrix
    Matrix<T> matrix(number_of_rows, number_of_columns);
    for (size_t row = 0; row < number_of_rows; ++row)
    {
      std::getline(in_file, line);
      Vector<std::string> elements = Split(line, '\t');
      for (Int column = 0; column < (Int)elements.size(); ++column)
      {
        matrix.SetElement(row, column, (T)StringToDouble(elements[column]));
      }
    }

    in_file.close();
    return matrix;
  }

  /**
   Constructs a matrix of all ones. Equivalent to Matlab's ones(N,M).
   */
  static Matrix Ones(
    size_t number_of_rows,
    size_t number_of_columns) noexcept
  {
    Matrix<T> matrix(number_of_rows, number_of_columns);
    for (size_t row = 0; row < number_of_rows; ++row)
    {
      for (size_t column = 0; column < number_of_columns; ++column) { matrix.SetElement(row, column, (T)1.0); }
    }
    return matrix;
  }

  /**
   Constructs a matrix of all ones. Equivalent to Matlab's ones(N).
   */
  static Matrix Ones(
    size_t matrix_dimension) noexcept { return Matrix<T>::Ones(matrix_dimension, matrix_dimension); }

private:
  // Outer is rows, inner is columns. Hence, data_[0] is the first column.
  Vector<Vector<T>> data_;
  size_t num_rows_;
  size_t num_columns_;
};
} // namespace mcl
