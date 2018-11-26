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
#include "mcltypes.h"
#include "matrix.hpp"
#include "comparisonop.h"
#include "basicop.h"
#include "elementaryop.h"

#if MCL_LOAD_EIGEN
  #include <Eigen/Dense>
#endif

namespace mcl {

// Forward declaration
Vector<std::string> Split(
  const std::string& string,
  char delim) noexcept;
  
  
template<class T>
void Print(
  const Matrix<T>& matrix) noexcept
{
  for (size_t i=0; i<matrix.num_rows(); ++i)
  {
    for (size_t j=0; j<matrix.num_columns(); ++j)
    {
      std::cout<<matrix.GetElement(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }
}

/** Transposes the matrix. Equivalent to Matlab's matrix' */
template<class T>
Matrix<T> Transpose(
  const Matrix<T>& matrix) noexcept
{
  Matrix<T> output(matrix.num_columns(), matrix.num_rows());
  for (size_t i=0; i<output.num_rows(); ++i)
  {
    for (size_t j=0; j<output.num_columns(); ++j)
    {
      output.SetElement(i, j, matrix.GetElement(j, i));
    }
  }
  return output;
}
  
/** 
 Multiplies all the elements of `matrix` by `value`. Equivalent
 to Matlabs' matrix.*value
 */
template<class T>
Matrix<T> Multiply(
  const Matrix<T>& matrix, T value) noexcept
{
  Matrix<T> output(matrix.num_rows(), matrix.num_columns());
  for (size_t i=0; i<output.num_rows(); ++i)
  {
    for (size_t j=0; j<output.num_columns(); ++j)
    {
      output.SetElement(i, j, matrix.GetElement(i, j)*value);
    }
  }
  return output;
}
 
  
/** Matrix multiplication. Equivalent to Matlabs' matrix_a*matrix_b */
template<class T>
Matrix<T> Multiply(
  const Matrix<T>& matrix_a,
  const Matrix<T>& matrix_b) noexcept
{
  ASSERT(matrix_a.num_columns() == matrix_b.num_rows());
  
  Matrix<T> output(matrix_a.num_rows(), matrix_b.num_columns());
  for (size_t i=0; i<output.num_rows(); ++i)
  {
    for (size_t j=0; j<output.num_columns(); ++j)
    {
      T output_value = (T) 0.0;
      for (size_t k=0; k<matrix_a.num_columns(); ++k)
      {
        output_value += matrix_a.GetElement(i, k) * matrix_b.GetElement(k, j);
      }
      output.SetElement(i, j, output_value);
    }
  }
  return output;
}
  
template<class T>
Vector<T> Multiply(
  const Matrix<T>& matrix_a,
  const Vector<T>& vector) noexcept
{
  ASSERT(matrix_a.num_columns() == vector.length());
  Matrix<T> temp_input((Int) vector.length(), 1);
  temp_input.SetColumn(0, vector);
  Matrix<T> temp_output = Multiply(matrix_a, temp_input);
  ASSERT(temp_output.num_columns() == 1);
  ASSERT(temp_output.num_rows() == vector.length());
  return temp_output.GetColumn(0);
}
  
  
/** 
 Extract the maximum value of the matrix. Equivalent to Matlab's
 max(max(matrix)) 
 */
template<class T>
T Max(
  const Matrix<T>& matrix) noexcept
{
  return Max<T>(matrix.Serial());
}
  
  
/** Contains eigenvalues and eigenvectors */
template<typename T>
struct EigOutput
{
  Vector<Complex<T>> eigen_values; /**< Eigenvalues */
  Vector<Vector<Complex<T>>> eigen_vectors; /**< Eigenvectors */
};
  
template<typename T>
Matrix<T> RealPart(const Matrix<Complex<T>>& input) noexcept
{
  Matrix<Real> output(input.num_rows(), input.num_columns());
  for (Int i=0; i<input.num_rows(); ++i)
  {
    for (Int j=0; j<input.num_columns(); ++j)
    {
      output.SetElement(i, j, input.GetElement(i, j).real());
    }
  }
  return output;
}
  
#if MCL_LOAD_EIGEN
template<typename T>
Eigen::MatrixXd ConvertToEigen(
  const Matrix<T>& input)
{
  Eigen::MatrixXd output(input.num_rows(), input.num_columns());
  for (Int i=0; i<input.num_rows(); ++i)
  {
    for (Int j=0; j<input.num_columns(); ++j)
    {
      output(i, j) = input.element(i, j);
    }
  }
  return output;
}

  

template<typename T>
EigOutput Eig(
  const Matrix<T>& matrix)
{
  ASSERT(matrix.num_columns() == matrix.num_rows());
  
  const Int N = matrix.num_columns();
  EigOutput output;
  output.eigen_values = Vector<Complex<T>>(N);
  output.eigen_vectors = std::vector<Vector<Complex<T>> >(N);
  
  // The following constructor triggers compute()
  Eigen::EigenSolver<Eigen::MatrixXd> es(ConvertToEigen(matrix));
  
  for (Int value_i=0; value_i<N; ++value_i)
  {
    output.eigen_values[value_i] = es.eigenvalues()[value_i];
    Vector<Complex<T>> eigen_vector(N);
    for (Int i=0; i<N; ++i)
    {
      eigen_vector[i] = es.eigenvectors()(i, value_i);
    }
    output.eigen_vectors[value_i] = eigen_vector;
  }
  return output;
}
#endif
  
  
template<class T>
bool IsEqual(
  const Matrix<T>& matrix_a,
  const Matrix<T>& matrix_b) noexcept
{
  if (matrix_a.num_rows() != matrix_b.num_rows() |
      matrix_a.num_columns() != matrix_b.num_columns())
  {
    return false;
  }
  
  for (size_t i=0; i<matrix_a.num_rows(); ++i)
  {
    for (size_t j=0; j<matrix_a.num_columns(); ++j)
    {
      if (!IsEqual(matrix_a.GetElement(i, j), matrix_b.GetElement(i, j))) {
        return  false;
      }
    }
  }
  return true;
}


template<typename T>
Matrix<T> Cov(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  Vector<Vector<T>> input(2);
  input[0] = x;
  input[1] = y;
  return Cov(input);
}

template<typename T>
Matrix<T> Cov(
  const Vector<Vector<T> >& input) noexcept
{
  const size_t N = input.length();
  Matrix<T> output(N, N);
  for (size_t i=0; i<N; ++i)
  {
    for (size_t j=0; j<N; ++j)
    {
      output.SetElement(i, j, CovElement(input[i], input[j]));
    }
  }
  return output;
}
  
template<typename T>
T CovElement(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  ASSERT(x.length() == y.length());
  const size_t N = x.length();
  
  T output = Sum(Multiply(AddScalar(x, -Mean(x)), AddScalar(y, -Mean(y))));
  // In case N>1 use the unbiased estimator of covariance.
  output = (N > 1) ? output/((T) (N-1)) : output/((T) (N));
  return output;
}


/** Writes the vector to a file. The separator is endline. */
template<typename T>
inline void Save(
  const Vector<T>& vector,
  const std::string& file_name,
  const mcl::Int precision = 5)
{
  mcl::Matrix<Real> matrix(vector.length(), 1);
  matrix.SetColumn(0, vector);
  matrix.Save(file_name, precision);
}

bool MatrixOpTest();
  
  
} // namespace mcl
