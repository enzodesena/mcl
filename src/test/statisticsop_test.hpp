/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "statisticsop.hpp"

namespace mcl
{

inline bool StatisticsOpTest()
{
  Vector<Real> vector_c(3);
  vector_c[0] = -0.3;
  vector_c[1] = 0.3;
  vector_c[2] = 2.4;

  Real vector_c_mean = Mean(vector_c);
  ASSERT(IsApproximatelyEqual(vector_c_mean, 0.8, VERY_SMALL));

  ASSERT(IsApproximatelyEqual(Sum(vector_c), -0.3+0.3+2.4, VERY_SMALL));


  Vector<Real> vector_v(3);
  vector_v[0] = 0.1;
  vector_v[1] = -0.5;
  vector_v[2] = 4.0;

  Vector<Real> colonop_b_cmp = Zeros<Real>(4);
  colonop_b_cmp[0] = -1.0;
  colonop_b_cmp[1] = 0.0;
  colonop_b_cmp[2] = 1.0;
  colonop_b_cmp[3] = 2.0;
  
  // Testing Std
  ASSERT(IsApproximatelyEqual(Std(vector_v), 2.443358344574123));
  ASSERT(IsApproximatelyEqual(Std(colonop_b_cmp), 1.290994448735806));

  // Testing var
  ASSERT(IsApproximatelyEqual(Var(vector_v), 5.96999999999999));



  Vector<Real> vector_ba(4);
  vector_ba[0] = -1.2;
  vector_ba[1] = 2.3;
  vector_ba[2] = 3.4;
  vector_ba[3] = -5.0;



  Vector<Real> weights_uniform_a = Multiply<Real>(Ones<Real>(4), 1.0/4.0);
  ASSERT(Mean(vector_ba) == Mean(vector_ba, weights_uniform_a));
  Vector<Real> weights_uniform_b = Multiply<Real>(Ones<Real>(4), 1.0);
  ASSERT(Mean(vector_ba) == Mean(vector_ba, weights_uniform_b));
  Vector<Real> weights_uniform_c = Zeros<Real>(4);
  weights_uniform_c[0] = 0.5;
  weights_uniform_c[2] = 0.5;
  ASSERT(Mean(vector_ba, weights_uniform_c) == 1.1);


  ASSERT(! IsNonNegative(vector_ba));
  Vector<Real> weights_ba_var = Zeros<Real>(4);
  weights_ba_var[0] = 0.2;
  weights_ba_var[1] = 0.3;
  weights_ba_var[2] = 0.6;
  weights_ba_var[3] = 0.5;
  ASSERT(IsApproximatelyEqual(Var(vector_ba, weights_ba_var), 13.319335937500000));


  Vector<Real> vector_aaa(4);
  vector_aaa[0] = 1;
  vector_aaa[1] = 2;
  vector_aaa[2] = 1.5;
  vector_aaa[3] = -1;

  Vector<Real> vector_bbb(4);
  vector_bbb[0] = 0.5;
  vector_bbb[1] = -1;
  vector_bbb[2] = 2;
  vector_bbb[3] = -3;

  ASSERT(IsApproximatelyEqual(Corr(vector_aaa, vector_bbb), 0.689797863572754));


  // Testing entropy
  Vector<Real> pdf_1(4);
  pdf_1[0] = 0.2;
  pdf_1[1] = 0.35;
  pdf_1[2] = 0.15;
  pdf_1[3] = 0.3;

  ASSERT(IsApproximatelyEqual(Entropy(pdf_1, exp(1)), 1.335085165092020));
  ASSERT(IsApproximatelyEqual(Entropy(pdf_1, 2.0), 1.926120746842681));

  Vector<Real> pdf_2(2);
  pdf_2[0] = 0.5;
  pdf_2[1] = 0.5;
  ASSERT(IsApproximatelyEqual(Entropy(pdf_2, exp(1)), 0.693147180559945));
  ASSERT(IsApproximatelyEqual(Entropy(pdf_2, 2.0), 1.0));


  // Testing covariance matrix
  Vector<Real> vector_e(4);
  vector_e[0] = -0.3;
  vector_e[1] = 30.3;
  vector_e[2] = 2.4;
  vector_e[3] = 12.4;

  Vector<Real> vector_f(4);
  vector_f[0] = 2.5;
  vector_f[1] = 1.3;
  vector_f[2] = -2.4;
  vector_f[3] = -1.0;

  ASSERT(! IsEqual(Cov(vector_e, vector_f), Cov(vector_f, vector_e)));
  Matrix<Real> cov_e_f = Cov(vector_e, vector_f);
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(0,0), Var(vector_e)));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(0,0), 191.9800000000000));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(1,1), 4.886666666666667));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(1,1), Var(vector_f)));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(0,0), Var(vector_e)));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(0,1), cov_e_f.GetElement(1,0)));
  ASSERT(IsApproximatelyEqual(cov_e_f.GetElement(0,1), 5.333333333333333));

  return true;
}

} // namespace
