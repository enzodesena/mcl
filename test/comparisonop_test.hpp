/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "comparisonop.hpp"
#include <cassert>

namespace mcl
{
  

inline bool ComparisonOpTest()
{
  ASSERT(IsEqual(0.2, 0.2));
  ASSERT(IsEqual(-0.25, -0.25));
  
  ASSERT(IsSmallerOrEqual(0.1, 0.5));
  ASSERT(IsSmallerOrEqual(-0.5, -0.1));
  
  Vector<Real> vector_a(3);
  vector_a[0] = 0.1;
  vector_a[1] = -0.1;
  vector_a[2] = -2.6;
  Vector<Real> vector_b(3);
  vector_b[0] = 0.1;
  vector_b[1] = -0.5;
  vector_b[2] = -3.0;
  
  ASSERT(AreAllSmallerOrEqual(vector_b, vector_a));
  ASSERT(!AreAllSmallerOrEqual(vector_a, vector_b));
  
  Vector<Real> vector_c(3);
  vector_c[0] = 0.1;
  vector_c[1] = -0.5;
  vector_c[2] = 2.1;
  
  ASSERT(! AreAllSmallerOrEqual(vector_c, vector_a));
  
  
  
  /* Returns true if num is nan */
  ASSERT(IsNan(NAN));
  ASSERT(!IsNan(1.0));
  ASSERT(!IsNan(-1.0));
  
//  /* Returns true if num is nan */
//  Vector<Real> input;
//  input.PushBack(1.0);
//  input.PushBack(NAN);
//  Vector<bool> output_a = IsNan(input);
//  ASSERT(!output_a[0]);
//  ASSERT(output_a[1]);
//  
//  /* Returns opposite bool as input */
//  Vector<bool> output_b = Not(output_a);
//  ASSERT(output_b[0]);
//  ASSERT(!output_b[1]);
//  
//  /* Returns true if all bools are true */
//  ASSERT(!AreAllTrue(output_b));
//  ASSERT(Any(output_b));
//  output_b[1] = true;
//  ASSERT(AreAllTrue(output_b));
//  ASSERT(Any(output_b));
//  ASSERT(!Any(Not(output_b)));
//  
//  Vector<bool> output_c;
//  output_c.PushBack(false);
//  output_c.PushBack(false);
//  output_c.PushBack(false);
//  ASSERT(None(output_c));
//  
//  output_c[0] = true;
//  ASSERT(!None(output_c));
//  output_c[1] = true;
//  output_c[2] = true;
//  ASSERT(!None(output_c));
//  
//  // Testing infs
//  ASSERT(IsInf(INFINITY));
//  ASSERT(IsInf(-INFINITY));
//  
//  Vector<Real> input_d;
//  input_d.PushBack(INFINITY);
//  input_d.PushBack(-INFINITY);
//  input_d.PushBack(0.1);
//  
//  Vector<bool> output_d = IsInf(input_d);
//  ASSERT(output_d[0] == true);
//  ASSERT(output_d[1] == true);
//  ASSERT(output_d[2] == false);
  
  
  return true;
}

  
} // namespace mcl
