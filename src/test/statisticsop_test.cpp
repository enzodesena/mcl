/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "transformop.h"
#include "comparisonop.h"
#include "statisticsop.h"
#include <vector>
#include <cassert>

namespace mcl {
  
bool StatisticsOpTest() {
  
  Vector<Real> vector_a(4);
  vector_a[0] = 1;
  vector_a[1] = 2;
  vector_a[2] = 1.5;
  vector_a[3] = -1;
  
  Vector<Real> vector_b(4);
  vector_b[0] = 0.5;
  vector_b[1] = -1;
  vector_b[2] = 2;
  vector_b[3] = -3;
  
  ASSERT(IsApproximatelyEqual(Corr(vector_a, vector_b), 0.689797863572754));
  
  return true;
}
  
} // namespace mcl
