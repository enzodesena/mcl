/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include "pointwiseop.h"
#include <vector>

using std::vector;

namespace mcl {
  
/** Returns the Pearson linear correlation between `vector_a` and `vector_b` */
template<typename T>
T Corr(
  const Vector<T>& x,
  const Vector<T>& y) noexcept
{
  
  T pearson_num_lin = Sum(mcl::Multiply(Add(x,-Mean(x)), Add(y,-Mean(y))));
  T pearson_den_lin = sqrt(Sum(Pow(Add(x, -Mean(x)),2.0)))*
  sqrt(Sum(Pow(Add(y, -Mean(y)),2.0)));
  return pearson_num_lin/pearson_den_lin;
}

  
bool StatisticsOpTest();
  
} /**< namespace mcl */

#endif
