/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#ifndef MCL_DATAOP_H
#define MCL_DATAOP_H

#include <vector>
#include "mcltypes.h"

namespace mcl {
  
/** Writes the vector to a file. The separator is endline. */
template<typename T>
void Save(
  const Vector<T>& vector,
  const std::string& file_name,
  const mcl::Int precision = 5)
{
  mcl::Matrix<Real> matrix(vector.length(), 1);
  matrix.SetColumn(0, vector);
  matrix.Save(file_name, precision);
}

} // namespace mcl



#endif
