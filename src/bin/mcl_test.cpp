/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include <iostream>
#include "matrixop.hpp"
#include "vectorop.hpp"
#include "vector.hpp"
#include "elementaryop.hpp"
#include "quaternion.hpp"
#include "transformop.hpp"
#include "firfilter.hpp"
#include "randomop.hpp"
#include "iirfilter.hpp"



int main (int argc, char * const argv[]) {
  using namespace mcl;

#ifndef NDEBUG
  FirFilterTest();
  QuaternionTest();
  ElementaryOpTest();
  BasicOpTest();
  VectorOpTest();
  PointWiseOpTest();
  TransformOpTest();
  MatrixOpTest();
  ComparisonOpTest();
  PointTest();
  IirFilterTest();
  RandomGenerator::Test();
  std::cout<<"All tests succeded!\n";
#else
  std::cout<<"Not running tests since NDEBUG is defined and asserts are ignored.\n";
#endif

//  FirFilter::SpeedTests();
  
  return 0;
}
