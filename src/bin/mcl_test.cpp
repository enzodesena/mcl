/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@me.com
 */


#include <iostream>
#include "matrixop.h"
#include "vectorop.h"
#include "transformop.h"
#include "statisticsop.h"
#include "firfilter.h"
#include "randomop.h"
#include "iirfilter.h"
#include "exception.h"



int main (int argc, char * const argv[]) {
  using namespace mcl;
  
  Quaternion::Test();
  ElementaryOpTest();
  BasicOpTest();
  VectorOpTest();
  PointWiseOpTest();
  TransformOpTest();
  MatrixOpTest();
  StatisticsOpTest();
  ComparisonOpTest();
  PointTest();
  
  IirFilter::Test();
  FirFilter::SpeedTests();
  FirFilter::Test();
  RandomGenerator::Test();
  
  Exception::ExceptionTest();
  
  
  std::cout<<"All tests succeded!\n";
  
  return 0;
}