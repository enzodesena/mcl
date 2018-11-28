/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include <iostream>
#include "vectorop.hpp"
#include "randomop.hpp"
#include "point.hpp"
#include "elementaryop.hpp"
#include "quaternion.hpp"
#include "basicop.hpp"
#include "matrixop.hpp"
#include "pointwiseop.hpp"
#include "mclintrinsics.hpp"
#include "statisticsop.hpp"
#include "transformop.hpp"
#include "firfilter.hpp"
#include "iirfilter.hpp"

#include "basicop_test.hpp"
#include "comparisonop_test.hpp"
#include "elementaryop_test.hpp"
#include "matrixop_test.hpp"
#include "mclintrinsics_test.hpp"
#include "pointwiseop_test.hpp"
#include "point_test.hpp"
#include "quaternion_test.hpp"
#include "random_test.hpp"
#include "statisticsop_test.hpp"
#include "transformop_test.hpp"
#include "vectorop_test.hpp"

#include "digitalfilter_test.hpp"

int main (int argc, char * const argv[]) {
#ifndef NDEBUG
  mcl::ComparisonOpTest();
  mcl::RandomGeneratorTest();
  mcl::MclIntrinsicsTest();
  mcl::FirFilterTest();
  mcl::QuaternionTest();
  mcl::ElementaryOpTest();
  mcl::BasicOpTest();
  mcl::VectorOpTest();
  mcl::PointWiseOpTest();
  mcl::TransformOpTest();
  mcl::MatrixOpTest();
  mcl::PointTest();
  mcl::IirFilterTest();
  mcl::StatisticsOpTest();
  std::cout<<"All tests succeded!\n";
#else
  std::cout<<"Not running tests since NDEBUG is defined and asserts are ignored.\n";
#endif

  mcl::FirFilterSpeedTests();
  
  return 0;
}
