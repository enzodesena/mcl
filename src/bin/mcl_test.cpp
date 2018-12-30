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
#include "digitalfilter.hpp"

#include "basicop_test.ipp"
#include "comparisonop_test.ipp"
#include "elementaryop_test.ipp"
#include "matrixop_test.ipp"
#include "mclintrinsics_test.ipp"
#include "pointwiseop_test.ipp"
#include "point_test.ipp"
#include "quaternion_test.ipp"
#include "random_test.ipp"
#include "statisticsop_test.ipp"
#include "transformop_test.ipp"
#include "vectorop_test.ipp"
#include "rampsmoother_test.ipp"
#include "digitalfilter_test.ipp"

int main (int argc, char * const argv[]) {
#ifndef NDEBUG
  mcl::ComparisonOpTest();
  mcl::RandomGeneratorTest();
  mcl::MclIntrinsicsTest();
  mcl::DigitalFilterTest();
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
  mcl::RampSmootherTest();
  std::cout<<"All tests succeded!\n";
#else
  std::cout<<"Not running tests since NDEBUG is defined and asserts are ignored.\n";
#endif

  mcl::DigitalFilterSpeedTests();
  mcl::IirFilterSpeedTests();
  
  return 0;
}
