/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "quaternion.hpp"
#include "comparisonop.hpp"
#include "mcltypes.hpp"
#include "point.hpp"
#include "vectorop.hpp"
#include "constants.hpp"
#include <vector>

namespace mcl
{

bool QuaternionTest()
{
  
  Quaternion q(0.1, 0.2, 0.3, 0.4);
  ASSERT(IsEqual(q.w(), 0.1));
  ASSERT(IsEqual(q.x(), 0.2));
  ASSERT(IsEqual(q.y(), 0.3));
  ASSERT(IsEqual(q.z(), 0.4));
  ASSERT(IsEqual(q, Quaternion(0.1, 0.2, 0.3, 0.4)));
  
  Quaternion p(1.0, 0.0, 1.0, 0.0);
  Quaternion r(1.0, 0.5, 0.5, 0.75);
  ASSERT(IsEqual(QuatMultiply(p, r), Quaternion(0.5, 1.25, 1.5, 0.25)));
  ASSERT(! IsEqual(QuatMultiply(r, p), Quaternion(0.5, 1.25, 1.5, 0.25)));
  Quaternion s(2.0, 1.0, 0.1, 0.1);
  ASSERT(IsEqual(QuatMultiply(p, s), Quaternion(1.9, 1.1, 2.1, -0.9)));
  ASSERT(! IsEqual(QuatMultiply(s, p), Quaternion(1.9, 1.1, 2.1, -0.9)));
  
  // Quaternion multiplication is distributive but not commutative
  ASSERT(IsApproximatelyEqual(QuatMultiply(QuatMultiply(q, p), r), QuatMultiply(q, QuatMultiply(p, r))));
  
  ASSERT(IsApproximatelyEqual(QuatRotate(Quaternion(0.5,0.5,0.5,0.5), Point(1.0, 2.0, 3.0), kRightHanded),
                 Point(3.0, 1.0, 2.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(p, Point(1.0,1.0,1.0), kLeftHanded), Point(-1.0,1.0,1.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(p, Point(2.0,3.0,4.0), kLeftHanded), Point(-4.0,3.0,2.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(Quaternion(1.0, 0.5, 0.3, 0.1), Point(1.0,1.0,1.0), kLeftHanded),
                 Point(0.8519,1.4741,0.3185), 1.0E-3));
  
  ASSERT(IsApproximatelyEqual(AxAng2Quat(1.0,0.0,0.0,PI/2.0),
                 Quaternion(sqrt(2.0)/2.0,sqrt(2.0)/2.0,0.0,0.0)));
  
  ASSERT(IsApproximatelyEqual(QuatRotate(AxAng2Quat(0.0,1.0,0.0,-PI/2.0), Point(1.0,0.0,0.0), kRightHanded),
                 Point(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(QuatConj(AxAng2Quat(0.0,1.0,0.0,PI/2.0)), Point(1.0,0.0,0.0), kRightHanded),
                 Point(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(AxAng2Quat(0.0,1.0,0.0,PI/2.0), Point(1.0,0.0,0.0), kLeftHanded),
                 Point(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(AxAng2Quat(0.0,1.0,0.0,PI/2.0), Point(-1.0,0.0,0.0), kRightHanded),
                 Point(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(QuatRotate(AxAng2Quat(0.0,1.0,0.0,-PI/2.0), Point(-1.0,0.0,0.0), kLeftHanded),
                 Point(0.0,0.0,1.0)));
  
  // Test identity quaternion
  ASSERT(IsApproximatelyEqual(QuatRotate(Quaternion<Real>::Identity(), Point(1.0, 2.0, 3.0)), Point(1.0, 2.0, 3.0)));
  ASSERT(IsApproximatelyEqual(QuatMultiply(Quaternion<Real>::Identity(), p), p));
  ASSERT(IsApproximatelyEqual(QuatMultiply(Quaternion<Real>::Identity(), r), r));
  ASSERT(IsApproximatelyEqual(QuatMultiply(Quaternion<Real>::Identity(), s), s));
  
  // Test QuatInverse
  ASSERT(IsApproximatelyEqual(QuatInverse(Quaternion(1.0,0.0,1.0,0.0)), Quaternion(0.5,0.0,-0.5,0.0)));
  
  AxAng axang1 = Quat2AxAng(AxAng2Quat(0.0,0.0,1.0,PI/2.0));
  ASSERT(IsEqual(axang1.x, 0.0));
  ASSERT(IsEqual(axang1.y, 0.0));
  ASSERT(IsEqual(axang1.z, 1.0));
  ASSERT(IsEqual(axang1.angle, PI/2.0));
  
  AxAng axang2 = Quat2AxAng(AxAng2Quat(1.0,0.0,0.0,-0.456));
  ASSERT(IsApproximatelyEqual(axang2.x, -1.0));
  ASSERT(IsApproximatelyEqual(axang2.y, 0.0));
  ASSERT(IsApproximatelyEqual(axang2.z, 0.0));
  ASSERT(IsApproximatelyEqual(axang2.angle, 0.456));
  
  AxAng axang3 = Quat2AxAng(AxAng2Quat(1.0,0.0,0.0,0.456));
  ASSERT(IsApproximatelyEqual(axang3.x, 1.0));
  ASSERT(IsApproximatelyEqual(axang3.y, 0.0));
  ASSERT(IsApproximatelyEqual(axang3.z, 0.0));
  ASSERT(IsApproximatelyEqual(axang3.angle, 0.456));
  
  return true;
}
  
} // namespace mcl
