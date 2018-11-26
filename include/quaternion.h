/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

// This file contains definitions of matrix operations and classes

#pragma once

#include "mcltypes.h"
#include "point.h"
#include "constants.h"
#include "elementaryop.h"
#include <cassert>

namespace mcl {
  
// Forward declarations
template<typename T>
struct Point;
// End of forward declarations
  
/** Enum describing the angles ordering convention for Euler angles */
enum EulerOrder {
  zxz, xyx, yzy, zyz, xzx, yxy,
  xyz, yzx, zxy, xzy, zyx, yxz
};
  
template<typename T>
struct AxAng {
  T x;
  T y;
  T z;
  T angle;
};

/** Quaternion class */
template<typename T>
class Quaternion
{
public:
  /** Constructs a quaternion, with the first element being the scalar
   component and the following three forming the vector component */
  Quaternion(
    const T w,
    const T x,
    const T y,
    const T z) noexcept : w_(w), x_(x), y_(y), z_(z)
  {
  }
  
  T w() const noexcept { return w_; }
  T x() const noexcept { return x_; }
  T y() const noexcept { return y_; }
  T z() const noexcept { return z_; }
  
  /** Constructs a quaternion that is neutral to rotations, 
   i.e. a multiplicative identity quaternion */
  static Quaternion<T> Identity()
  {
    return Quaternion(1.0, 0.0, 0.0, 0.0);
  }
  
private:
  // q = w + x*i + y*j + z*k where i² = j² = k² = i*j*k = -1
  T w_;
  T x_;
  T y_;
  T z_;
};
  
/** Returns a Quaternion using a given axis-angle representation */
template<typename T>
Quaternion<T> AxAng2Quat(
  const T x,
  const T y,
  const T z,
  const T angle) noexcept
{
  T norm = sqrt(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0));
  return Quaternion(
    cos(angle/2.0),
    sin(angle/2.0)*x/norm,
    sin(angle/2.0)*y/norm,
    sin(angle/2.0)*z/norm);
}

template<typename T>
AxAng<T> Quat2AxAng(
  const Quaternion<T>& input) noexcept
{
  AxAng<T> output;
  
  Quaternion<T> q = (Norm(input)>1.0) ? Quaternion<T>(
    input.w()/Norm(input),
    input.x()/Norm(input),
    input.y()/Norm(input),
    input.z()/Norm(input)) : input;
  
  Point q_r(q.x(), q.y(), q.z());
  if (q_r.norm() == 0.0) // TODO: verify whether this should be an approx equal
  {
    output.x = 0;
    output.y = 0;
    output.z = 1;
    output.angle = 0;
    return output;
  }
  
  T n_inv = 1.0/Sqrt(1.0-q.w()*q.w());
  output.angle = 2.0*acos(q.w());
  output.x = q.x()*n_inv;
  output.y = q.y()*n_inv;
  output.z = q.z()*n_inv;
  return output;
}

template<typename T>
Quaternion<T> QuatConj(
  const Quaternion<T>& q) noexcept
{
  return Quaternion(q.w(), -q.x(), -q.y(), -q.z());
}

/** Returns the norm of a quaternion (defined the same as the Eucledian 
 norm in R^4) */
template<typename T>
T Norm(const Quaternion<T>& q) noexcept
{
  return sqrt(pow(q.w(),2.0)+pow(q.x(),2.0)+pow(q.y(),2.0)+pow(q.z(),2.0));
}
  

/** Implements the (Hamilton) quaternion multiplication **/
template<typename T>
Quaternion<T> QuatMultiply(
  const Quaternion<T>& q,
  const Quaternion<T>& r) noexcept
{
return Quaternion(
  r.w()*q.w()-r.x()*q.x()-r.y()*q.y()-r.z()*q.z(),
  r.w()*q.x()+r.x()*q.w()-r.y()*q.z()+r.z()*q.y(),
  r.w()*q.y()+r.x()*q.z()+r.y()*q.w()-r.z()*q.x(),
  r.w()*q.z()-r.x()*q.y()+r.y()*q.x()+r.z()*q.w());
}
  
template<typename T>
Point<T> QuatRotate(
  const Quaternion<T>& q,
  const Point<T>& r,
  const Handedness handedness = kRightHanded) noexcept
{
  T norm = Norm(q);
  Quaternion q_norm(q.w()/norm, q.x()/norm, q.y()/norm, q.z()/norm);
  
  Quaternion p(0.0, r.x(), r.y(), r.z());
  
  if (handedness == kRightHanded)
  {
    Quaternion result = QuatMultiply(QuatMultiply(q_norm, p), QuatConj(q_norm));
    return Point(result.x(), result.y(), result.z());
  }
  else
  {
    Quaternion result = QuatMultiply(QuatMultiply(QuatConj(q_norm), p), q_norm);
    return Point(result.x(), result.y(), result.z());
  }
}
  
/** Converts Euler angles with a given convention to a Quaternion. 
 Each input angle corresponds to the associated ordering.
 E.g. for zyx convention (which is the default), the first angle is for a
 rotation around the z axis, second for y, etc.
 */
template<typename T>
Quaternion<T> Eul2Quat(
  const T angle_1,
  const T angle_2,
  const T angle_3,
  const EulerOrder order = zyx) noexcept
{
  Quaternion<T> rotation_1 = Quaternion<T>::Identity();
  Quaternion<T> rotation_2 = Quaternion<T>::Identity();
  Quaternion<T> rotation_3 = Quaternion<T>::Identity();

  switch (order)
  {
    case xyx: case xzx: case xyz: case xzy:
      rotation_1 = Quaternion(cos(angle_1/2.0), sin(angle_1/2.0), 0.0, 0.0);
      break;
    case yzy: case yxy: case yzx: case yxz:
      rotation_1 = Quaternion(cos(angle_1/2.0), 0.0, sin(angle_1/2.0), 0.0);
      break;
    case zxz: case zyz: case zxy: case zyx:
      rotation_1 = Quaternion(cos(angle_1/2.0), 0.0, 0.0, sin(angle_1/2.0));
      break;
    default:
      ASSERT(false);
      break;
  }
  
  switch (order)
  {
    case yxy: case yxz: case zxy: case zxz:
      rotation_2 = Quaternion(cos(angle_2/2.0), sin(angle_2/2.0), 0.0, 0.0);
      break;
    case zyz: case xyx: case xyz: case zyx:
      rotation_2 = Quaternion(cos(angle_2/2.0), 0.0, sin(angle_2/2.0), 0.0);
      break;
    case yzy: case yzx: case xzx: case xzy:
      rotation_2 = Quaternion(cos(angle_2/2.0), 0.0, 0.0, sin(angle_2/2.0));
      break;
    default:
      ASSERT(false);
      break;
  }
  
  switch (order)
  {
    case xyx: case zyx: case yzx: case xzx:
      rotation_3 = Quaternion(cos(angle_3/2.0), sin(angle_3/2.0), 0.0, 0.0);
      break;
    case zxy: case yxy: case yzy: case xzy:
      rotation_3 = Quaternion(cos(angle_3/2.0), 0.0, sin(angle_3/2.0), 0.0);
      break;
    case yxz: case zxz: case zyz: case xyz:
      rotation_3 = Quaternion(cos(angle_3/2.0), 0.0, 0.0, sin(angle_3/2.0));
      break;
    default:
      ASSERT(false);
      break;
  }
  
  return QuatMultiply(QuatMultiply(rotation_1, rotation_2), rotation_3);
}
  
/** Returns the Euler angle around the x-axis associated to a given quaternion
 and for a given Euler rotation convention */
template<typename T>
T Quat2EulX(
  const Quaternion<T> q,
  const EulerOrder order = zyx) noexcept
{
  switch (order)
  {
    case zyx:
      return atan2(-2.0*q.y()*q.z()+2.0*q.w()*q.x(),
                   pow(q.w(),2.0)+pow(q.z(),2.0)-pow(q.y(),2.0)-pow(q.x(),2.0));
      break;
    default:
      ASSERT(false);
      return NAN;
      break;
  }
}
  
/** Returns the Euler angle around the y-axis associated to a given quaternion
 and for a given Euler rotation convention */
template<typename T>
T Quat2EulY(
  const Quaternion<T> q,
  const EulerOrder order = zyx) noexcept
{
  switch (order)
  {
    case zyx:
      return asin(2.0*q.x()*q.z()+2.0*q.w()*q.y());
      break;
    default:
      ASSERT(false);
      return NAN;
      break;
  }
}
  
/** Returns the Euler angle around the z-axis associated to a given quaternion
 and for a given Euler rotation convention */
template<typename T>
T Quat2EulZ(
  const Quaternion<T> q,
  const EulerOrder order = zyx) noexcept
{
  switch (order)
  {
    case zyx:
      return atan2(-2.0*q.x()*q.y()+2.0*q.w()*q.z(),
                   pow(q.w(), 2.0)+pow(q.x(), 2.0)
                   -pow(q.y(), 2.0)-pow(q.z(), 2.0));
      break;
    default:
      ASSERT(false);
      return NAN;
      break;
  }
}
  
template<typename T>
Quaternion<T> QuatInverse(
  const Quaternion<T> q) noexcept
{
  T norm_sq = pow(Norm(q), 2.0);
  Quaternion conj = QuatConj(q);
  return Quaternion(
    conj.w()/norm_sq,
    conj.x()/norm_sq,
    conj.y()/norm_sq,
    conj.z()/norm_sq);
}



template<typename T>
inline bool IsEqual(
  const Quaternion<T>& q_a,
  const Quaternion<T>& q_b) noexcept
{
  return q_a.w() == q_b.w() && q_a.x() == q_b.x() &
    q_a.y() == q_b.y() && q_a.z() == q_b.z();
}




template<typename T>
inline bool IsApproximatelyEqual(
  const Quaternion<T>& quat_a,
  const Quaternion<T>& quat_b,
  const T precision = VERY_SMALL)
{
  return
    IsApproximatelyEqual(quat_a.w(), quat_b.w()) &&
    IsApproximatelyEqual(quat_a.x(), quat_b.x()) &&
    IsApproximatelyEqual(quat_a.y(), quat_b.y()) &&
    IsApproximatelyEqual(quat_a.z(), quat_b.z());
}

bool QuaternionTest();
  
} // namespace mcl
