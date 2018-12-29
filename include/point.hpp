/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "comparisonop.hpp"

namespace mcl
{
template<typename T>
class Triplet
{
public:
  /** Constructs a `Point` with all coordinates set to zero. */
  Triplet() noexcept
  {
  }


  /** Constructor with explicit definition of all coordinates. */
  Triplet(
    const T x,
    const T y,
    const T z) noexcept
    : x_(x)
    , y_(y)
    , z_(z)
  {
  }


  // Getter methods.
  T x() const noexcept
  {
    return x_;
  }


  T y() const noexcept
  {
    return y_;
  }


  T z() const noexcept
  {
    return z_;
  }


  void SetX(
    const T value) noexcept
  {
    x_ = value;
  }


  void SetY(
    const T value) noexcept
  {
    y_ = value;
  }


  void SetZ(
    const T value) noexcept
  {
    z_ = value;
  }


  /** Returns the norm of the vector, or, in other words, the distance
   of the point from the origin (0,0,0) */
  T norm() const noexcept
  {
    return sqrt(pow(x(), 2) + pow(y(), 2) + pow(z(), 2));
  }


  /** Returns the angle formed with the z-axis */
  T theta() const noexcept
  {
    return (T)acos(z() / norm());
  }


  /** Returns the angle formed between the projection on the x-y axis and
   the x-axis */
  T phi() const noexcept
  {
    return (Real)atan2(y(), x());
  }


  /** Modifies the point (i.e. vector) such that its norm is equal to 1. */
  void Normalize() noexcept
  {
    T vector_norm = norm();
    x_ = x_ / vector_norm;
    y_ = y_ / vector_norm;
    z_ = z_ / vector_norm;
  }


private:
  T x_;
  T y_;
  T z_;
};


template<typename T>
struct Point : public Triplet<T>
{
  Point() noexcept
  {
  }


  /** Constructor with explicit definition of all coordinates. */
  Point(
    const T x,
    const T y,
    const T z) noexcept
    : Triplet<T>(x, y, z)
  {
  }
};


/**
 Rotates the reference system about the x-axis with the right-hand rule.
 e.g. RotateAboutX(Point(0.0,1.0,0.0), pi/2) == Point(0.0,0.0,1.0)
 */
template<typename T>
Point<T> RotateAboutX(
  const Point<T>& point,
  const T angle) noexcept
{
  T cos_angle = cos(angle);
  T sin_angle = sin(angle);

  return Point
  (
    point.x(),
    point.y() * cos_angle + point.z() * (-sin_angle),
    point.y() * (sin_angle) + point.z() * cos_angle);
}


/**
 Rotates the reference system about the y-axis with the right-hand rule.
 e.g. RotateAboutY(Point(1.0,0.0,0.0), pi/2) == Point(0.0,0.0,-1.0)
 */
template<typename T>
Point<T> RotateAboutY(
  const Point<T>& point,
  const T angle) noexcept
{
  T cos_angle = cos(angle);
  T sin_angle = sin(angle);

  return Point
  (
    point.x() * cos_angle + point.z() * sin_angle,
    point.y(),
    point.x() * (-sin_angle) + point.z() * cos_angle);
}


/**
 Rotates the reference system about the z-axis with the right-hand rule.
 e.g. RotateAboutZ(Point(0.0,1.0,0.0), pi/2) == Point(-1.0,0.0,0.0)
 */
template<typename T>
Point<T> RotateAboutZ(
  const Point<T>& point,
  const T angle) noexcept
{
  T cos_angle = cos(angle);
  T sin_angle = sin(angle);
  return Point
  (
    point.x() * cos_angle + point.y() * (-sin_angle),
    point.x() * sin_angle + point.y() * cos_angle,
    point.z());
}


/**
 Rotates the reference system with euler angles. Convention is ZYX with
 angles phi, theta and psi, respectively.
 */
template<typename T>
Point<T> Rotate(
  const Point<T>& point,
  const T theta_1,
  const T theta_2,
  const T theta_3) noexcept
{
  Point rotated_about_z = RotateAboutZ(point, theta_1);
  Point rotated_about_y = RotateAboutY(rotated_about_z, theta_2);
  Point rotated_about_x = RotateAboutX(rotated_about_y, theta_3);
  return rotated_about_x;
}


template<typename T>
T DotProduct(
  const Point<T>& point_a,
  const Point<T>& point_b) noexcept
{
  return point_a.x() * point_b.x() + point_a.y() * point_b.y() + point_a.z() *
    point_b.z();
}


template<typename T>
T Distance(
  const Point<T>& point_a,
  const Point<T>& point_b) noexcept
{
  Point point
  (
    point_a.x() - point_b.x(),
    point_a.y() - point_b.y(),
    point_a.z() - point_b.z());
  return point.norm();
}


template<typename T>
T Theta(
  const Point<T>& point_a,
  const Point<T>& point_b) noexcept
{
  // point_a is taken as the new centre of the reference system
  Point point
  (
    point_b.x() - point_a.x(),
    point_b.y() - point_a.y(),
    point_b.z() - point_a.z());
  return point.theta();
}


template<typename T>
T Phi(
  const Point<T>& point_a,
  const Point<T>& point_b) noexcept
{
  // point_a is taken as the new centre of the reference system
  Point point
  (
    point_b.x() - point_a.x(),
    point_b.y() - point_a.y(),
    point_b.z() - point_a.z());
  return point.phi();
}


template<typename T>
T AngleBetweenDirections(
  T theta_a,
  T phi_a,
  T theta_b,
  T phi_b) noexcept
{
  Point point_a(
    sin(theta_a) * cos(phi_a), sin(theta_a) * sin(phi_a), cos(theta_a));
  Point point_b(
    sin(theta_b) * cos(phi_b), sin(theta_b) * sin(phi_b), cos(theta_b));
  return acos(DotProduct(point_a, point_b));
}


template<typename T>
T AngleBetweenPoints(
  Point<T> point_a,
  Point<T> point_b) noexcept
{
  point_a.Normalize();
  point_b.Normalize();
  return acos(DotProduct(point_a, point_b));
}


/**
 Contructs a point from spherical coordinates, with (r, 0, 0) corresponding
 to the z-axis, and (r, pi/2, 0) corresponding to x-axis. Right-hand rule.
 */
template<typename T>
Point<T> PointSpherical(
  const T r,
  const T theta,
  const T phi) noexcept
{
  T x = r * cos(phi) * sin(theta);
  T y = r * sin(phi) * sin(theta);
  T z = r * cos(theta);
  return Point(x, y, z);
}


/**
 This returns the point on the line between `point_a` and `point_b` which
 has a distance of `distance` from `point_a`
 */
template<typename T>
Point<T> PointOnLine(
  const Point<T> point_a,
  const Point<T> point_b,
  const T distance) noexcept
{
  Point point_centered = Subtract(point_b, point_a);
  Point out_point_centered = PointSpherical(
    distance, point_centered.theta(), point_centered.phi());
  return Sum(out_point_centered, point_a);
}


/** Sums the coordinates of `point_a` and `point_b` */
template<typename T>
Point<T> Sum(
  const Point<T> point_a,
  const Point<T> point_b) noexcept
{
  return Point
  (
    point_a.x() + point_b.x(),
    point_a.y() + point_b.y(),
    point_a.z() + point_b.z());
}


/** Subtracts the coordinates of `point_a` from `point_b` (point_a-point_b) */
template<typename T>
Point<T> Subtract(
  const Point<T> point_a,
  const Point<T> point_b) noexcept
{
  return Point
  (
    point_a.x() - point_b.x(),
    point_a.y() - point_b.y(),
    point_a.z() - point_b.z());
}


/**
 Multiplies all coordinates by given constant. Has the effect of changing
 of changing the length of the vector.
 */
template<typename T>
Point<T> Multiply(
  const Point<T> point,
  const T constant) noexcept
{
  return Point
  (
    point.x() * constant,
    point.y() * constant,
    point.z() * constant);
}


/**
 Constructs a vector that is the the projection of the input `vector`
 on the plane (passing through the origin) identified by the vector
 normal to the plane `plane_normal_vector`.
 */
template<typename T>
Point<T> Projection(
  const Point<T>& vector,
  const Point<T>& plane_normal_vector) noexcept
{
  Point normalised_plane_normal_vector = Normalized(plane_normal_vector);
  return Subtract
  (
    vector,
    Multiply
    (
      normalised_plane_normal_vector,
      DotProduct(vector, normalised_plane_normal_vector)));
}


/**
 Returns a new point that is a normalized (norm == 1) version of `point`.
 */
template<typename T>
Point<T> Normalized(
  Point<T> point) noexcept
{
  point.Normalize();
  return point;
}


/**
 Returns whther or not an intersection between a plane and a line exists. 
 The line is identified by a point on a line, line_point, 
 and the direction of the line,
 line_direction (every point on the line can be expressed as
 p=d line_point+line_direction, with d any scalar).
 The plane is identified as any point on the plane, plane_point, and
 the normal to the plane, plane_normal (every point on the plane can be
 expressed as (p-plane_point, plane_normal)=0 where (x,y) is scalar product).
 */
template<typename T>
bool IntersectionPlaneLineExists(
  const Point<T>& line_point,
  const Point<T>& line_direction,
  const Point<T>& plane_point,
  const Point<T>& plane_normal) noexcept
{
  // TODO: find a way to avoid rewriting this calculation in the next function
  T numerator = DotProduct
  (
    Subtract(plane_point, line_point),
    plane_normal);
  T denominator = DotProduct(line_direction, plane_normal);

  bool zero_numerator = IsEqual(numerator, 0.0);
  bool zero_denominator = IsEqual(denominator, 0.0);

  // If denominator = 0 then the line and plane are parallel.
  // There will be two cases:
  // if numerator = 0 then the line is contained in the plane, that is, the line
  // intersects the plane at each point of the line (and thus an intersection
  // still exists).
  // Otherwise, the line and plane have no intersection.
  // So, the case of no intersection is when the denominator = 0 and the
  // numerator is not = 0.
  return ! (zero_denominator & ! zero_numerator);
}


/**
 Returns the intersection point between a plane and a line. The line is
 identified by a point on a line, line_point, and the direction of the line,
 line_direction (every point on the line can be expressed as
 p=d line_point+line_direction, with d any scalar).
 The plane is identified as any point on the plane, plane_point, and
 the normal to the plane, plane_normal (every point on the plane can be
 expressed as (p-plane_point, plane_normal)=0 where (x,y) is scalar product).
 If the line is contained in the plane, the function returns line_point.
 If there is no intersection, then the function returns a point of NANs.
 The user should first check whether an intersection exists using
 IntersectionPlaneLineExists.
 */
template<typename T>
Point<T> IntersectionPlaneLine(
  const Point<T>& line_point,
  const Point<T>& line_direction,
  const Point<T>& plane_point,
  const Point<T>& plane_normal) noexcept
{
  if (! IntersectionPlaneLineExists
    (
      line_point,
      line_direction,
      plane_point,
      plane_normal))
  {
    return Point<double>(NAN, NAN, NAN);
  }

  T d = DotProduct
    (
      Subtract(plane_point, line_point),
      plane_normal) /
    DotProduct(line_direction, plane_normal);

  // if line and plane are parallel, and line is contained in plane
  if (IsNan(d))
  {
    return line_point;
  }
  return Sum(Multiply(line_direction, d), line_point);
}


template<typename T>
bool IsApproximatelyEqual(
  const Point<T>& point_a,
  const Point<T>& point_b,
  const T precision = T(VERY_SMALL)) noexcept
{
  return
    mcl::IsApproximatelyEqual(point_a.x(), point_b.x(), precision) &&
    mcl::IsApproximatelyEqual(point_a.y(), point_b.y(), precision) &&
    mcl::IsApproximatelyEqual(point_a.z(), point_b.z(), precision);
}


template<typename T>
bool IsEqual(
  const Point<T>& point_a,
  const Point<T>& point_b)
{
  return
    point_a.x() == point_b.x() &&
    point_a.y() == point_b.y() &&
    point_a.z() == point_b.z();
}


template<typename T>
bool IsEqual(
  const Vector<Point<T>>& points_a,
  const Vector<Point<T>>& points_b) noexcept
{
  const Int num_points = (Int)points_a.size();
  if (num_points != (Int)points_b.size())
  {
    return false;
  }
  for (Int i = 0; i < num_points; ++i)
  {
    if (! IsEqual(points_a[i], points_b[i]))
    {
      return false;
    }
  }
  return true;
}
} // namespace mcl
