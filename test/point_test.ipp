/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "point.hpp"
#include <cmath>
#include <cassert>

namespace mcl
{

inline bool PointTest()
{
  const Real pi = 4*atan(1);
  using mcl::IsApproximatelyEqual;
  
  Point<Real> point_a(0.0,0.0,1.0);
  ASSERT(IsApproximatelyEqual(point_a,Point<Real>(0.0,0.0,1.0)));
  ASSERT(! IsApproximatelyEqual(point_a, Point<Real>(1.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(point_a.norm(), 1.0));
  ASSERT(IsApproximatelyEqual(point_a.theta(),0));
  
  Point<Real> point_b(0.0,0.0,0.0);
  ASSERT(IsApproximatelyEqual(Distance(point_a, point_b),1.0));
  ASSERT(IsApproximatelyEqual(point_b.norm(), 0.0));
  
  Point<Real> point_c(1.0,2.0,-1.0);
  ASSERT(IsApproximatelyEqual(Distance(point_a,point_c), Distance(point_c,point_a)));
  ASSERT(IsApproximatelyEqual(Distance(point_a,point_c), sqrt(4.0+1.0+4.0)));
  ASSERT(IsApproximatelyEqual(DotProduct(point_a,point_c), -1.0));
  ASSERT(IsApproximatelyEqual(point_c.norm(), sqrt(1.0+pow(2.0,2)+1.0)));
  
  Point<Real> point_d(1.0,0.0,1.0);
  ASSERT(IsApproximatelyEqual(point_d.phi(), 0.0));
  ASSERT(IsApproximatelyEqual(point_d.theta(), pi/4.0));
  ASSERT(IsApproximatelyEqual(point_d.norm(), sqrt(2.0)));
  ASSERT(IsApproximatelyEqual(Theta(point_b,point_d), pi/4.0));
  
  Point<Real> point_e(0.0,1.0,0.0);
  ASSERT(IsApproximatelyEqual(point_e.phi(), pi/2.0));
  ASSERT(IsApproximatelyEqual(point_e.theta(), pi/2.0));
  ASSERT(IsApproximatelyEqual(Theta(point_e,point_b), pi/2.0));
  ASSERT(IsApproximatelyEqual(Theta(point_b,point_e), pi/2.0));
  ASSERT(IsApproximatelyEqual(Phi(point_b,point_e), pi/2.0));
  ASSERT(IsApproximatelyEqual(Phi(point_e,point_b), -pi/2.0));
  
  ASSERT(IsApproximatelyEqual(AngleBetweenDirections(point_a.theta(), point_a.phi(),
                                        point_d.theta(), point_d.phi()),
                 pi/4.0));
  ASSERT(IsApproximatelyEqual(AngleBetweenDirections(point_d.theta(), point_d.phi(),
                                        point_a.theta(), point_a.phi()),
                 pi/4.0));
  
  ASSERT(IsApproximatelyEqual(RotateAboutX(point_e, 0.0), point_e));
  ASSERT(IsApproximatelyEqual(RotateAboutY(point_e, 0.0), point_e));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(point_e, 0.0), point_e));
  
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(1.0,0.0,0.0), pi/2.0), Point<Real>(1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(1.0,0.0,0.0), -pi/2.0), Point<Real>(1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(0.0,1.0,0.0), pi/2.0), Point<Real>(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(0.0,1.0,0.0), -pi/2.0), Point<Real>(0.0,0.0,-1.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(0.0,0.0,1.0), pi/2.0), Point<Real>(0.0,-1.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(0.0,0.0,1.0), -pi/2.0), Point<Real>(0.0,1.0,0.0)));
  
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(1.0,0.0,0.0), pi/2.0), Point<Real>(0.0,0.0,-1.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(1.0,0.0,0.0), -pi/2.0), Point<Real>(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(0.0,1.0,0.0), pi/2.0), Point<Real>(0.0,1.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(0.0,1.0,0.0), -pi/2.0), Point<Real>(0.0,1.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(0.0,0.0,1.0), pi/2.0), Point<Real>(1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(0.0,0.0,1.0), -pi/2.0), Point<Real>(-1.0,0.0,0.0)));
  
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(1.0,0.0,0.0), pi/2.0), Point<Real>(0.0,1.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(1.0,0.0,0.0), -pi/2.0), Point<Real>(0.0,-1.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(0.0,1.0,0.0), pi/2.0), Point<Real>(-1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(0.0,1.0,0.0), -pi/2.0), Point<Real>(1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(0.0,0.0,1.0), pi/2.0), Point<Real>(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(0.0,0.0,1.0), -pi/2.0), Point<Real>(0.0,0.0,1.0)));
  
  ASSERT(IsApproximatelyEqual(RotateAboutX(Point<Real>(1.5,0.0,0.0), pi/2.0), Point<Real>(1.5,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutY(Point<Real>(0.0,1.5,0.0), pi/2.0), Point<Real>(0.0,1.5,0.0)));
  ASSERT(IsApproximatelyEqual(RotateAboutZ(Point<Real>(0.0,0.0,1.5), pi/2.0), Point<Real>(0.0,0.0,1.5)));
  
  Point<Real> point_f = RotateAboutX(RotateAboutY(
                         RotateAboutZ(Point<Real>(1.5,-1.0,0.5),pi/5.0),pi/6.0),pi/7.0);
  AVOID_UNUSED_WARNING(point_f);
  
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(1.5,-1.0,0.5), pi/5.0, pi/6.0, pi/7.0), point_f));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,0.0,1.0), pi/2.0, 0.0, 0.0), Point<Real>(0.0,0.0,1.0)));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,0.0,1.0), 0.0, pi/2.0, 0.0), Point<Real>(1.0,0.0,0.0)));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,0.0,1.0), 0.0, 0.0, pi/2.0), Point<Real>(0.0,-1.0,0.0)));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(1.0,0.0,0.0), pi/2.0, 0.0, 0.0), Point<Real>(0.0,1.0,0.0)));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,1.0,0.0), 0.0, pi/2.0, 0.0), Point<Real>(0.0,1.0,0.0)));
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,1.0,0.0), -PI/2.0, -PI/2.0, 0.0), Point<Real>(0.0,0.0,1.0)));

  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,1.0,0.0), PI/2.0, PI/2.0, 0.0), Point<Real>(0.0,0.0,1.0)));
  
  // Proof that is extrinsic, not intrinsic
  ASSERT(!IsApproximatelyEqual(Rotate(Point<Real>(0.0,1.0,0.0), PI/2.0, PI/2.0, 0.0), Point<Real>(-1.0,0.0,0.0))); // Intrinsic rotation
  ASSERT(IsApproximatelyEqual(Rotate(Point<Real>(0.0,1.0,0.0), PI/2.0, PI/2.0, 0.0), Point<Real>(0.0,0.0,1.0))); // Extrinsic rotation

  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, PI/2.0, 0.0), Point<Real>(1.0, 0.0, 0.0)));
  ASSERT(IsApproximatelyEqual(PointSpherical(2.0, PI/2.0, 0.0), Point<Real>(2.0, 0.0, 0.0)));
  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, 0.0, 0.0), Point<Real>(0.0, 0.0, 1.0)));
  ASSERT(IsApproximatelyEqual(PointSpherical(2.0, 0.0, 0.0), Point<Real>(0.0, 0.0, 2.0)));
  ASSERT(IsApproximatelyEqual(PointSpherical(1.5, PI/2.0, PI/2.0), Point<Real>(0.0, 1.5, 0.0)));
  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, -PI/2.0, 0.0), Point<Real>(-1.0, 0.0, 0.0)));

  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, PI/4.0, 0.0), Point<Real>(1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))));
  
  
  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, PI/4.0, PI/2.0), Point<Real>(0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0))));
  
  
  ASSERT(IsApproximatelyEqual(Sum(Point<Real>(1.0,-1.0,2.0), Point<Real>(0.0,0.0,0.0)), Point<Real>(1.0,-1.0,2.0)));
  ASSERT(IsApproximatelyEqual(Sum(Point<Real>(0.0,0.0,0.0), Point<Real>(1.0,-1.0,2.0)), Point<Real>(1.0,-1.0,2.0)));
  ASSERT(IsApproximatelyEqual(Sum(Point<Real>(0.5,1.0,-0.0), Point<Real>(1.0,-1.0,2.0)), Point<Real>(1.5,0.0,2.0)));
  
  ASSERT(IsApproximatelyEqual(Subtract(Point<Real>(1.0,-1.0,2.0), Point<Real>(0.0,0.0,0.0)), Point<Real>(1.0,-1.0,2.0)));
  ASSERT(IsApproximatelyEqual(Subtract(Point<Real>(0.0,0.0,0.0), Point<Real>(1.0,-1.0,2.0)), Point<Real>(-1.0,1.0,-2.0)));
  ASSERT(IsApproximatelyEqual(Subtract(Point<Real>(0.5,1.0,-0.0), Point<Real>(1.0,-1.0,2.0)), Point<Real>(-0.5,2.0,-2.0)));
  
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,0,0), Point<Real>(1,0,0), 0.3), Point<Real>(0.3,0,0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(1,0,0), Point<Real>(0,0,0), 0.3), Point<Real>(0.7,0,0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,0,0), Point<Real>(0,1,0), 0.3), Point<Real>(0,0.3,0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,1,0), Point<Real>(0,0,0), 0.3), Point<Real>(0,0.7,0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,0,0), Point<Real>(2,2,0), 1.0), Point<Real>(1.0/sqrt(2.0),1.0/sqrt(2.0),0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,0,0), Point<Real>(0,2,2), 1.0), Point<Real>(0,1.0/sqrt(2.0),1.0/sqrt(2.0))));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(0,2,2), Point<Real>(0,0,0), 1.0), Point<Real>(0,2-1.0/sqrt(2.0),2-1.0/sqrt(2.0))));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(1,1,0), Point<Real>(2,2,0), 1.0), Point<Real>(1.0+1.0/sqrt(2.0),1.0+1.0/sqrt(2.0),0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(1,1,0), Point<Real>(1,3,0), 1.0), Point<Real>(1,2.0,0)));
  ASSERT(IsApproximatelyEqual(PointOnLine(Point<Real>(1,1,0), Point<Real>(3,1,0), 1.0), Point<Real>(2.0,1,0)));
  
  ASSERT(IsApproximatelyEqual(PointSpherical(1.0, PI/4.0, PI/2.0), Point<Real>(0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0))));
  
  Point<Real> point_g(1.0, 2.0, 1.5);
  AVOID_UNUSED_WARNING(point_g);
  Point<Real> point_h = Normalized(point_g);
  AVOID_UNUSED_WARNING(point_h);
  ASSERT(IsApproximatelyEqual(point_g, Point<Real>(1.0, 2.0, 1.5)));
  ASSERT(IsApproximatelyEqual(point_h.norm(), 1.0));
  ASSERT(IsApproximatelyEqual(point_h, Point<Real>(0.371390676354104,
                                       0.742781352708207,
                                       0.557086014531156), VERY_SMALL));
  
  // Test projections
  Point<Real> point_i(2.0, 0.0, 0.0);
  ASSERT(IsApproximatelyEqual(Projection(point_i, Point<Real>(1.0, 0.0, 0.0)), Point<Real>(0.0, 0.0, 0.0)));
  ASSERT(IsApproximatelyEqual(Projection(point_i, Point<Real>(0.0, 1.0, 0.0)), point_i));
  ASSERT(IsApproximatelyEqual(Projection(point_i, Point<Real>(0.0, 2.0, 0.0)), point_i));
  ASSERT(IsApproximatelyEqual(Projection(point_i, Point<Real>(0.0, 0.0, 1.0)), point_i));
  ASSERT(IsApproximatelyEqual(Projection(point_i, Point<Real>(0.0, 0.0, 2.0)), point_i));
  
  Point<Real> point_l(1.5, 2.0, 3.0);
  Point<Real> point_l_cmp =
          PointSpherical(point_l.norm()*cos(PI/2.0-point_l.theta()), PI/2.0, point_l.phi());
  AVOID_UNUSED_WARNING(point_l_cmp);
  ASSERT(IsApproximatelyEqual(Projection(point_l, Point<Real>(0,0,1)), point_l_cmp));
  
  // Test IntersectionPlaneLine
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(1,0,0),
                                       Point<Real>(1,0,0), Point<Real>(1,0,0)),
                 Point<Real>(1,0,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(2,0,0),
                                       Point<Real>(1,0,0), Point<Real>(1,0,0)),
                 Point<Real>(1,0,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(1,0,0),
                                       Point<Real>(2,0,0), Point<Real>(1,0,0)),
                 Point<Real>(2,0,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(1,0,0),
                                       Point<Real>(2,0,0), Point<Real>(1,0,0)),
                 Point<Real>(2,0,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,1,0), Point<Real>(1,0,0),
                                       Point<Real>(2,0,0), Point<Real>(1,0,0)),
                 Point<Real>(2,1,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(1,0,0),
                                       Point<Real>(2,0,0), Point<Real>(1,0,0)),
                 Point<Real>(2,0,0)));
  
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(1,1,-1), Point<Real>(0,0,1),
                                       Point<Real>(0,0,1), Point<Real>(0,0,1)),
                 Point<Real>(1,1,1)));
  
  // Case line parallel to plane: returning line_point
  ASSERT(IsApproximatelyEqual(IntersectionPlaneLine(Point<Real>(0,0,0), Point<Real>(0,1,0),
                                       Point<Real>(2,0,0), Point<Real>(0,0,1)),
                 Point<Real>(0,0,0)));
  ASSERT(IntersectionPlaneLineExists(Point<Real>(0,0,0), Point<Real>(0,1,0),
                                     Point<Real>(2,0,0), Point<Real>(0,0,1)));
  
  ASSERT(! IntersectionPlaneLineExists(Point<Real>(0,0,1), Point<Real>(0,1,0),
                                       Point<Real>(2,0,0), Point<Real>(0,0,1)));
  
  
  
  return true;
}

} // namespace sal
