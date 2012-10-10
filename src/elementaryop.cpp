/*
 mathfunctions.h
 Matlab Cpp Library (MCL)
 
 Authors: Enzo De Sena, enzodesena@me.com
 
 Last committed:     $Revision: 95 $
 Last changed date:  $Date: 2012-06-07 20:07:36 +0100 (Thu, 07 Jun 2012) $
 */

#include "elementaryop.h"
#include "equalityop.h"

namespace mcl {



Int Sign(const Real scalar) {
  if (IsEqual(scalar, 0.0)) { return 0; }
  else if (scalar > 0.0) { return 1; }
  else { return -1; }
}

Int Fix(const Real scalar) {
  if (scalar >= 0.0) 
    return floor(scalar);
  else 
    return ceil(scalar);
}


Real Rem(const Real& x, const Real& y) {
  if (IsEqual(y, 0)) { return NAN; }
  if (IsEqual(x, y)) { return 0.0; }
  Int n = Fix(x/y);
  return x - ((Real) n)*y; 
}

Real Mod(const Real& x, const Real& y) {
  if (IsEqual(y, 0)) { return x; }
  if (IsEqual(x, y)) { return 0.0; }
  Real rem(Rem(x, y));
  Int signum(Sign(x/y));
  if (signum == 1 || signum == 0) { return rem; }
  else { return (rem + y); }
}

  


Int Mod(const Int& x, const Int& y) {
  if (y == 0) { return x; }
  if (x == y || x == -y) { return 0; }
  Int n = Fix(((Real) x)/((Real) y));
  Int rem(x - n*y);
  Int signum(Sign(((Real) x)/((Real) y)));
  if (signum == 1 || signum == 0) { return rem; }
  else { return (rem + y); }
}

  
Complex Conj(Complex scalar) {
  return Complex(scalar.real(), -scalar.imag());
}

Real RealPart(Complex scalar) {
  return scalar.real();
}
  
UInt NextPow2(Real input) { return ceil(log2(fabs(input))); }

  
double StringToDouble(const std::string& s) {
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
}
  
Real Abs(Real input) {
  return fabs(input);
}
  
} // namespace mcl