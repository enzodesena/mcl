/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "pointwiseop.h"
#include "comparisonop.h"
#include "vectorop.h"
#include "mcltypes.h"
#include <vector>

namespace mcl {
  
bool PointWiseOpTest() {
  
  Vector<Complex> vector_a(3);
  vector_a[0] = Complex(1.0, 0.0);
  vector_a[1] = Complex(0.0, 1.0);
  vector_a[2] = Complex(1.0, 0.5);
  
  
  ASSERT(! IsReal(vector_a));
  
  Vector<Real> abs_vector_a = Abs(vector_a);
  Vector<Real> abs_vector_a_cmp(3);
  abs_vector_a_cmp[0] = 1.0;
  abs_vector_a_cmp[1] = 1.0;
  abs_vector_a_cmp[2] = 1.118033988749895;
  
  ASSERT(IsEqual(abs_vector_a, abs_vector_a_cmp));
  
  
  Vector<Complex> vector_b(3);
  vector_b[0] = Complex(0.5, 1.0);
  vector_b[1] = Complex(0.2, -1.0);
  vector_b[2] = Complex(-1.0, -0.5);
  
  
  Vector<Complex> mult_vector = Multiply(vector_a, vector_b);
  Vector<Complex> mult_vector_cmp(3);
  mult_vector_cmp[0] = Complex(0.5, 1.0);
  mult_vector_cmp[1] = Complex(1.0, 0.2);
  mult_vector_cmp[2] = Complex(-0.75, -1.0);
  
  ASSERT(IsEqual(mult_vector, mult_vector_cmp));
  
  Vector<Complex> add_vector = Add(vector_a, vector_b);
  Vector<Complex> add_vector_cmp(3);
  add_vector_cmp[0] = Complex(1.5, 1.0);
  add_vector_cmp[1] = Complex(0.2, 0.0);
  add_vector_cmp[2] = Complex(0.0, 0.0);
  
  ASSERT(IsEqual(add_vector, add_vector_cmp));
  
  Vector<Complex> sub_vector = Subtract(vector_a, vector_b);
  Vector<Complex> sub_vector_cmp(3);
  sub_vector_cmp[0] = Complex(0.5, -1.0);
  sub_vector_cmp[1] = Complex(-0.2, +2.0);
  sub_vector_cmp[2] = Complex(2.0, 1.0);
  
  ASSERT(IsEqual(sub_vector, sub_vector_cmp));

  
  
  Vector<Real> vector_c(3);
  vector_c[0] = -0.3;
  vector_c[1] = 0.3;
  vector_c[2] = 2.4;
  
  Vector<Real> vector_c_inv = Inverse(vector_c);
  Vector<Real> vector_c_inv_cmp(3);
  vector_c_inv_cmp[0] = 1.0/(-0.3);
  vector_c_inv_cmp[1] = 1.0/(0.3);
  vector_c_inv_cmp[2] = 1.0/(2.4);
  ASSERT(IsEqual(vector_c_inv, vector_c_inv_cmp));
  
  
  Vector<Complex> vector_cc = ComplexVector(vector_c);
  Vector<Complex> vector_cc_cmp(3);
  vector_cc_cmp[0] = Complex(-0.3, 0.0);
  vector_cc_cmp[1] = Complex(0.3, 0.0);
  vector_cc_cmp[2] = Complex(2.4, 0.0);
  
  ASSERT(IsReal(vector_cc));
  
  ASSERT(IsEqual(vector_cc, vector_cc_cmp));

  
  
  Vector<Real> vector_e(4);
  vector_e[0] = -0.3;
  vector_e[1] = 30.3;
  vector_e[2] = 2.4;
  vector_e[3] = 12.4;
  
  
  // Testing Pow
  Vector<Real> pow_vector_e = Pow(vector_e, 3.0);
  
  Vector<Real> pow_vector_e_cmp(4);
  pow_vector_e_cmp[0] = pow(vector_e[0], 3.0);
  pow_vector_e_cmp[1] = pow(vector_e[1], 3.0);
  pow_vector_e_cmp[2] = pow(vector_e[2], 3.0);
  pow_vector_e_cmp[3] = pow(vector_e[3], 3.0);
  
  ASSERT(IsEqual(pow_vector_e, pow_vector_e_cmp));
  
  
  
  Vector<Real> abs_vector_e = Abs(vector_e);
  Vector<Real> abs_vector_e_cmp(4);
  abs_vector_e_cmp[0] = 0.3;
  abs_vector_e_cmp[1] = 30.3;
  abs_vector_e_cmp[2] = 2.4;
  abs_vector_e_cmp[3] = 12.4;
  
  ASSERT(IsEqual(abs_vector_e, abs_vector_e_cmp));
  
  Vector<Real> abs_pow_vector_e = Pow(abs_vector_e, -1.0/3.0);
  
  Vector<Real> abs_pow_vector_e_cmp(4);
  abs_pow_vector_e_cmp[0] = pow(abs_vector_e[0], -1.0/3.0);
  abs_pow_vector_e_cmp[1] = pow(abs_vector_e[1], -1.0/3.0);
  abs_pow_vector_e_cmp[2] = pow(abs_vector_e[2], -1.0/3.0);
  abs_pow_vector_e_cmp[3] = pow(abs_vector_e[3], -1.0/3.0);
  
  ASSERT(IsEqual(abs_pow_vector_e, abs_pow_vector_e_cmp));
  
  Vector<Real> hw_vector_e = HalfWave(vector_e);
  Vector<Real> hw_vector_e_cmp(4);
  hw_vector_e_cmp[0] = 0.0;
  hw_vector_e_cmp[1] = 30.3;
  hw_vector_e_cmp[2] = 2.4;
  hw_vector_e_cmp[3] = 12.4;
  
  ASSERT(IsEqual(hw_vector_e, hw_vector_e_cmp));
  
  
  
  Vector<double, 3> vector_o(Zeros<double,3>());
  vector_o[0] = 1.0;
  vector_o[1] = 2.5;
  vector_o[2] = 4.2;
  
  Vector<Real> log_vector_o = Log(vector_o);
  Vector<Real> log_vector_o_cmp(3);
  log_vector_o_cmp[0] = 0.0;
  log_vector_o_cmp[1] = 0.916290731874155;
  log_vector_o_cmp[2] = 1.435084525289323;
  ASSERT(IsEqual(log_vector_o_cmp, log_vector_o));
  
  Vector<Real> log10_vector_o = Log10(vector_o);
  Vector<Real> log10_vector_o_cmp(3);
  log10_vector_o_cmp[0] = 0.0;
  log10_vector_o_cmp[1] = 0.397940008672038;
  log10_vector_o_cmp[2] = 0.623249290397900;
  ASSERT(IsEqual(log10_vector_o_cmp, log10_vector_o));
  
  Vector<Complex> exp_vector_a = Exp(vector_a);
  Vector<Complex> exp_vector_a_cmp(3);
  exp_vector_a_cmp[0] = Complex(2.718281828459046, 0.0);
  exp_vector_a_cmp[1] = Complex(0.5403023058681398, 0.8414709848078965);
  exp_vector_a_cmp[2] = Complex(2.385516730959136, 1.303213729686996);
  ASSERT(IsEqual(exp_vector_a_cmp, exp_vector_a));

  Vector<Real> colonop_a = ColonOperator<Real>(2, 4);
  Vector<Real> colonop_a_cmp = Zeros<Real>(3);
  colonop_a_cmp[0] = 2.0;
  colonop_a_cmp[1] = 3.0;
  colonop_a_cmp[2] = 4.0;
  
  Vector<Real> cosvector = Cos(colonop_a_cmp);
  Vector<Real> cosvector_cmp(3);
  cosvector_cmp[0] = cos(2.0);
  cosvector_cmp[1] = cos(3.0);
  cosvector_cmp[2] = cos(4.0);
  ASSERT(IsEqual(cosvector_cmp, cosvector));

  Vector<Real> sinvector = Sin(colonop_a_cmp);
  Vector<Real> sinvector_cmp(3);
  sinvector_cmp[0] = sin(2.0);
  sinvector_cmp[1] = sin(3.0);
  sinvector_cmp[2] = sin(4.0);
  ASSERT(IsEqual(sinvector_cmp, sinvector));
  
  // Testing Divide
//  Vector<Real> vector_o = Zeros<Real>(3);
//  vector_o[0] = 1.0;
//  vector_o[1] = 2.5;
//  vector_o[2] = 4.2;
  Vector<Real> vector_p(3);
  vector_p[0] = -1.4;
  vector_p[1] = 2.3;
  vector_p[2] = 4.2;
  Vector<Real> division_o_p = Divide(vector_o, vector_p);
  Vector<Real> division_o_p_cmp(3);
  division_o_p_cmp[0] = 1.0/(-1.4);
  division_o_p_cmp[1] = 2.5/2.3;
  division_o_p_cmp[2] = 4.2/4.2;
  ASSERT(IsEqual(division_o_p_cmp, division_o_p));
  
  return true;
}
} // namespace mcl
