/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mclintrinsics.hpp"

namespace mcl
{

inline bool MclIntrinsicsTest()
{

  Vector<Real> vector_z = ColonOperator(0.0, 2.0, 4.0);
  
  // Testing summation
  Vector<Real> vector_zb = AddScalar(vector_z, 1.5);
  ASSERT(vector_zb.size() == 3);
  Vector<Real> vector_zb_cmp(3);
  vector_zb_cmp[0] = 1.5;
  vector_zb_cmp[1] = 3.5;
  vector_zb_cmp[2] = 5.5;
  ASSERT(IsEqual(vector_zb, vector_zb_cmp));

  Vector<Real> vector_zc = AddScalar(vector_z, -1.0);
  ASSERT(vector_zc.size() == 3);
  Vector<Real> vector_zc_cmp(3);
  vector_zc_cmp[0] = -1.0;
  vector_zc_cmp[1] = 1.0;
  vector_zc_cmp[2] = 3.0;
  ASSERT(IsEqual(vector_zc, vector_zc_cmp));



  Vector<Real> vector_ba(4);
  vector_ba[0] = -1.2;
  vector_ba[1] = 2.3;
  vector_ba[2] = 3.4;
  vector_ba[3] = -5.0;
  
  Vector<Real> vector_bc(4);
  vector_bc[0] = 0.0;
  vector_bc[1] = 1.0;
  vector_bc[2] = -0.5;
  vector_bc[3] = -0.0;
  Vector<Real> vector_bb_result(4);
  MultiplyAdd(vector_ba, 0.5, vector_bc, vector_bb_result);
  Vector<Real> vector_bb_result_cmp(4);
  vector_bb_result_cmp[0] = -1.2*0.5+0.0;
  vector_bb_result_cmp[1] = 2.3*0.5+1.0;
  vector_bb_result_cmp[2] = 3.4*0.5-0.5;
  vector_bb_result_cmp[3] = -5.0*0.5+0.0;
  ASSERT(IsEqual(vector_bb_result, vector_bb_result_cmp));




  Vector<Real> vector_q(2);
  vector_q[0] = -1.2;
  vector_q[1] = 4.5;
  Vector<Real> vector_o = Zeros<Real>(3);
  vector_o[0] = 1.0;
  vector_o[1] = 2.5;
  vector_o[2] = 4.2;
  Vector<Real> vector_m = LinSpace(0.0, -2.0, 3);

  Vector<Vector<Real>> vectors(3);
  vectors[0] = vector_q; // -1.2, 4.5
  vectors[1] = vector_o; // 1.0, 2.5, 4.2
  vectors[2] = vector_m; // 0.0, -1.0, -2.0

  Vector<Real> add_vectors = AddVectors(vectors);
  ASSERT(add_vectors.size() == 3);
  Vector<Real> add_vectors_cmp(3);
  add_vectors_cmp[0] = -1.2+1.0;
  add_vectors_cmp[1] = 4.5+2.5-1.0;
  add_vectors_cmp[2] = 4.2-2.0;
  ASSERT(IsEqual(add_vectors_cmp, add_vectors));





  Vector<Complex<Real>> vector_a(3);
  vector_a[0] = Complex<Real>(1.0, 0.0);
  vector_a[1] = Complex<Real>(0.0, 1.0);
  vector_a[2] = Complex<Real>(1.0, 0.5);

  Vector<Complex<Real>> vector_b(3);
  vector_b[0] = Complex<Real>(0.5, 1.0);
  vector_b[1] = Complex<Real>(0.2, -1.0);
  vector_b[2] = Complex<Real>(-1.0, -0.5);


  Vector<Complex<Real>> mult_vector = Multiply(vector_a, vector_b);
  Vector<Complex<Real>> mult_vector_cmp(3);
  mult_vector_cmp[0] = Complex<Real>(0.5, 1.0);
  mult_vector_cmp[1] = Complex<Real>(1.0, 0.2);
  mult_vector_cmp[2] = Complex<Real>(-0.75, -1.0);

  ASSERT(IsEqual(mult_vector, mult_vector_cmp));

  Vector<Complex<Real>> add_vector = Add(vector_a, vector_b);
  Vector<Complex<Real>> add_vector_cmp(3);
  add_vector_cmp[0] = Complex<Real>(1.5, 1.0);
  add_vector_cmp[1] = Complex<Real>(0.2, 0.0);
  add_vector_cmp[2] = Complex<Real>(0.0, 0.0);

  ASSERT(IsEqual(add_vector, add_vector_cmp));

  Vector<Complex<Real>> sub_vector = Subtract(vector_a, vector_b);
  Vector<Complex<Real>> sub_vector_cmp(3);
  sub_vector_cmp[0] = Complex<Real>(0.5, -1.0);
  sub_vector_cmp[1] = Complex<Real>(-0.2, +2.0);
  sub_vector_cmp[2] = Complex<Real>(2.0, 1.0);

  ASSERT(IsEqual(sub_vector, sub_vector_cmp));
  
  return false;
}

}

