/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "vectorop.h"
#include "vector.h"
#include "comparisonop.h"
//#include "mcltypes.h"

namespace mcl {


  
bool VectorOpTest() {

  Vector<double, 3> myvector_a;
  myvector_a[0] = 0.1;
  myvector_a[1] = 0.2;
  myvector_a[2] = 0.3;
  ASSERT(myvector_a[0] == 0.1);
  ASSERT(myvector_a[1] == 0.2);
  ASSERT(myvector_a[2] == 0.3);
  ASSERT(myvector_a[2] == 0.3);
  myvector_a[2] = -0.3;
  ASSERT(myvector_a[2] == -0.3);
  
  Vector<double, kDynamicLength> myvector_b;
  myvector_b.PushBack(0.1);
  myvector_b.PushBack(0.2);
  myvector_b.PushBack(0.3);
  ASSERT(myvector_b[0] == 0.1);
  ASSERT(myvector_b[1] == 0.2);
  ASSERT(myvector_b[2] == 0.3);
  
  myvector_b[2] = -0.3;
  ASSERT(myvector_b[2] == -0.3);

  Vector<Complex<double>, 3> vector_a;
  vector_a[0] = Complex<double>(1.0, 0.0);
  vector_a[1] = Complex<double>(0.0, 1.0);
  vector_a[2] = Complex<double>(1.0, 0.5);


  ASSERT(Length(vector_a) == 3);

  Vector<Complex<double>, 3> flip_vector_a = Flip(vector_a);
  Vector<Complex<double>, 3> flip_vector_a_cmp(3);
  flip_vector_a_cmp[0] = Complex<double>(1.0, 0.5);
  flip_vector_a_cmp[1] = Complex<double>(0.0, 1.0);
  flip_vector_a_cmp[2] = Complex<double>(1.0, 0.0);
  ASSERT(IsEqual(flip_vector_a_cmp, flip_vector_a));

  Vector<Real> vector_bb(2);
  vector_bb[0] = 1.0;
  vector_bb[1] = -1.0;
  Vector<Real> flip_vector_bb = Flip(vector_bb);
  Vector<Real> flip_vector_bb_cmp(2);
  flip_vector_bb_cmp[0] = -1.0;
  flip_vector_bb_cmp[1] = 1.0;
  ASSERT(IsEqual(flip_vector_bb, flip_vector_bb_cmp));
  ASSERT(IsEqual(UnaryVector(2.2), Flip(UnaryVector(2.2))));

  Vector<Real> vector_l = LinSpace(0.0, 2.0, 3);
  ASSERT(vector_l.length() == 3);
  Vector<Real> vector_l_cmp(3);
  vector_l_cmp[0] = 0.0;
  vector_l_cmp[1] = 1.0;
  vector_l_cmp[2] = 2.0;
  ASSERT(IsEqual(vector_l_cmp, vector_l));

  Vector<Real> vector_m = LinSpace(0.0, -2.0, 3);
  ASSERT(vector_m.length() == 3);
  Vector<Real> vector_m_cmp(3);
  vector_m_cmp[0] = 0.0;
  vector_m_cmp[1] = -1.0;
  vector_m_cmp[2] = -2.0;
  ASSERT(IsEqual(vector_m_cmp, vector_m));

  Vector<Real> vector_n = LinSpace(1.0, 2.0, 4);
  ASSERT(vector_n.length() == 4);
  Vector<Real> vector_n_cmp(4);
  vector_n_cmp[0] = 1.0;
  vector_n_cmp[1] = 1.33333333333333333333;
  vector_n_cmp[2] = 1.66666666666665666667;
  vector_n_cmp[3] = 2.0;
  bool aa = IsVectorApproximatelyEqual(vector_n_cmp, vector_n, VERY_SMALL);
  ASSERT(aa);

  Vector<Real> vector_c(3);
  vector_c[0] = -0.3;
  vector_c[1] = 0.3;
  vector_c[2] = 2.4;

  Real vector_c_mean = Mean(vector_c);
  ASSERT(IsApproximatelyEqual(vector_c_mean, 0.8, VERY_SMALL));

  ASSERT(IsApproximatelyEqual(Sum(vector_c), -0.3+0.3+2.4, VERY_SMALL));

  Vector<Real, 5> pad_vector_c;
  ZeroPad(vector_c, pad_vector_c);
  Vector<Real> pad_vector_c_cmp(5);
  pad_vector_c_cmp[0] = -0.3;
  pad_vector_c_cmp[1] = 0.3;
  pad_vector_c_cmp[2] = 2.4;
  pad_vector_c_cmp[3] = 0.0;
  pad_vector_c_cmp[4] = 0.0;
  ASSERT(Length(pad_vector_c) == 5);
  ASSERT(IsEqual(pad_vector_c, pad_vector_c_cmp));


  Vector<Real> circ_2_vector_n = CircShift(vector_n, 2);
  Vector<Real> circ_2_vector_n_cmp(4);
  circ_2_vector_n_cmp[0] = vector_n[2];
  circ_2_vector_n_cmp[1] = vector_n[3];
  circ_2_vector_n_cmp[2] = vector_n[0];
  circ_2_vector_n_cmp[3] = vector_n[1];
  ASSERT(IsEqual(circ_2_vector_n_cmp, circ_2_vector_n));

  ASSERT(IsEqual(CircShift(vector_n, -2), CircShift(vector_n, 2)));

  Vector<Real> vector_d(4);
  vector_d[0] = -0.3;
  vector_d[1] = 0.3;
  vector_d[2] = 2.4;
  vector_d[3] = -12.4;

  Vector<Real> flip_vector_d = Flip(vector_d);
  Vector<Real> flip_vector_d_cmp(4);
  flip_vector_d_cmp[0] = -12.4;
  flip_vector_d_cmp[1] = 2.4;
  flip_vector_d_cmp[2] = 0.3;
  flip_vector_d_cmp[3] = -0.3;
  ASSERT(IsEqual(flip_vector_d_cmp, flip_vector_d));

  Vector<Real> vector_e(4);
  vector_e[0] = -0.3;
  vector_e[1] = 30.3;
  vector_e[2] = 2.4;
  vector_e[3] = 12.4;

  Vector<Real> vector_f(4);
  vector_f[0] = 2.5;
  vector_f[1] = 1.3;
  vector_f[2] = -2.4;
  vector_f[3] = -1.0;

  ASSERT(Sum(vector_e) == -0.3+30.3+2.4+12.4);
  
  ASSERT(Downsample(vector_f, 2).length() == 2);
  ASSERT(Downsample(vector_f, 3).length() == 2);
  ASSERT(Downsample(vector_f, 4).length() == 1);
  ASSERT(Downsample(vector_f, 10).length() == 1);
  ASSERT(Downsample(pad_vector_c_cmp, 2).length() == 3);
  ASSERT(Downsample(pad_vector_c_cmp, 3).length() == 2);
  ASSERT(Downsample(pad_vector_c_cmp, 4).length() == 2);
  ASSERT(Downsample(pad_vector_c_cmp, 5).length() == 1);

  Vector<Real> cmp_downsample_f(2);
  cmp_downsample_f[0] = 2.5;
  cmp_downsample_f[1] = -2.4;
  ASSERT(IsEqual(cmp_downsample_f, Downsample(vector_f, 2)));

  Vector<Real> cmp_downsample_f_3(2);
  cmp_downsample_f_3[0] = 2.5;
  cmp_downsample_f_3[1] = -1.0;
  ASSERT(IsEqual(cmp_downsample_f_3, Downsample(vector_f, 3)));

  Vector<Real> vector_f_sub_0_2 = Subset(vector_f, 0, 2);
  ASSERT(vector_f_sub_0_2.length() == 3);
  Vector<Real> vector_f_sub_0_2_cmp(3);
  vector_f_sub_0_2_cmp[0] = 2.5;
  vector_f_sub_0_2_cmp[1] = 1.3;
  vector_f_sub_0_2_cmp[2] = -2.4;
  ASSERT(IsEqual(vector_f_sub_0_2_cmp, vector_f_sub_0_2));

  Vector<Real> vector_f_sub_1_2 = Subset(vector_f, 1, 2);
  ASSERT(vector_f_sub_1_2.length() == 2);
  Vector<Real> vector_f_sub_1_2_cmp(2);
  vector_f_sub_1_2_cmp[0] = 1.3;
  vector_f_sub_1_2_cmp[1] = -2.4;
  ASSERT(IsEqual(vector_f_sub_1_2_cmp, vector_f_sub_1_2));


  Vector<Real> vector_g(3);
  vector_g[0] = 2.5;
  vector_g[1] = 0.0;
  vector_g[2] = -2.4;


  Vector<Real> vector_h(3);
  vector_h[0] = -2.5;
  vector_h[1] = 1.3;
  vector_h[2] = 2.4;

  Vector<Real> vector_g_h = Concatenate(vector_g, vector_h);
  ASSERT(vector_g[0] == 2.5);
  ASSERT(vector_h[1] == 1.3);
  ASSERT(vector_g_h.length() == 6);
  ASSERT(vector_g_h[0] == 2.5);
  ASSERT(vector_g_h[1] == 0.0);
  ASSERT(vector_g_h[3] == -2.5);



  Vector<Real> unary_vector = UnaryVector((Real) 2.0);
  ASSERT(Length(unary_vector) == 1);
  ASSERT(unary_vector[0] == 2.0);



//  Vector<Real> ones_a = Ones(3);
//  ASSERT(Length(ones_a) == 3);
//  ASSERT(IsEqual(ones_a[0], 1.0));
//  ASSERT(IsEqual(ones_a[1], 1.0));
//  ASSERT(IsEqual(ones_a[2], 1.0));
//
//  Vector<Real> vector_o = Zeros<Real>(3);
//  vector_o[0] = 1.0;
//  vector_o[1] = 2.5;
//  vector_o[2] = 4.2;
//
//
//  ASSERT(IsEqual(vector_o, CircShift(vector_o, 0)));
//  ASSERT(IsEqual(vector_o, CircShift(vector_o, 3)));
//  ASSERT(IsEqual(vector_o, CircShift(vector_o, -3)));
//  ASSERT(IsEqual(CircShift(vector_o, 1), CircShift(vector_o, 4)));
//  ASSERT(IsEqual(CircShift(vector_o, -1), CircShift(vector_o, -4)));
//  Vector<Real> circ_1_vector_o = CircShift(vector_o, 1);
//  Vector<Real> circ_1_vector_o_cmp(3);
//  circ_1_vector_o_cmp[0] = vector_o[2];
//  circ_1_vector_o_cmp[1] = vector_o[0];
//  circ_1_vector_o_cmp[2] = vector_o[1];
//  ASSERT(IsEqual(circ_1_vector_o_cmp, circ_1_vector_o));
//
//  Vector<Real> circ_m1_vector_o = CircShift(vector_o, -1);
//  Vector<Real> circ_m1_vector_o_cmp(3);
//  circ_m1_vector_o_cmp[0] = vector_o[1];
//  circ_m1_vector_o_cmp[1] = vector_o[2];
//  circ_m1_vector_o_cmp[2] = vector_o[0];
//  ASSERT(IsEqual(circ_m1_vector_o_cmp, circ_m1_vector_o));
//
//  Vector<Real> vector_p = LinSpace(0.0, 2.0, 3);
//  Vector<Real> vector_q(2);
//  vector_q[0] = -1.2;
//  vector_q[1] = 4.5;
//
//  ASSERT(IsEqual(Conv(vector_p, vector_q), Conv(vector_q,vector_p)));
//
//  Vector<Real> conv_p_q = Conv(vector_p, vector_q);
//  ASSERT(conv_p_q.length() == 4);
//  Vector<Real> conv_p_q_cmp(4);
//  conv_p_q_cmp[0] = 0.0;
//  conv_p_q_cmp[1] = -1.2;
//  conv_p_q_cmp[2] = 2.1;
//  conv_p_q_cmp[3] = 9;
//
//  ASSERT(IsEqual(conv_p_q_cmp, conv_p_q));
//
//  Vector<Real> vector_r(2);
//  vector_r[0] = 3.0;
//  vector_r[1] = 6.0;
//
//  Vector<Real> conv_q_r = Conv(vector_r, vector_q);
//  ASSERT(conv_q_r.length() == 3);
//  Vector<Real> conv_q_r_cmp(3);
//  conv_q_r_cmp[0] = -3.6;
//  conv_q_r_cmp[1] = 6.3;
//  conv_q_r_cmp[2] = 27.0;
//  ASSERT(IsEqual(conv_q_r, conv_q_r_cmp));
//
//
//  Vector<Real> vector_q_padded = Concatenate(Zeros<Real>(0),
//                                                  vector_q);
//  ASSERT(IsEqual(vector_q_padded, vector_q));
//
//
//
//
//  std::vector<Vector<Real> > vectors;
//  vectors.push_back(vector_q); // -1.2, 4.5
//  vectors.push_back(vector_o); // 1.0, 2.5, 4.2
//  vectors.push_back(vector_m); // 0.0, -1.0, -2.0
//
//  Vector<Real> add_vectors = AddVectors(vectors);
//  ASSERT(add_vectors.length() == 3);
//  Vector<Real> add_vectors_cmp(3);
//  add_vectors_cmp[0] = -1.2+1.0;
//  add_vectors_cmp[1] = 4.5+2.5-1.0;
//  add_vectors_cmp[2] = 4.2-2.0;
//  ASSERT(IsEqual(add_vectors_cmp, add_vectors));
//
//  Vector<Real> colonop_a = ColonOperator<Real>(2, 4);
//  Vector<Real> colonop_a_cmp = Zeros<Real>(3);
//  colonop_a_cmp[0] = 2.0;
//  colonop_a_cmp[1] = 3.0;
//  colonop_a_cmp[2] = 4.0;
//  ASSERT(IsEqual(colonop_a, colonop_a_cmp));
//
//  Vector<Real> colonop_b = ColonOperator<Real>(-1, 2);
//  Vector<Real> colonop_b_cmp = Zeros<Real>(4);
//  colonop_b_cmp[0] = -1.0;
//  colonop_b_cmp[1] = 0.0;
//  colonop_b_cmp[2] = 1.0;
//  colonop_b_cmp[3] = 2.0;
//  ASSERT(IsEqual(colonop_b, colonop_b_cmp));
//
//
//
//
//  std::vector<Int> colonop_c = ColonOperator<Int>(2, 4);
//  std::vector<Int> colonop_c_cmp = Zeros<Int>(3);
//  colonop_c_cmp[0] = 2;
//  colonop_c_cmp[1] = 3;
//  colonop_c_cmp[2] = 4;
//  ASSERT(IsEqual(colonop_c, colonop_c_cmp));
//
//
//  // Test Poly
//
//  std::vector<Complex> poly_a = Poly(ComplexVector(Ones(3)));
//  std::vector<Complex> poly_a_cmp = Zeros<Complex>(4);
//  poly_a_cmp[0] = Complex(1.0, 0.0);
//  poly_a_cmp[1] = Complex(-3.0, 0.0);
//  poly_a_cmp[2] = Complex(3.0, 0.0);
//  poly_a_cmp[3] = Complex(-1.0, 0.0);
//  ASSERT(IsEqual(poly_a_cmp, poly_a));
//
//  std::vector<Complex> roots_b(4);
//  roots_b[0] = Complex(1.3, 2.0);
//  roots_b[1] = Complex(-1.4, 0.0);
//  roots_b[2] = Complex(2.5, -1.0);
//  roots_b[3] = Complex(0.2, 0.0);
//  std::vector<Complex> poly_b = Poly(roots_b);
//  std::vector<Complex> poly_b_cmp = Zeros<Complex>(5);
//  poly_b_cmp[0] = Complex(1.0, 0.0);
//  poly_b_cmp[1] = Complex(-2.600000000000001, -1.000000000000000);
//  poly_b_cmp[2] = Complex(0.410000000000001, 2.500000000000000);
//  poly_b_cmp[3] = Complex(7.364000000000000, 4.720000000000000);
//  poly_b_cmp[4] = Complex(-1.470000000000000, -1.036000000000000);
//  ASSERT(IsEqual(poly_b, poly_b_cmp));
//
//  ASSERT(IsEqual(Poly(mcl::Zeros<Real>(0)), UnaryVector(Complex(1.0,0.0))));
//  ASSERT(IsEqual(Poly(mcl::Zeros<Complex>(0)), UnaryVector(Complex(1.0,0.0))));
//
//
//  // Testing CopyFrom
//
//  Vector<Real> vector_i(3);
//  vector_i[0] = 0.1;
//  vector_i[1] = -0.5;
//  vector_i[2] = 4.0;
//
//  Vector<Real> vector_i_restr_1 = Elements(vector_i, 1, 2);
//  ASSERT(vector_i_restr_1.length() == 2);
//  Vector<Real> vector_i_restr_1_cmp(2);
//  vector_i_restr_1_cmp[0] = -0.5;
//  vector_i_restr_1_cmp[1] = 4.0;
//  ASSERT(IsEqual(vector_i_restr_1, vector_i_restr_1_cmp));
//
//  Vector<Real> vector_v(3);
//  vector_v[0] = 0.1;
//  vector_v[1] = -0.5;
//  vector_v[2] = 4.0;
//
//  Vector<Real> vector_v_restr_1 = Elements(vector_v, 1, 2);
//  ASSERT(vector_v_restr_1.length() == 2);
//  Vector<Real> vector_v_restr_1_cmp(2);
//  vector_v_restr_1_cmp[0] = -0.5;
//  vector_v_restr_1_cmp[1] = 4.0;
//  ASSERT(IsEqual(vector_v_restr_1, vector_v_restr_1_cmp));
//
//  Vector<Real> vector_v_restr_2 = Elements(vector_v, 0, 0);
//  ASSERT(vector_v_restr_2.length() == 1);
//  Vector<Real> vector_v_restr_2_cmp(1);
//  vector_v_restr_2_cmp[0] = 0.1;
//  ASSERT(IsEqual(vector_v_restr_2, vector_v_restr_2_cmp));
//
//  // Testing GetSegment
//
//  //  Vector<Real> vector_g(3);
//  //  vector_g[0] = 2.5;
//  //  vector_g[1] = 0.0;
//  //  vector_g[2] = -2.4;
//
//  ASSERT(IsEqual(GetSegment(vector_g, 0, 1), UnaryVector((Real) 2.5)));
//  ASSERT(IsEqual(GetSegment(vector_g, 1, 1), UnaryVector((Real) 0.0)));
//  ASSERT(IsEqual(GetSegment(vector_g, 2, 1), UnaryVector((Real) -2.4)));
//
//  Vector<Real> vector_g_frame_0(2);
//  vector_g_frame_0[0] = 2.5;
//  vector_g_frame_0[1] = 0.0;
//  ASSERT(IsEqual(GetSegment(vector_g, 0, 2), vector_g_frame_0));
//  ASSERT(IsEqual(GetSegment(vector_g, 1, 2), UnaryVector((Real) -2.4)));
//  ASSERT(IsEqual<Real>(GetSegment(vector_g, 1, 2, true), BinaryVector<Real>(-2.4, 0.0)));
//  ASSERT(IsEqual(GetSegment(vector_g, 2, 2, false), Vector<Real>()));
//  ASSERT(IsEqual<Real>(GetSegment(vector_g, 2, 2, true), BinaryVector<Real>(0.0, 0.0)));
//
//  // Testing prod()
//  ASSERT(IsEqual(Prod(vector_g_frame_0), 0.0));
//  ASSERT(IsEqual(Prod(vector_i), -0.2));
//  ASSERT(IsEqual(Prod(vector_a), Complex(-0.5, 1.0)));
//
//  // Testing dot product
//  ASSERT(IsEqual(Dot(vector_g, vector_v), -9.350));
//
//
//  // Testing Std
//  ASSERT(IsEqual(Std(vector_v), 2.443358344574123));
//  ASSERT(IsEqual(Std(colonop_b_cmp), 1.290994448735806));
//
//  // Testing var
//  ASSERT(IsEqual(Var(vector_v), 5.96999999999999));
//
//  // Testing colon operator
//  Vector<Real> vector_z = ColonOperator(0.0, 2.0, 4.0);
//  Vector<Real> vector_z_cmp(3);
//  vector_z_cmp[0] = 0.0;
//  vector_z_cmp[1] = 2.0;
//  vector_z_cmp[2] = 4.0;
//  ASSERT(IsEqual(vector_z, vector_z_cmp));
//
//  Vector<Real> vector_za = ColonOperator(0.0, 3.0, 4.0);
//  Vector<Real> vector_za_cmp(2);
//  vector_za_cmp[0] = 0.0;
//  vector_za_cmp[1] = 3.0;
//  ASSERT(IsEqual(vector_za, vector_za_cmp));
//
//  Vector<Real> vector_aa = ColonOperator(-3.5, 3.0, 3.0);
//  Vector<Real> vector_aa_cmp(3);
//  vector_aa_cmp[0] = -3.5;
//  vector_aa_cmp[1] = -0.5;
//  vector_aa_cmp[2] = 2.5;
//  ASSERT(IsEqual(vector_aa, vector_aa_cmp));
//
//  Vector<Real> vector_ab = mcl::ColonOperator(-0.001, 0.00025, 0.001);
//  Vector<Real> vector_ab_cmp(9);
//  vector_ab_cmp[0] = -0.001;
//  vector_ab_cmp[1] = -0.00075;
//  vector_ab_cmp[2] = -0.0005;
//  vector_ab_cmp[3] = -0.00025;
//  vector_ab_cmp[4] = 0;
//  vector_ab_cmp[5] = 0.00025;
//  vector_ab_cmp[6] = 0.0005;
//  vector_ab_cmp[7] = 0.00075;
//  vector_ab_cmp[8] = 0.001;
//  ASSERT(IsEqual(vector_ab, vector_ab_cmp));
//
//  // Testing summation
//  Vector<Real> vector_zb = Add(vector_z, (Real) 1.5);
//  ASSERT(vector_zb.length() == 3);
//  Vector<Real> vector_zb_cmp(3);
//  vector_zb_cmp[0] = 1.5;
//  vector_zb_cmp[1] = 3.5;
//  vector_zb_cmp[2] = 5.5;
//  ASSERT(IsEqual(vector_zb, vector_zb_cmp));
//
//  Vector<Real> vector_zc = Add(vector_z, (Real) -1.0);
//  ASSERT(vector_zc.length() == 3);
//  Vector<Real> vector_zc_cmp(3);
//  vector_zc_cmp[0] = -1.0;
//  vector_zc_cmp[1] = 1.0;
//  vector_zc_cmp[2] = 3.0;
//  ASSERT(IsEqual(vector_zc, vector_zc_cmp));
//
//  // Testing UnaryVector
//  Vector<Real> vector_zd = UnaryVector<Real>(-1.0);
//  ASSERT(vector_zd.length() == 1);
//  ASSERT(IsEqual(vector_zd[0], -1.0));
//
//  std::vector<Int> vector_ze = UnaryVector<Int>(-2);
//  ASSERT(vector_ze.length() == 1);
//  ASSERT(vector_ze[0] == -2);
//
//  // Testing BinaryVector
//  Vector<Real> vector_zf = BinaryVector<Real>(-1.0, 2.0);
//  ASSERT(vector_zf.length() == 2);
//  ASSERT(IsEqual(vector_zf[0], -1.0));
//  ASSERT(IsEqual(vector_zf[1], 2.0));
//
//
//  // Testing Hanning window
//  Vector<Real> vector_hann_3 = Hann(3);
//  Vector<Real> vector_hann_3_cmp(3);
//  vector_hann_3_cmp[0] = 0.0;
//  vector_hann_3_cmp[1] = 1.0;
//  vector_hann_3_cmp[2] = 0.0;
//  ASSERT(IsEqual(vector_hann_3, vector_hann_3_cmp));
//
//  Vector<Real> vector_hann_4 = Hann(4);
//  Vector<Real> vector_hann_4_cmp(4);
//  vector_hann_4_cmp[0] = 0.0;
//  vector_hann_4_cmp[1] = 0.75;
//  vector_hann_4_cmp[2] = 0.75;
//  vector_hann_4_cmp[3] = 0.0;
//  ASSERT(IsEqual(vector_hann_4, vector_hann_4_cmp));
//
//  Vector<Real> vector_hann_5 = Hann(5);
//  Vector<Real> vector_hann_5_cmp(5);
//  vector_hann_5_cmp[0] = 0.0;
//  vector_hann_5_cmp[1] = 0.5;
//  vector_hann_5_cmp[2] = 1.0;
//  vector_hann_5_cmp[3] = 0.5;
//  vector_hann_5_cmp[4] = 0.0;
//  ASSERT(IsEqual(vector_hann_5, vector_hann_5_cmp));
//
//  ASSERT(IsNonNegative(vector_hann_5_cmp));
//
//  // Testing Hamming window
//  Vector<Real> vector_hamming_3 = Hamming(3);
//  Vector<Real> vector_hamming_3_cmp(3);
//  vector_hamming_3_cmp[0] = 0.08;
//  vector_hamming_3_cmp[1] = 1.0;
//  vector_hamming_3_cmp[2] = 0.08;
//  ASSERT(IsEqual(vector_hamming_3, vector_hamming_3_cmp));
//
//  Vector<Real> vector_hamming_4 = Hamming(4);
//  Vector<Real> vector_hamming_4_cmp(4);
//  vector_hamming_4_cmp[0] = 0.08;
//  vector_hamming_4_cmp[1] = 0.77;
//  vector_hamming_4_cmp[2] = 0.77;
//  vector_hamming_4_cmp[3] = 0.08;
//  ASSERT(IsEqual(vector_hamming_4, vector_hamming_4_cmp));
//
//  Vector<Real> vector_hamming_5 = Hann(5);
//  Vector<Real> vector_hamming_5_cmp(5);
//  vector_hamming_5_cmp[0] = 0.0;
//  vector_hamming_5_cmp[1] = 0.5;
//  vector_hamming_5_cmp[2] = 1.0;
//  vector_hamming_5_cmp[3] = 0.5;
//  vector_hamming_5_cmp[4] = 0.0;
//  ASSERT(IsEqual(vector_hamming_5, vector_hamming_5_cmp));
//
//  ASSERT(IsNonNegative(vector_hamming_5));
//
//  // Testing Tukey Window
//  ASSERT(IsEqual(TukeyWin(4, 1.0), vector_hann_4_cmp));
//  Vector<Real> vector_tukey_1 = TukeyWin(5,0.6);
//  Vector<Real> vector_tukey_1_cmp(5);
//  vector_tukey_1_cmp[0] = 0.0;
//  vector_tukey_1_cmp[1] = 0.933012701892219;
//  vector_tukey_1_cmp[2] = 1.0;
//  vector_tukey_1_cmp[3] = 0.933012701892219;
//  vector_tukey_1_cmp[4] = 0.0;
//  ASSERT(IsEqual(vector_tukey_1, vector_tukey_1_cmp));
//
//  Vector<Real>  aa = TukeyWin(1, 0.6);
//  ASSERT(TukeyWin(1, 0.6).length() == 1);
//  ASSERT(TukeyWin(0, 0.6).length() == 0);
//  ASSERT(IsEqual(TukeyWin(1, 0.6), UnaryVector<Real>(1.0)));
//
//  ASSERT(IsEqual(TukeyWin(6, 0.0), Ones(6)));
//  ASSERT(IsEqual(TukeyWin(6, -2.0), Ones(6)));
//
//  // Testing norm
//  Vector<Real> vector_ba(4);
//  vector_ba[0] = -1.2;
//  vector_ba[1] = 2.3;
//  vector_ba[2] = 3.4;
//  vector_ba[3] = -5.0;
//  ASSERT(IsEqual(Norm(vector_ba, 2.0), 6.579513659838392));
//  ASSERT(IsEqual(Norm(vector_ba, 1.0), 11.899999999999999));
//  ASSERT(IsEqual(Norm(vector_ba, 2.4), 6.056130782634900));
//
//
//  Vector<Real> vector_bc(4);
//  vector_bc[0] = 0.0;
//  vector_bc[1] = 1.0;
//  vector_bc[2] = -0.5;
//  vector_bc[3] = -0.0;
//  Vector<Real> vector_bb_result(4);
//  MultiplyAdd(vector_ba.data(), 0.5, vector_bc.data(), 4,
//              vector_bb_result.data());
//  Vector<Real> vector_bb_result_cmp(4);
//  vector_bb_result_cmp[0] = -1.2*0.5+0.0;
//  vector_bb_result_cmp[1] = 2.3*0.5+1.0;
//  vector_bb_result_cmp[2] = 3.4*0.5-0.5;
//  vector_bb_result_cmp[3] = -5.0*0.5+0.0;
//  ASSERT(IsEqual(vector_bb_result, vector_bb_result_cmp));
//
//  Vector<Real> weights_uniform_a = Multiply<Real>(Ones(4), 1.0/4.0);
//  ASSERT(IsEqual(Mean(vector_ba), Mean(vector_ba, weights_uniform_a)));
//  Vector<Real> weights_uniform_b = Multiply<Real>(Ones(4), 1.0);
//  ASSERT(IsEqual(Mean(vector_ba), Mean(vector_ba, weights_uniform_b)));
//  Vector<Real> weights_uniform_c = Zeros<Real>(4);
//  weights_uniform_c[0] = 0.5;
//  weights_uniform_c[2] = 0.5;
//  ASSERT(IsEqual(Mean(vector_ba, weights_uniform_c), 1.1));
//
//
//  ASSERT(! IsNonNegative(vector_ba));
//  Vector<Real> weights_ba_var = Zeros<Real>(4);
//  weights_ba_var[0] = 0.2;
//  weights_ba_var[1] = 0.3;
//  weights_ba_var[2] = 0.6;
//  weights_ba_var[3] = 0.5;
//  ASSERT(IsEqual(Var(vector_ba, weights_ba_var), 13.319335937500000));
//
//  // Testing covariance matrix
//
////  Vector<Real> vector_e(4);
////  vector_e[0] = -0.3;
////  vector_e[1] = 30.3;
////  vector_e[2] = 2.4;
////  vector_e[3] = 12.4;
////
////  Vector<Real> vector_f(4);
////  vector_f[0] = 2.5;
////  vector_f[1] = 1.3;
////  vector_f[2] = -2.4;
////  vector_f[3] = -1.0;
//  ASSERT(! IsEqual(Cov(vector_e, vector_f), Cov(vector_f, vector_e)));
//  Matrix<Real> cov_e_f = Cov(vector_e, vector_f);
//  ASSERT(IsEqual(cov_e_f.GetElement(0,0), Var(vector_e)));
//  ASSERT(IsEqual(cov_e_f.GetElement(0,0), 191.9800000000000));
//  ASSERT(IsEqual(cov_e_f.GetElement(1,1), 4.886666666666667));
//  ASSERT(IsEqual(cov_e_f.GetElement(1,1), Var(vector_f)));
//  ASSERT(IsEqual(cov_e_f.GetElement(0,0), Var(vector_e)));
//  ASSERT(IsEqual(cov_e_f.GetElement(0,1), cov_e_f.GetElement(1,0)));
//  ASSERT(IsEqual(cov_e_f.GetElement(0,1), 5.333333333333333));
//
//  // Cumsum
//  //  Vector<Real> vector_e(4);
//  //  vector_e[0] = -0.3;
//  //  vector_e[1] = 30.3;
//  //  vector_e[2] = 2.4;
//  //  vector_e[3] = 12.4;
//
//  Vector<Real> vector_cumsum_e = CumSum(vector_e);
//  Vector<Real> vector_cumsum_e_cmp(4);
//  vector_cumsum_e_cmp[0] = -0.300000000000000;
//  vector_cumsum_e_cmp[1] = 30.000000000000000;
//  vector_cumsum_e_cmp[2] = 32.399999999999999;
//  vector_cumsum_e_cmp[3] = 44.79999999999999;
//  ASSERT(IsEqual(vector_cumsum_e, vector_cumsum_e_cmp));
//
//  // Geomean
//  ASSERT(IsEqual(mcl::Geomean(vector_o), 2.189759569943945));
//
//
//  // Inteleaves
//  //  Vector<Real> vector_e(4);
//  //  vector_e[0] = -0.3;
//  //  vector_e[1] = 30.3;
//  //  vector_e[2] = 2.4;
//  //  vector_e[3] = 12.4;
//  //
//  //  Vector<Real> vector_f(4);
//  //  vector_f[0] = 2.5;
//  //  vector_f[1] = 1.3;
//  //  vector_f[2] = -2.4;
//  //  vector_f[3] = -1.0;
//  Vector<Real> interleaves_e_f_cmp;
//  interleaves_e_f_cmp.push_back(vector_e[0]);
//  interleaves_e_f_cmp.push_back(vector_f[0]);
//  interleaves_e_f_cmp.push_back(vector_e[1]);
//  interleaves_e_f_cmp.push_back(vector_f[1]);
//  interleaves_e_f_cmp.push_back(vector_e[2]);
//  interleaves_e_f_cmp.push_back(vector_f[2]);
//  interleaves_e_f_cmp.push_back(vector_e[3]);
//  interleaves_e_f_cmp.push_back(vector_f[3]);
//  ASSERT(IsEqual(Interleave(vector_e, vector_f), interleaves_e_f_cmp));
//
//
//  // Test Enframe
//  Vector<Real> vector = LinSpace(1, 10, 10);
//  Vector<Real> window = {0.2, 0.4};
//  std::vector<Vector<Real> > output_3 = Enframe(vector, window, 3);
//  ASSERT(output_3.length() == 3);
//  Vector<Real> output_3_0_cmp = {0.2, 0.8};
//  Vector<Real> output_3_1_cmp = {0.8, 2.0};
//  Vector<Real> output_3_2_cmp = {1.4, 3.2};
//  ASSERT(IsEqual(output_3[0], output_3_0_cmp));
//  ASSERT(IsEqual(output_3[1], output_3_1_cmp));
//  ASSERT(IsEqual(output_3[2], output_3_2_cmp));
//
//
//  std::vector<Vector<Real> > output_full = Enframe(Ones(4), window, 2);
//  ASSERT(output_full.length()==2);
//  Vector<Real> output_full_0_cmp = {0.2, 0.4};
//  Vector<Real> output_full_1_cmp = {0.2, 0.4};
//  ASSERT(IsEqual(output_full[0], output_full_0_cmp));
//  ASSERT(IsEqual(output_full[1], output_full_1_cmp));
//
//
//  Vector<Real> window_b;
//  window_b.push_back(0.2);
//  window_b.push_back(0.4);
//  window_b.push_back(-0.3);
//  std::vector<Vector<Real> > output_1 = Enframe(vector, window_b, 1);
//  ASSERT(output_1.length() == 8);
//  Vector<Real> output_1_0_cmp = {0.2000,   0.8000,   -0.9000};
//  Vector<Real> output_1_1_cmp = {0.4000,   1.2000,   -1.2000};
//  Vector<Real> output_1_2_cmp = {0.6000,   1.6000,   -1.5000};
//  Vector<Real> output_1_3_cmp = {0.8000,   2.0000,   -1.8000};
//  Vector<Real> output_1_4_cmp = {1.0000,   2.4000,   -2.1000};
//  Vector<Real> output_1_5_cmp = {1.2000,   2.8000,   -2.4000};
//  Vector<Real> output_1_6_cmp = {1.4000,   3.2000,   -2.7000};
//  Vector<Real> output_1_7_cmp = {1.6000,   3.6000,   -3.0000};
//  ASSERT(IsEqual(output_1[0], output_1_0_cmp));
//  ASSERT(IsEqual(output_1[1], output_1_1_cmp));
//  ASSERT(IsEqual(output_1[2], output_1_2_cmp));
//  ASSERT(IsEqual(output_1[3], output_1_3_cmp));
//  ASSERT(IsEqual(output_1[4], output_1_4_cmp));
//  ASSERT(IsEqual(output_1[5], output_1_5_cmp));
//  ASSERT(IsEqual(output_1[6], output_1_6_cmp));
//  ASSERT(IsEqual(output_1[7], output_1_7_cmp));
//
//  // Test overlapadd
//  std::vector<Vector<Real> > frames_a;
//  Vector<Real> frames_a_0 = {1,2,3};
//  Vector<Real> frames_a_1 = {4,5,6};
//  Vector<Real> frames_a_2 = {7,8,9};
//  frames_a.push_back(frames_a_0);
//  frames_a.push_back(frames_a_1);
//  frames_a.push_back(frames_a_2);
//  Vector<Real> window_d = {0.2,-0.4, 0.3};
//
//  Vector<Real> frames_a_1_cmp = {0.2000, 0, 0.3000, -1.4000, 2.7000};
//  ASSERT(IsEqual(OverlapAdd(frames_a, window_d, 1), frames_a_1_cmp));
//  Vector<Real> frames_a_2_cmp = {0.2000, -0.8000, 1.7000,  -2.0000,
//                                    3.2000, -3.2000, 2.7000};
//  ASSERT(IsEqual(OverlapAdd(frames_a, window_d, 2), frames_a_2_cmp));
//  Vector<Real> frames_a_4_cmp = {0.2000, -0.8000, 0.9000, 0, 0.8000,
//                                      -2.0000, 1.8000, 0, 1.4000, -3.2000,
//                                      2.7000};
//  ASSERT(IsEqual(OverlapAdd(frames_a, window_d, 4), frames_a_4_cmp));
  
  return true;
}
  
  
} // namespace mcl
