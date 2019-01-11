/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "vectorop.hpp"
#include "comparisonop.hpp"

namespace mcl
{


  
inline bool VectorOpTest()
{

  Vector<double> myvector_a(3);
  myvector_a[0] = 0.1;
  myvector_a[1] = 0.2;
  myvector_a[2] = 0.3;
  ASSERT(myvector_a[0] == 0.1);
  ASSERT(myvector_a[1] == 0.2);
  ASSERT(myvector_a[2] == 0.3);
  ASSERT(myvector_a[2] == 0.3);
  myvector_a[2] = -0.3;
  ASSERT(myvector_a[2] == -0.3);
  
  
  ASSERT(myvector_a.OwnsData());
  Vector<double> myvector_a_copy(myvector_a);
  ASSERT(myvector_a_copy.OwnsData());
  
  Vector<double> myvector_a_reference = MakeReference(myvector_a, 0, 3);
  ASSERT(! myvector_a_reference.OwnsData());
  myvector_a_reference[1] = 1.0;
  ASSERT(myvector_a_reference[1] == 1.0);
  ASSERT(myvector_a[1] == 1.0);
  ASSERT(myvector_a_copy[1] == 0.2); // Should not have changed
  myvector_a[1] = 0.2; // Put it back
  
  Vector<double> myvector_b(3);
  myvector_b[0] = 0.1;
  myvector_b[1] = 0.2;
  myvector_b[2] = 0.3;
  ASSERT(myvector_b[0] == 0.1);
  ASSERT(myvector_b[1] == 0.2);
  ASSERT(myvector_b[2] == 0.3);
  
  myvector_b[2] = -0.3;
  ASSERT(myvector_b[2] == -0.3);

  Vector<Complex<double>> vector_a(3);
  vector_a[0] = Complex<double>(1.0, 0.0);
  vector_a[1] = Complex<double>(0.0, 1.0);
  vector_a[2] = Complex<double>(1.0, 0.5);


  ASSERT(Length(vector_a) == 3);

  Vector<Complex<double>> flip_vector_a = Flip(vector_a);
  Vector<Complex<double>> flip_vector_a_cmp(3);
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
  ASSERT(vector_l.size() == 3);
  Vector<Real> vector_l_cmp(3);
  vector_l_cmp[0] = 0.0;
  vector_l_cmp[1] = 1.0;
  vector_l_cmp[2] = 2.0;
  ASSERT(IsEqual(vector_l_cmp, vector_l));

  Vector<Real> vector_m = LinSpace(0.0, -2.0, 3);
  ASSERT(vector_m.size() == 3);
  Vector<Real> vector_m_cmp(3);
  vector_m_cmp[0] = 0.0;
  vector_m_cmp[1] = -1.0;
  vector_m_cmp[2] = -2.0;
  ASSERT(IsEqual(vector_m_cmp, vector_m));

  Vector<Real> vector_n = LinSpace(1.0, 2.0, 4);
  ASSERT(vector_n.size() == 4);
  Vector<Real> vector_n_cmp(4);
  vector_n_cmp[0] = 1.0;
  vector_n_cmp[1] = 1.33333333333333333333;
  vector_n_cmp[2] = 1.66666666666665666667;
  vector_n_cmp[3] = 2.0;
  ASSERT(IsApproximatelyEqual(vector_n_cmp, vector_n, VERY_SMALL));

  Vector<Real> vector_c(3);
  vector_c[0] = -0.3;
  vector_c[1] = 0.3;
  vector_c[2] = 2.4;

  ASSERT(IsApproximatelyEqual(Sum(vector_c), -0.3+0.3+2.4, VERY_SMALL));

  Vector<Real> pad_vector_c(5);
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
  
//  ASSERT(Downsample(vector_f, 2).size() == 2);
//  ASSERT(Downsample(vector_f, 3).size() == 2);
//  ASSERT(Downsample(vector_f, 4).size() == 1);
//  ASSERT(Downsample(vector_f, 10).size() == 1);
//  ASSERT(Downsample(pad_vector_c_cmp, 2).size() == 3);
//  ASSERT(Downsample(pad_vector_c_cmp, 3).size() == 2);
//  ASSERT(Downsample(pad_vector_c_cmp, 4).size() == 2);
//  ASSERT(Downsample(pad_vector_c_cmp, 5).size() == 1);
//
//  Vector<Real> cmp_downsample_f(2);
//  cmp_downsample_f[0] = 2.5;
//  cmp_downsample_f[1] = -2.4;
//  ASSERT(IsEqual(cmp_downsample_f, Downsample(vector_f, 2)));
//
//  Vector<Real> cmp_downsample_f_3(2);
//  cmp_downsample_f_3[0] = 2.5;
//  cmp_downsample_f_3[1] = -1.0;
//  ASSERT(IsEqual(cmp_downsample_f_3, Downsample(vector_f, 3)));

  Vector<Real> vector_f_sub_0_2 = Subset(vector_f, 0, 2);
  ASSERT(vector_f_sub_0_2.size() == 3);
  Vector<Real> vector_f_sub_0_2_cmp(3);
  vector_f_sub_0_2_cmp[0] = 2.5;
  vector_f_sub_0_2_cmp[1] = 1.3;
  vector_f_sub_0_2_cmp[2] = -2.4;
  ASSERT(IsEqual(vector_f_sub_0_2_cmp, vector_f_sub_0_2));

  Vector<Real> vector_f_sub_1_2 = Subset(vector_f, 1, 2);
  ASSERT(vector_f_sub_1_2.size() == 2);
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
  ASSERT(vector_g_h.size() == 6);
  ASSERT(vector_g_h[0] == 2.5);
  ASSERT(vector_g_h[1] == 0.0);
  ASSERT(vector_g_h[3] == -2.5);



  Vector<Real> unary_vector = UnaryVector((Real) 2.0);
  ASSERT(Length(unary_vector) == 1);
  ASSERT(unary_vector[0] == 2.0);



  Vector<Real> ones_a = Ones<Real>(3);
  ASSERT(Length(ones_a) == 3);
  ASSERT(ones_a[0] == 1.0);
  ASSERT(ones_a[1] == 1.0);
  ASSERT(ones_a[2] == 1.0);

  Vector<Real> vector_o = Zeros<Real>(3);
  vector_o[0] = 1.0;
  vector_o[1] = 2.5;
  vector_o[2] = 4.2;


  ASSERT(IsEqual(vector_o, CircShift(vector_o, 0)));
  ASSERT(IsEqual(vector_o, CircShift(vector_o, 3)));
  ASSERT(IsEqual(vector_o, CircShift(vector_o, -3)));
  ASSERT(IsEqual(CircShift(vector_o, 1), CircShift(vector_o, 4)));
  ASSERT(IsEqual(CircShift(vector_o, -1), CircShift(vector_o, -4)));
  Vector<Real> circ_1_vector_o = CircShift(vector_o, 1);
  Vector<Real> circ_1_vector_o_cmp(3);
  circ_1_vector_o_cmp[0] = vector_o[2];
  circ_1_vector_o_cmp[1] = vector_o[0];
  circ_1_vector_o_cmp[2] = vector_o[1];
  ASSERT(IsEqual(circ_1_vector_o_cmp, circ_1_vector_o));

  Vector<Real> circ_m1_vector_o = CircShift(vector_o, -1);
  Vector<Real> circ_m1_vector_o_cmp(3);
  circ_m1_vector_o_cmp[0] = vector_o[1];
  circ_m1_vector_o_cmp[1] = vector_o[2];
  circ_m1_vector_o_cmp[2] = vector_o[0];
  ASSERT(IsEqual(circ_m1_vector_o_cmp, circ_m1_vector_o));

  Vector<Real> vector_p = LinSpace(0.0, 2.0, 3);
  Vector<Real> vector_q(2);
  vector_q[0] = -1.2;
  vector_q[1] = 4.5;

  ASSERT(IsEqual(Conv(vector_p, vector_q), Conv(vector_q,vector_p)));

  Vector<Real> conv_p_q = Conv(vector_p, vector_q);
  ASSERT(conv_p_q.size() == 4);
  Vector<Real> conv_p_q_cmp(4);
  conv_p_q_cmp[0] = 0.0;
  conv_p_q_cmp[1] = -1.2;
  conv_p_q_cmp[2] = 2.1;
  conv_p_q_cmp[3] = 9;

  ASSERT(IsEqual(conv_p_q_cmp, conv_p_q));

  Vector<Real> vector_r(2);
  vector_r[0] = 3.0;
  vector_r[1] = 6.0;

  Vector<Real> conv_q_r = Conv(vector_r, vector_q);
  ASSERT(conv_q_r.size() == 3);
  Vector<Real> conv_q_r_cmp(3);
  conv_q_r_cmp[0] = -3.6;
  conv_q_r_cmp[1] = 6.3;
  conv_q_r_cmp[2] = 27.0;
  ASSERT(IsApproximatelyEqual(conv_q_r, conv_q_r_cmp, VERY_SMALL));


  Vector<Real> vector_q_padded = Concatenate(Zeros<Real>(0),
                                                  vector_q);
  ASSERT(IsEqual(vector_q_padded, vector_q));




  Vector<Real> colonop_a = ColonOperator<Real>(2, 4);
  Vector<Real> colonop_a_cmp = Zeros<Real>(3);
  colonop_a_cmp[0] = 2.0;
  colonop_a_cmp[1] = 3.0;
  colonop_a_cmp[2] = 4.0;
  ASSERT(IsEqual(colonop_a, colonop_a_cmp));

  Vector<Real> colonop_b = ColonOperator<Real>(-1, 2);
  Vector<Real> colonop_b_cmp = Zeros<Real>(4);
  colonop_b_cmp[0] = -1.0;
  colonop_b_cmp[1] = 0.0;
  colonop_b_cmp[2] = 1.0;
  colonop_b_cmp[3] = 2.0;
  ASSERT(IsEqual(colonop_b, colonop_b_cmp));




  Vector<Int> colonop_c = ColonOperator<Int>(2, 4);
  Vector<Int> colonop_c_cmp = Zeros<Int>(3);
  colonop_c_cmp[0] = 2;
  colonop_c_cmp[1] = 3;
  colonop_c_cmp[2] = 4;
  ASSERT(IsEqual(colonop_c, colonop_c_cmp));


  // Test Poly

  Vector<Complex<Real>> poly_a = Poly(CastToComplex(Ones<Real>(3)));
  Vector<Complex<Real>> poly_a_cmp = Zeros<Complex<Real>>(4);
  poly_a_cmp[0] = Complex<Real>(1.0, 0.0);
  poly_a_cmp[1] = Complex<Real>(-3.0, 0.0);
  poly_a_cmp[2] = Complex<Real>(3.0, 0.0);
  poly_a_cmp[3] = Complex<Real>(-1.0, 0.0);
  ASSERT(IsEqual(poly_a_cmp, poly_a));

  Vector<Complex<Real>> roots_b(4);
  roots_b[0] = Complex<Real>(1.3, 2.0);
  roots_b[1] = Complex<Real>(-1.4, 0.0);
  roots_b[2] = Complex<Real>(2.5, -1.0);
  roots_b[3] = Complex<Real>(0.2, 0.0);
  Vector<Complex<Real>> poly_b = Poly(roots_b);
  Vector<Complex<Real>> poly_b_cmp = Zeros<Complex<Real>>(5);
  poly_b_cmp[0] = Complex<Real>(1.0, 0.0);
  poly_b_cmp[1] = Complex<Real>(-2.600000000000001, -1.000000000000000);
  poly_b_cmp[2] = Complex<Real>(0.410000000000001, 2.500000000000000);
  poly_b_cmp[3] = Complex<Real>(7.364000000000000, 4.720000000000000);
  poly_b_cmp[4] = Complex<Real>(-1.470000000000000, -1.036000000000000);
  ASSERT(IsApproximatelyEqual(poly_b, poly_b_cmp, 0.01));

  ASSERT(IsEqual(Poly(
    mcl::Zeros<Real>(0)),
    UnaryVector(Complex<Real>(1.0,0.0))));
  ASSERT(IsEqual(Poly(
    mcl::Zeros<Complex<Real>>(0)),
    UnaryVector(Complex<Real>(1.0,0.0))));


  // Testing CopyFrom

  Vector<Real> vector_i(3);
  vector_i[0] = 0.1;
  vector_i[1] = -0.5;
  vector_i[2] = 4.0;

  Vector<Real> vector_i_restr_1 = Elements(vector_i, 1, 2);
  ASSERT(vector_i_restr_1.size() == 2);
  Vector<Real> vector_i_restr_1_cmp(2);
  vector_i_restr_1_cmp[0] = -0.5;
  vector_i_restr_1_cmp[1] = 4.0;
  ASSERT(IsEqual(vector_i_restr_1, vector_i_restr_1_cmp));

  Vector<Real> vector_v(3);
  vector_v[0] = 0.1;
  vector_v[1] = -0.5;
  vector_v[2] = 4.0;

  Vector<Real> vector_v_restr_1 = Elements(vector_v, 1, 2);
  ASSERT(vector_v_restr_1.size() == 2);
  Vector<Real> vector_v_restr_1_cmp(2);
  vector_v_restr_1_cmp[0] = -0.5;
  vector_v_restr_1_cmp[1] = 4.0;
  ASSERT(IsEqual(vector_v_restr_1, vector_v_restr_1_cmp));

  Vector<Real> vector_v_restr_2 = Elements(vector_v, 0, 0);
  ASSERT(vector_v_restr_2.size() == 1);
  Vector<Real> vector_v_restr_2_cmp(1);
  vector_v_restr_2_cmp[0] = 0.1;
  ASSERT(IsEqual(vector_v_restr_2, vector_v_restr_2_cmp));

  // Testing GetSegment

  //  Vector<Real> vector_g(3);
  //  vector_g[0] = 2.5;
  //  vector_g[1] = 0.0;
  //  vector_g[2] = -2.4;

  ASSERT(IsEqual(GetSegment(vector_g, 0, 1), UnaryVector((Real) 2.5)));
  ASSERT(IsEqual(GetSegment(vector_g, 1, 1), UnaryVector((Real) 0.0)));
  ASSERT(IsEqual(GetSegment(vector_g, 2, 1), UnaryVector((Real) -2.4)));

  Vector<Real> vector_g_frame_0(2);
  vector_g_frame_0[0] = 2.5;
  vector_g_frame_0[1] = 0.0;
  ASSERT(IsEqual(GetSegment(vector_g, 0, 2), vector_g_frame_0));
  ASSERT(IsEqual(GetSegment(vector_g, 1, 2), UnaryVector((Real) -2.4)));
  ASSERT(IsEqual<Real>(GetSegment(vector_g, 1, 2, true), BinaryVector<Real>(-2.4, 0.0)));
  Vector<Real> hello = GetSegment(vector_g, 2, 2, false);
  ASSERT(IsEqual(GetSegment(vector_g, 2, 2, false), Vector<Real>()));
  ASSERT(IsEqual<Real>(GetSegment(vector_g, 2, 2, true), BinaryVector<Real>(0.0, 0.0)));

  // Testing prod()
  ASSERT(Prod(vector_g_frame_0) == 0.0);
  ASSERT(Prod(vector_i) == -0.2);
  ASSERT(Prod(vector_a) == Complex<double>(-0.5, 1.0));

  // Testing dot product
  ASSERT(Dot(vector_g, vector_v) == -9.350);



  // Testing colon operator
  Vector<Real> vector_z = ColonOperator(0.0, 2.0, 4.0);
  Vector<Real> vector_z_cmp(3);
  vector_z_cmp[0] = 0.0;
  vector_z_cmp[1] = 2.0;
  vector_z_cmp[2] = 4.0;
  ASSERT(IsEqual(vector_z, vector_z_cmp));

  Vector<Real> vector_za = ColonOperator(0.0, 3.0, 4.0);
  Vector<Real> vector_za_cmp(2);
  vector_za_cmp[0] = 0.0;
  vector_za_cmp[1] = 3.0;
  ASSERT(IsEqual(vector_za, vector_za_cmp));

  Vector<Real> vector_aa = ColonOperator(-3.5, 3.0, 3.0);
  Vector<Real> vector_aa_cmp(3);
  vector_aa_cmp[0] = -3.5;
  vector_aa_cmp[1] = -0.5;
  vector_aa_cmp[2] = 2.5;
  ASSERT(IsEqual(vector_aa, vector_aa_cmp));

  Vector<Real> vector_ab = mcl::ColonOperator(-0.001, 0.00025, 0.001);
  Vector<Real> vector_ab_cmp(9);
  vector_ab_cmp[0] = -0.001;
  vector_ab_cmp[1] = -0.00075;
  vector_ab_cmp[2] = -0.0005;
  vector_ab_cmp[3] = -0.00025;
  vector_ab_cmp[4] = 0;
  vector_ab_cmp[5] = 0.00025;
  vector_ab_cmp[6] = 0.0005;
  vector_ab_cmp[7] = 0.00075;
  vector_ab_cmp[8] = 0.001;
  ASSERT(IsEqual(vector_ab, vector_ab_cmp));

  // Testing UnaryVector
  Vector<Real> vector_zd = UnaryVector<Real>(-1.0);
  ASSERT(vector_zd.size() == 1);
  ASSERT(vector_zd[0] == -1.0);

  Vector<Int> vector_ze = UnaryVector<Int>(-2);
  ASSERT(vector_ze.size() == 1);
  ASSERT(vector_ze[0] == -2);

  // Testing BinaryVector
  Vector<Real> vector_zf = BinaryVector<Real>(-1.0, 2.0);
  ASSERT(vector_zf.size() == 2);
  ASSERT(vector_zf[0] == -1.0);
  ASSERT(vector_zf[1] == 2.0);


  // Testing Hanning window
  Vector<Real> vector_hann_3 = Hann<Real>(3);
  Vector<Real> vector_hann_3_cmp(3);
  vector_hann_3_cmp[0] = 0.0;
  vector_hann_3_cmp[1] = 1.0;
  vector_hann_3_cmp[2] = 0.0;
  ASSERT(IsEqual(vector_hann_3, vector_hann_3_cmp));

  Vector<Real> vector_hann_4 = Hann<Real>(4);
  Vector<Real> vector_hann_4_cmp(4);
  vector_hann_4_cmp[0] = 0.0;
  vector_hann_4_cmp[1] = 0.75;
  vector_hann_4_cmp[2] = 0.75;
  vector_hann_4_cmp[3] = 0.0;
  ASSERT(IsApproximatelyEqual(vector_hann_4, vector_hann_4_cmp, VERY_SMALL));

  Vector<Real> vector_hann_5 = Hann<Real>(5);
  Vector<Real> vector_hann_5_cmp(5);
  vector_hann_5_cmp[0] = 0.0;
  vector_hann_5_cmp[1] = 0.5;
  vector_hann_5_cmp[2] = 1.0;
  vector_hann_5_cmp[3] = 0.5;
  vector_hann_5_cmp[4] = 0.0;
  ASSERT(IsApproximatelyEqual(vector_hann_5, vector_hann_5_cmp, VERY_SMALL));

//  ASSERT(IsNonNegative(vector_hann_5_cmp));

  // Testing Hamming window
  Vector<Real> vector_hamming_3 = Hamming<Real>(3);
  Vector<Real> vector_hamming_3_cmp(3);
  vector_hamming_3_cmp[0] = 0.08;
  vector_hamming_3_cmp[1] = 1.0;
  vector_hamming_3_cmp[2] = 0.08;
  ASSERT(IsApproximatelyEqual(vector_hamming_3, vector_hamming_3_cmp, VERY_SMALL));

  Vector<Real> vector_hamming_4 = Hamming<Real>(4);
  Vector<Real> vector_hamming_4_cmp(4);
  vector_hamming_4_cmp[0] = 0.08;
  vector_hamming_4_cmp[1] = 0.77;
  vector_hamming_4_cmp[2] = 0.77;
  vector_hamming_4_cmp[3] = 0.08;
  ASSERT(IsApproximatelyEqual(vector_hamming_4, vector_hamming_4_cmp, VERY_SMALL));

  Vector<Real> vector_hamming_5 = Hann<Real>(5);
  Vector<Real> vector_hamming_5_cmp(5);
  vector_hamming_5_cmp[0] = 0.0;
  vector_hamming_5_cmp[1] = 0.5;
  vector_hamming_5_cmp[2] = 1.0;
  vector_hamming_5_cmp[3] = 0.5;
  vector_hamming_5_cmp[4] = 0.0;
  ASSERT(IsApproximatelyEqual(vector_hamming_5, vector_hamming_5_cmp, VERY_SMALL));

//  ASSERT(IsNonNegative(vector_hamming_5));

  // Testing Tukey Window
  ASSERT(IsApproximatelyEqual(TukeyWin(4, 1.0), vector_hann_4_cmp, VERY_SMALL));
  Vector<Real> vector_tukey_1 = TukeyWin(5,0.6);
  Vector<Real> vector_tukey_1_cmp(5);
  vector_tukey_1_cmp[0] = 0.0;
  vector_tukey_1_cmp[1] = 0.933012701892219;
  vector_tukey_1_cmp[2] = 1.0;
  vector_tukey_1_cmp[3] = 0.933012701892219;
  vector_tukey_1_cmp[4] = 0.0;
  ASSERT(IsApproximatelyEqual(vector_tukey_1, vector_tukey_1_cmp, VERY_SMALL));

//  Vector<Real> aa = TukeyWin<Real>(1, 0.6);
  ASSERT(TukeyWin(1, 0.6).size() == 1);
  ASSERT(TukeyWin(0, 0.6).size() == 0);
  ASSERT(IsEqual(TukeyWin(1, 0.6), UnaryVector<Real>(1.0)));

  ASSERT(IsEqual(TukeyWin(6, 0.0), Ones<Real>(6)));
  ASSERT(IsEqual(TukeyWin(6, -2.0), Ones<Real>(6)));

  // Testing norm
  Vector<Real> vector_ba(4);
  vector_ba[0] = -1.2;
  vector_ba[1] = 2.3;
  vector_ba[2] = 3.4;
  vector_ba[3] = -5.0;
  ASSERT(IsApproximatelyEqual(Norm(vector_ba, 2.0), 6.579513659838392));
  ASSERT(IsApproximatelyEqual(Norm(vector_ba, 1.0), 11.899999999999999));
  ASSERT(IsApproximatelyEqual(Norm(vector_ba, 2.4), 6.056130782634900));




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
//  interleaves_e_f_cmp[] = vector_e[0]);
//  interleaves_e_f_cmp[] = vector_f[0]);
//  interleaves_e_f_cmp[] = vector_e[1]);
//  interleaves_e_f_cmp[] = vector_f[1]);
//  interleaves_e_f_cmp[] = vector_e[2]);
//  interleaves_e_f_cmp[] = vector_f[2]);
//  interleaves_e_f_cmp[] = vector_e[3]);
//  interleaves_e_f_cmp[] = vector_f[3]);
//  ASSERT(IsEqual(Interleave(vector_e, vector_f), interleaves_e_f_cmp));
//
//
//  // Test Enframe
//  Vector<Real> vector = LinSpace(1, 10, 10);
//  Vector<Real> window = {0.2, 0.4};
//  Vector<Vector<Real> > output_3 = Enframe(vector, window, 3);
//  ASSERT(output_3.size() == 3);
//  Vector<Real> output_3_0_cmp = {0.2, 0.8};
//  Vector<Real> output_3_1_cmp = {0.8, 2.0};
//  Vector<Real> output_3_2_cmp = {1.4, 3.2};
//  ASSERT(IsEqual(output_3[0], output_3_0_cmp));
//  ASSERT(IsEqual(output_3[1], output_3_1_cmp));
//  ASSERT(IsEqual(output_3[2], output_3_2_cmp));
//
//
//  Vector<Vector<Real> > output_full = Enframe(Ones<Real>(4), window, 2);
//  ASSERT(output_full.size()==2);
//  Vector<Real> output_full_0_cmp = {0.2, 0.4};
//  Vector<Real> output_full_1_cmp = {0.2, 0.4};
//  ASSERT(IsEqual(output_full[0], output_full_0_cmp));
//  ASSERT(IsEqual(output_full[1], output_full_1_cmp));
//
//
//  Vector<Real> window_b;
//  window_b[] = 0.2);
//  window_b[] = 0.4);
//  window_b[] = -0.3);
//  Vector<Vector<Real> > output_1 = Enframe(vector, window_b, 1);
//  ASSERT(output_1.size() == 8);
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
//  Vector<Vector<Real> > frames_a;
//  Vector<Real> frames_a_0 = {1,2,3};
//  Vector<Real> frames_a_1 = {4,5,6};
//  Vector<Real> frames_a_2 = {7,8,9};
//  frames_a[] = frames_a_0);
//  frames_a[] = frames_a_1);
//  frames_a[] = frames_a_2);
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
  
  
  // Testing vector reference
  Vector<Real> referenced(3, 0.0);
  referenced[1] = 2.0;
  ASSERT(referenced.size() == 3);
  ASSERT(referenced[0] == 0.0);
  ASSERT(referenced[1] == 2.0);
  ASSERT(referenced[2] == 0.0);
  Vector<Real> reference_a = MakeReference(referenced, 0, 3);
  reference_a[0] = 1.0;
  reference_a[1] = 3.0;
  reference_a[2] = 5.0;
  ASSERT(reference_a[0] == 1.0);
  ASSERT(reference_a[1] == 3.0);
  ASSERT(reference_a[2] == 5.0);
  ASSERT(referenced[0] == 1.0);
  ASSERT(referenced[1] == 3.0);
  ASSERT(referenced[2] == 5.0);
  ASSERT(reference_a.size() == 3);
  
  Vector<Real> reference_b = MakeReference(referenced, 0, 2);
  ASSERT(reference_b.size() == 2);
  ASSERT(reference_b[0] == 1.0);
  ASSERT(reference_b[1] == 3.0);
  
  Vector<Real> reference_c = MakeReference(referenced, 1, 2);
  ASSERT(reference_c.size() == 2);
  ASSERT(reference_c[0] == 3.0);
  ASSERT(reference_c[1] == 5.0);
  
  Vector<Real> reference_d = MakeReference(referenced, 2, 1);
  ASSERT(reference_d.size() == 1);
  ASSERT(reference_d[0] == 5.0);
  
  referenced[2] = 10.0;
  ASSERT(referenced[2] == 10.0);
  ASSERT(reference_c[1] == 10.0);
  ASSERT(reference_d[0] == 10.0);
  
  reference_d[0] = -5.0;
  ASSERT(referenced[2] == -5.0);
  ASSERT(reference_c[1] == -5.0);
  ASSERT(reference_d[0] == -5.0);
  
  auto iter_d = reference_d.begin();
  ASSERT(*iter_d == -5.0);
  AVOID_UNUSED_WARNING(iter_d);
  
  auto iter_a = reference_a.begin() + 2;
  ASSERT(*iter_a == -5.0);
  *iter_a = -2.0;
  ASSERT(reference_d[0] == -2.0);
  
  
  
  // Testing general forward iterator
  Vector<Real> referenced_b(3, 1.0);
  referenced_b[0] = 2.0;
  auto iter_fwd = referenced_b.GetFwdIteratorBegin();
  auto iter_fwd_end = referenced_b.GetFwdIteratorEnd();
  ASSERT(iter_fwd != iter_fwd_end);
  ASSERT(*iter_fwd == 2.0);
  
  ++iter_fwd;
  ASSERT(*iter_fwd == 1.0);
  ASSERT(iter_fwd != iter_fwd_end);
  
  ++iter_fwd;
  ASSERT(*iter_fwd == 1.0);
  ASSERT(++iter_fwd == iter_fwd_end);
  return true;
}


template<typename Iterator>
void FwdIteratorSpeedRoutine(Iterator iter, Iterator end)
{
  int k = 0;
  while (iter != end)
  {
    *iter++ = k++;
  }
}

inline void FwdIteratorSpeedTests()
{
  const size_t vector_length = 1000000;
  Vector<Real> vec_a(vector_length, 1.0);
  
  auto iter_fwd = vec_a.GetFwdIteratorBegin();
  auto iter_fwd_end = vec_a.GetFwdIteratorEnd();
  
  
  clock_t launch = clock();
  FwdIteratorSpeedRoutine
  (
    vec_a.begin(),
    vec_a.end());
  clock_t done = clock();
  std::cout<<"Iterator (std::vector): "<<(done - launch) / ((Real) CLOCKS_PER_SEC)<<"s \n";
  
  launch = clock();
  FwdIteratorSpeedRoutine
  (
    vec_a.GetFwdIteratorBegin(),
    vec_a.GetFwdIteratorEnd());
  done = clock();
  
  std::cout<<"Iterator (MCL): "<<(done - launch) / ((Real) CLOCKS_PER_SEC)<<"s \n";
  
  
  
}
  
} // namespace mcl
