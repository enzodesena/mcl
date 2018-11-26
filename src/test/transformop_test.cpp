/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "comparisonop.hpp"
#include "transformop.hpp"
#include "vectorop.hpp"
#include "vector.hpp"
#include <cassert>

namespace mcl {
  
bool TransformOpTest() {
  Vector<Complex<Real>> vector_a(3);
  vector_a[0] = Complex<Real>(1.0, 0.0);
  vector_a[1] = Complex<Real>(0.0, 1.0);
  vector_a[2] = Complex<Real>(1.0, 0.5);
  
  Vector<Complex<Real>> vector_b(3);
  vector_b[0] = Complex<Real>(0.5, 1.0);
  vector_b[1] = Complex<Real>(0.2, -1.0);
  vector_b[2] = Complex<Real>(-1.0, -0.5);
  
  Vector<Complex<Real>> fft_vector_a = Fft(vector_a, 3);
  Vector<Complex<Real>> fft_vector_cmp_a(3);
  fft_vector_cmp_a[0] = Complex<Real>(2.0, 1.5);
  fft_vector_cmp_a[1] = Complex<Real>(0.933012701892219, 0.116025403784439);
  fft_vector_cmp_a[2] = Complex<Real>(0.066987298107781, -1.616025403784439);
  ASSERT(IsApproximatelyEqual(fft_vector_a, fft_vector_cmp_a, VERY_SMALL));
  
  Vector<Complex<Real>> fft_vector_b = Fft(vector_b, 3);
  Vector<Complex<Real>> fft_vector_cmp_b(3);
  fft_vector_cmp_b[0] = Complex<Real>(-0.300000000000000, -0.500000000000000);
  fft_vector_cmp_b[1] = Complex<Real>(0.466987298107781, 0.710769515458674);
  fft_vector_cmp_b[2] = Complex<Real>(1.333012701892219, 2.789230484541326);
  
  ASSERT(IsApproximatelyEqual(fft_vector_b, fft_vector_cmp_b, VERY_SMALL));
  ASSERT(IsApproximatelyEqual(vector_a, Ifft(Fft(vector_a,3),3), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(vector_b, Ifft(Fft(vector_b,3),3), VERY_SMALL));
  
  // Check that when n_points is larger than the size of the vector,
  // then the result is the fft of a zero padded version of the input vector.
  Vector<Complex<Real>> fft_vector_a_4 = Fft(vector_a, 4);
  ASSERT(fft_vector_a_4.length() == 4);
  Vector<Complex<Real>> fft_vector_a_4_cmp(4);
  fft_vector_a_4_cmp[0] = Complex<Real>(2.0000, 1.5000);
  fft_vector_a_4_cmp[1] = Complex<Real>(1.0000, -0.5000);
  fft_vector_a_4_cmp[2] = Complex<Real>(2.0000, -0.5000);
  fft_vector_a_4_cmp[3] = Complex<Real>(-1.0000,-0.5000);
  ASSERT(IsApproximatelyEqual(fft_vector_a_4, fft_vector_a_4_cmp, VERY_SMALL));
  
  // Check that if n_points is smaller than the size of the vector,
  // then the result is the fft of the cut vector
  ASSERT(IsApproximatelyEqual(Fft(vector_a, 2), Fft(Elements(vector_a, 0, 1), 2), VERY_SMALL));
  
  // Test real fft
  Vector<Real> vector_h = {-1.0, 2.0, 3.0};
  ASSERT(IsApproximatelyEqual(Rfft(vector_h, 3),
                 Elements(Fft(CastToComplex<Real>(vector_h), 3), 0, 1), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(Rfft(vector_h, 4),
                 Elements(Fft(CastToComplex<Real>(vector_h), 4), 0, 2), VERY_SMALL));
  
  Vector<Real> vector_i = {-1.0, 2.0, 3.0, 4.0};
  ASSERT(IsApproximatelyEqual(Rfft(vector_h, 4),
                 Elements(Fft(CastToComplex<Real>(vector_h), 4), 0, 2), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(Rfft(vector_h, 5),
                 Elements(Fft(CastToComplex<Real>(vector_h), 5), 0, 2), VERY_SMALL));
  
  // Test Irfft
  ASSERT(IsApproximatelyEqual(Irfft(Rfft(vector_h, 3), 3), vector_h, VERY_SMALL));
  ASSERT(IsApproximatelyEqual(Irfft(Rfft(vector_i, 4), 4), vector_i, VERY_SMALL));
  
  // Test Irfft in strange cases
  Vector<Complex<Real>> vector_l(3);
  vector_l[0] = Complex<Real>(1.0, 0.0);
  vector_l[1] = Complex<Real>(0.0, 1.0);
  vector_l[2] = Complex<Real>(2.0, 0.0);
  Vector<Real> vector_l_3_cmp = {0.333333333333333,
                                      -0.244016935856292,
                                      0.910683602522959};
  ASSERT(IsApproximatelyEqual(Irfft(vector_l, 3), vector_l_3_cmp, VERY_SMALL));
  Vector<Real> vector_l_6_cmp = {0.833333333333333, -0.455341801261480,
                                      -0.455341801261479, 0.833333333333333,
                                      0.122008467928146, 0.122008467928146};
  ASSERT(IsApproximatelyEqual(Irfft(vector_l, 6), vector_l_6_cmp, VERY_SMALL));
  
  Vector<Complex<Real>> vector_m(vector_l); // Try out with even input
  vector_m.PushBack(Complex<Real>(3.0, 0.0));
  Vector<Real> vector_m_4_cmp = {0.75, -0.75, 0.75, 0.25};
  ASSERT(IsApproximatelyEqual(Irfft(vector_m, 4), vector_m_4_cmp, VERY_SMALL));
  Vector<Real> vector_m_3_cmp = {0.333333333333333, -0.244016935856292,
                                      0.910683602522959};
  ASSERT(IsApproximatelyEqual(Irfft(vector_m, 3), vector_m_3_cmp, VERY_SMALL));
  Vector<Real> vector_m_8_cmp = {1.375000000000000, -0.582106781186547,
                                      -0.625000000000000, 0.478553390593274,
                                      -0.125000000000000, 0.832106781186547,
                                      -0.125000000000000, -0.228553390593274};
  ASSERT(IsApproximatelyEqual(Irfft(vector_m, 8), vector_m_8_cmp, VERY_SMALL));
  Vector<Real> vector_m_9_cmp = {1.222222222222222, -0.287886945411706,
                                      -0.858709554352006, 0.363105465825680,
                                      0.042237498424953, 0.194246451014139,
                                      0.748005645285431, -0.421017219679913,
                                      -0.002203563328800};
  ASSERT(IsApproximatelyEqual(Irfft(vector_m, 9), vector_m_9_cmp, VERY_SMALL));
  
  
  // Test hilbert transform
  Vector<Real> vector_c(3);
  vector_c[0] = -0.3;
  vector_c[1] = 0.3;
  vector_c[2] = 2.4;
  
  
  Vector<Complex<Real>> hilbert_vector_c = Hilbert(vector_c);
  Vector<Complex<Real>> hilbert_vector_cmp_c(3);
  hilbert_vector_cmp_c[0] = Complex<Real>(-0.3000000000000000, 1.212435565298214);
  hilbert_vector_cmp_c[1] = Complex<Real>(0.300000000000000, -1.558845726811989);
  hilbert_vector_cmp_c[2] = Complex<Real>(2.400000000000000, 0.346410161513775);
  
  ASSERT(IsApproximatelyEqual(hilbert_vector_c, hilbert_vector_cmp_c, VERY_SMALL));
  
  
  Vector<Real> vector_d(4);
  vector_d[0] = -0.3;
  vector_d[1] = 0.3;
  vector_d[2] = 2.4;
  vector_d[3] = -12.4;
  
  
  Vector<Complex<Real>> hilbert_vector_d = Hilbert(vector_d);
  Vector<Complex<Real>> hilbert_vector_cmp_d(4);
  hilbert_vector_cmp_d[0] = Complex<Real>(-0.300000000000000, -6.350000000000001);
  hilbert_vector_cmp_d[1] = Complex<Real>(0.300000000000001, -1.350000000000000);
  hilbert_vector_cmp_d[2] = Complex<Real>(2.399999999999999, 6.350000000000001);
  hilbert_vector_cmp_d[3] = Complex<Real>(-12.400000000000000, 1.350000000000000);
  
  ASSERT(IsApproximatelyEqual(hilbert_vector_d, hilbert_vector_cmp_d, VERY_SMALL));

  
  Vector<Real> vector_f(4);
  vector_f[0] = 2.5;
  vector_f[1] = 1.3;
  vector_f[2] = -2.4;
  vector_f[3] = -1.0;
  
  
  Vector<Real> rceps_vector_f = RCeps(vector_f);
  Vector<Real> rceps_vector_f_cmp(4);
  rceps_vector_f_cmp[0] = 2.129647179286916e-01;
  rceps_vector_f_cmp[1] = 1.732867951399864e-01;
  rceps_vector_f_cmp[2] = -1.475829040082819e+00;
  rceps_vector_f_cmp[3] = 1.732867951399864e-01;
  ASSERT(IsApproximatelyEqual(rceps_vector_f_cmp, rceps_vector_f, VERY_SMALL));
  
  Vector<Real> mphase_vector_f = MinPhase(vector_f);
  Vector<Real> mphase_vector_f_cmp(4);
  mphase_vector_f_cmp[0] = 2.695552489391163e+00;
  mphase_vector_f_cmp[1] = 9.693272125606057e-01;
  mphase_vector_f_cmp[2] = -2.395552489391163e+00;
  mphase_vector_f_cmp[3] = -8.693272125606056e-01;
  ASSERT(IsApproximatelyEqual(mphase_vector_f_cmp, mphase_vector_f, VERY_SMALL));
  
  
  Vector<Real> vector_g(3);
  vector_g[0] = 2.5;
  vector_g[1] = 0.0;
  vector_g[2] = -2.4;
  
  Vector<Real> rceps_vector_g = RCeps(vector_g);
  Vector<Real> rceps_vector_g_cmp(3);
  rceps_vector_g_cmp[0] = 1.961140220646066e-01;
  rceps_vector_g_cmp[1] = -1.249349557529326e+00;
  rceps_vector_g_cmp[2] = -1.249349557529326e+00;
  ASSERT(IsApproximatelyEqual(rceps_vector_g_cmp, rceps_vector_g, VERY_SMALL));
  
  Vector<Real> mphase_vector_g = MinPhase(vector_g);
  Vector<Real> mphase_vector_g_cmp(3);
  mphase_vector_g_cmp[0] = -1.548105780039157e+00;
  mphase_vector_g_cmp[1] = -1.207601874293259e+00;
  mphase_vector_g_cmp[2] = 2.855707654332416e+00;
  ASSERT(IsApproximatelyEqual(mphase_vector_g_cmp, mphase_vector_g, VERY_SMALL));

  
  Vector<Real> vector_ef(4);
  vector_ef[0] = -0.3;
  vector_ef[1] = 30.3;
  vector_ef[2] = 2.4;
  vector_ef[3] = 12.4;
  
  Vector<Real> vector_ff(4);
  vector_ff[0] = 2.5;
  vector_ff[1] = 1.3;
  vector_ff[2] = -2.4;
  vector_ff[3] = -1.0;
  
  
  Vector<Real> xcorr_vector_f = XCorr(vector_ef, vector_ff);
  Vector<Real> xcorr_vector_f_cmp(7);
  xcorr_vector_f_cmp[0] = 0.3;
  xcorr_vector_f_cmp[1] = -29.58;
  xcorr_vector_f_cmp[2] = -75.51;
  xcorr_vector_f_cmp[3] = 20.48;
  xcorr_vector_f_cmp[4] = 49.110;
  xcorr_vector_f_cmp[5] = 22.12;
  xcorr_vector_f_cmp[6] = 31.0;
  ASSERT(IsApproximatelyEqual(xcorr_vector_f, xcorr_vector_f_cmp, VERY_SMALL));
  
  
  Vector<Real> vector_gf(3);
  vector_gf[0] = 2.5;
  vector_gf[1] = 0.0;
  vector_gf[2] = -2.4;
  
  
  Vector<Real> vector_hf(3);
  vector_hf[0] = -2.5;
  vector_hf[1] = 1.3;
  vector_hf[2] = 2.4;
  
  Vector<Real> xcorr_vector_g_h = XCorr(vector_gf, vector_hf);
  Vector<Real> xcorr_vector_g_h_cmp(5);
  xcorr_vector_g_h_cmp[0] = 6.0;
  xcorr_vector_g_h_cmp[1] = 3.25;
  xcorr_vector_g_h_cmp[2] = -12.01;
  xcorr_vector_g_h_cmp[3] = -3.12;
  xcorr_vector_g_h_cmp[4] = 6.0;
  ASSERT(IsApproximatelyEqual(xcorr_vector_g_h, xcorr_vector_g_h_cmp, VERY_SMALL));
  
  Vector<Real> xcorr_vector_h_g = XCorr(vector_hf, vector_gf);
  Vector<Real> xcorr_vector_h_g_cmp(5);
  xcorr_vector_h_g_cmp[0] = 6.0;
  xcorr_vector_h_g_cmp[1] = -3.12;
  xcorr_vector_h_g_cmp[2] = -12.01;
  xcorr_vector_h_g_cmp[3] = 3.25;
  xcorr_vector_h_g_cmp[4] = 6.0;
  ASSERT(IsApproximatelyEqual(xcorr_vector_h_g, xcorr_vector_h_g_cmp, VERY_SMALL));
  return true;
}
  
} // namespace mcl
