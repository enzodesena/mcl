/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "firfilter.hpp"
#include "butter.hpp"
#include "randomop.hpp"

namespace mcl
{

//inline bool IirFilterTest()
//{
//  DigitalFilter<Real> filter_a(GainFilter<Real>(1.0));
//  Vector<Real> output_value(1);
//  filter_a.Filter(UnaryVector<Real>(1.2), output_value);
//  ASSERT(IsApproximatelyEqual(output_value[0], 1.2));
//
//  DigitalFilter<Real> filter_b(GainFilter<Real>(0.76));
//  Vector<Real> output_value_b(1);
//  filter_b.Filter(UnaryVector<Real>(1.2), output_value_b);
//  ASSERT(IsApproximatelyEqual(output_value_b[0], 0.912));
//
//  Vector<Real> B(4);
//  Vector<Real> A(4);
//  B[0] = 0.76;
//  B[1] = -3.06;
//  B[2] = 1.76;
//  B[3] = 1.76;
//  A[0] = 1.0;
//  A[1] = -0.5;
//  A[2] = 0.23;
//  A[3] = 0.75;
//
//  IirFilter<Real> filter_c(B,A);
//
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(1.2), 0.912));
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(0.2), -3.064));
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(0.5), 0.13824));
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(-1.2), 0.11184));
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(0.0), 7.226124799999999999999999999999));
//  ASSERT(IsApproximatelyEqual(filter_c.Filter(0.0), 2.251659199999999999999999999999));
//
//  IirFilter<Real> filter_d = WallFilter<Real>(kCarpetPile, 44100.0);
//  ASSERT(IsApproximatelyEqual(filter_d.Filter(1.0), 0.562666833756030));
//  ASSERT(IsApproximatelyEqual(filter_d.Filter(0.5), 0.315580130841020));
//
//  IirFilter<Real> filter_e = WallFilter<Real>(kCarpetCotton, 44100.0);
//  ASSERT(IsApproximatelyEqual(filter_e.Filter(1.0), 0.687580695329600));
//  ASSERT(IsApproximatelyEqual(filter_e.Filter(0.5), 0.322032066558733));
//
//  IirFilter<Real> filter_f = WallFilter<Real>(kWallBricks, 44100.0);
//  ASSERT(IsApproximatelyEqual(filter_f.Filter(1.0), 0.978495798553620));
//  ASSERT(IsApproximatelyEqual(filter_f.Filter(0.5), 0.490594444043047));
//
//  IirFilter<Real> filter_g = WallFilter<Real>(kCeilingTile, 44100.0);
//  ASSERT(IsApproximatelyEqual(filter_g.Filter(1.0), 0.168413736374283));
//  ASSERT(IsApproximatelyEqual(filter_g.Filter(0.5), 0.151668254946940));
//
//
//  Vector<Real> input_a(4);
//  input_a[0] = 0.6;
//  input_a[1] = -3.5;
//  input_a[2] = 5.6;
//  input_a[3] = 2.3;
//
//  // Testing pinkifier filter
//  IirFilter<Real> pinkifier = PinkifierFilter<Real>();
//  Vector<Real> output_e = pinkifier.Filter(input_a);
//  Vector<Real> output_e_cmp(input_a.size());
//  output_e_cmp[0] = 0.600000000000000;
//  output_e_cmp[1] = -3.152233220000000;
//  output_e_cmp[2] = 3.815449359516707;
//  output_e_cmp[3] = 4.322130531286722;
//  ASSERT(IsApproximatelyEqual(output_e, output_e_cmp, VERY_SMALL));
//
//
//
//  Vector<Real> B_d(3);
//  B_d[0] = 1.0;
//  B_d[1] = -2.0;
//  B_d[2] = 1.0;
//  Vector<Real> A_d(3);
//  A_d[0] = 1.005844676087000;
//  A_d[1] = -1.999977314492666;
//  A_d[2] = 0.994178009420333;
//  B_d = Multiply(B_d, 1/A_d[0]);
//  A_d = Multiply(A_d, 1/A_d[0]);
//  IirFilter<Real> filter_l(B_d, A_d);
//  ASSERT(IsApproximatelyEqual(B_d, filter_l.B(), VERY_SMALL));
//  ASSERT(IsApproximatelyEqual(A_d, filter_l.A(), VERY_SMALL));
//
//  Vector<Real> signal_d = Zeros<Real>(4);
//  signal_d[0] = 0.989949493661167;
//  Vector<Real> signal_d_out_cmp(4);
//  signal_d_out_cmp[0] = 0.984197179938686;
//  signal_d_out_cmp[1] = -0.011459974617699;
//  signal_d_out_cmp[2] = -0.011370929428126;
//  signal_d_out_cmp[3] = -0.011282404149780;
//
//  DigitalFilter<Real> filter_poly(filter_l);
//  Vector<Real> output_d(signal_d.size());
//  filter_poly.Filter(signal_d, output_d);
//  ASSERT(IsApproximatelyEqual(output_d, signal_d_out_cmp, VERY_SMALL));
//
//  // Testing Reset()
//  filter_l.Reset();
//  ASSERT(IsApproximatelyEqual(filter_l.Filter(0.0), 0.0, VERY_SMALL));
//
//  Vector<Real> impulse_resp_2(3);
//  impulse_resp_2[0] = 0.2;
//  impulse_resp_2[1] = -0.1;
//  impulse_resp_2[2] = 2.5;
////
////  DigitalFilter<Real> filter_m(impulse_resp_2);
////  ASSERT(! IsApproximatelyEqual(filter_m.Filter(1.0), 0.0));
////  filter_m.Reset();
////  ASSERT(IsApproximatelyEqual(filter_m.Filter(0.0), 0.0));
////
////
//  // Testing butterworth filter
//  IirFilter<Real> butter_a = Butter<Real>(3, 0.2, 0.45);
//  Vector<Real> butter_a_num_cmp = Zeros<Real>(7);
//  butter_a_num_cmp[0] = 0.031689343849711;
//  butter_a_num_cmp[2] = -0.095068031549133;
//  butter_a_num_cmp[4] = 0.095068031549133;
//  butter_a_num_cmp[6] = -0.031689343849711;
//
//  Vector<Real> butter_a_den_cmp(7);
//  butter_a_den_cmp[0] = 1.000000000000000;
//  butter_a_den_cmp[1] = -2.521796622886441;
//  butter_a_den_cmp[2] = 3.643067063269880;
//  butter_a_den_cmp[3] = -3.325285581665978;
//  butter_a_den_cmp[4] = 2.149206132889376;
//  butter_a_den_cmp[5] = -0.850496842492471;
//  butter_a_den_cmp[6] = 0.197825187264320;
//
//  ASSERT(IsApproximatelyEqual(butter_a.B(), butter_a_num_cmp, VERY_SMALL));
//  ASSERT(IsApproximatelyEqual(butter_a.A(), butter_a_den_cmp, VERY_SMALL));
//
//
//  IirFilter<Real> butter_b = Butter<Real>(2, 0.12, 0.79);
//  Vector<Real> butter_b_num_cmp = Zeros<Real>(5);
//  butter_b_num_cmp[0] = 0.469043625796947;
//  butter_b_num_cmp[2] = -0.938087251593893;
//  butter_b_num_cmp[4] = 0.469043625796947;
//
//
//  Vector<Real> butter_b_den_cmp(5);
//  butter_b_den_cmp[0] = 1.000000000000000;
//  butter_b_den_cmp[1] = -0.388787442245741;
//  butter_b_den_cmp[2] = -0.583519141064213;
//  butter_b_den_cmp[3] = 0.041607774454425;
//  butter_b_den_cmp[4] = 0.243288940651677;
//
//  ASSERT(IsApproximatelyEqual(butter_b.B(), butter_b_num_cmp, VERY_SMALL));
//  ASSERT(IsApproximatelyEqual(butter_b.A(), butter_b_den_cmp, VERY_SMALL));
//
//  IirFilter<Real> filter_i;
//  ASSERT(IsApproximatelyEqual(filter_i.Filter(1.2), 1.2));
//  ASSERT(IsApproximatelyEqual(filter_i.Filter(-0.2), -0.2));
//
//  // Testing octave filter
//  IirFilter<Real> octave_a = OctaveFilter(3, 4000.0, 44100.0);
//  Vector<Real> octave_a_num_cmp = Zeros<Real>(7);
//  octave_a_num_cmp[0] = 0.005020133201471;
//  octave_a_num_cmp[2] = -0.015060399604412;
//  octave_a_num_cmp[4] = 0.015060399604412;
//  octave_a_num_cmp[6] = -0.005020133201471;
//
//  Vector<Real> octave_a_den_cmp(7);
//  octave_a_den_cmp[0] = 1.000000000000000;
//  octave_a_den_cmp[1] = -4.397041740651781;
//  octave_a_den_cmp[2] = 8.729527405676500;
//  octave_a_den_cmp[3] = -9.889119467962011;
//  octave_a_den_cmp[4] = 6.737413381809715;
//  octave_a_den_cmp[5] = -2.619423015108258;
//  octave_a_den_cmp[6] = 0.460896610043675;
//
//  ASSERT(IsApproximatelyEqual(octave_a.B(), octave_a_num_cmp, VERY_SMALL));
//  ASSERT(IsApproximatelyEqual(octave_a.A(), octave_a_den_cmp, VERY_SMALL));
//
//
//  // Testing iir filter bank
//  IirFilterBank octave_bank_a = OctaveFilterBank(3, 1, 4000.0, 44100.0);
//  ASSERT(IsApproximatelyEqual(octave_bank_a.Filter(1.25)[0], octave_a.Filter(1.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_a.Filter(0.25)[0], octave_a.Filter(0.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_a.Filter(5.0)[0], octave_a.Filter(5.0)));
//  ASSERT(octave_bank_a.Filter(1.25).size() == 1);
//  ASSERT(octave_bank_a.num_filters() == 1);
//
//  IirFilterBank octave_bank_b = OctaveFilterBank(3, 2, 2000.0, 44100.0);
//  octave_a.Reset();
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(1.25)[1], octave_a.Filter(1.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(0.25)[1], octave_a.Filter(0.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(5.0)[1], octave_a.Filter(5.0)));
//  ASSERT(octave_bank_b.Filter(1.25).size() == 2);
//  ASSERT(octave_bank_b.num_filters() == 2);
//
//  octave_a.Reset();
//  octave_bank_b.Reset();
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(1.25)[1], octave_a.Filter(1.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(0.25)[1], octave_a.Filter(0.25)));
//  ASSERT(IsApproximatelyEqual(octave_bank_b.Filter(5.0)[1], octave_a.Filter(5.0)));
//  ASSERT(octave_bank_b.Filter(1.25).size() == 2);
//
//  return true;
//}


//inline void IirFilterSpeedTests() {
//
//  constexpr size_t filter_order = 3;
//  constexpr Real input_seconds = 10.0;
//  constexpr Real sampling_frequency = 44100;
//  
//  const size_t num_input_samples = round(input_seconds*sampling_frequency);
//  
//  RandomGenerator random_generator;
//
//  Vector<Real> B = random_generator.Rand(filter_order);
//  Vector<Real> A = random_generator.Rand(filter_order);
//  A[0] = 1.0;
//  
//  IirFilter<Real> iir_filter(B, A);
//  Vector<Real> input = random_generator.Rand(num_input_samples);
//  Vector<Real> output(num_input_samples);
//
//  clock_t launch=clock();
//  for (size_t i=0; i<input.size(); ++i) {
//    output[i] = iir_filter.Filter(input[i]);
//  }
//  clock_t done=clock();
//
//  std::cout<<"Iir filter speed (sequential; filter length is not power of 2): "<<
//  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_seconds*100<<"% \n";
//}
  
  
inline bool DigitalFilterTest()
{
  Vector<Real>  ir = {0.1, 0.2, 0.3};
  DigitalFilter<Real> filter_lasplita(ir);
  Vector<Real> input = {1, 2, 3, 4, 5, 6, 7};
  Vector<Real> output_lasplita_a(3,0.0);
  Vector<Real> cmp_lasplita_a = { 0.1, 0.4, 1.0 };
#ifdef MCL_APPLE_ACCELERATE
//  filter_lasplita.FilterAppleDsp(&input.data()[0], 3, output_lasplita_a);
//  ASSERT(IsEqual(cmp_lasplita_a, output_lasplita_a));
#endif

  Vector<Real> output_lasplita_b(4, 0.0);
  Vector<Real> cmp_lasplita_b = { 1.6, 2.2, 2.8, 3.4 };
#ifdef MCL_APPLE_ACCELERATE
//  filter_lasplita.FilterAppleDsp(&input.data()[3], 4, output_lasplita_b);
//  ASSERT(IsEqual(cmp_lasplita_b, output_lasplita_b));
#endif

//  filter_lasplita.Reset();
//  Vector<Real,kReference> refa(input, 0, 3);
//  Vector<Real,kReference> refb(input, 0, 3);
//  filter_lasplita.Filter(refa, output_lasplita_a);
//  ASSERT(IsEqual(cmp_lasplita_a, output_lasplita_a));
//  filter_lasplita.Filter(refb, output_lasplita_b);
//  ASSERT(IsEqual(cmp_lasplita_b, output_lasplita_b));
  
  

  Vector<Real> impulse_resp(3);
  impulse_resp[0] = 0.2;
  impulse_resp[1] = -0.1;
  impulse_resp[2] = 2.5;

  Vector<Real> impulse(3);
  impulse[0] = 1.0;
  impulse[1] = 0.0;
  impulse[2] = 0.0;

  DigitalFilter<Real> filter_a(impulse_resp);
  Vector<Real> output_aa(impulse.size());
  filter_a.Filter(impulse, output_aa);
  ASSERT(IsEqual(output_aa, impulse_resp));

  DigitalFilter<Real> filter_b(impulse_resp);
  Vector<Real> input_a(4);
  input_a[0] = 0.6;
  input_a[1] = -3.5;
  input_a[2] = 5.6;
  input_a[3] = 2.3;

  Vector<Real> output_a_cmp(4);
  output_a_cmp[0] = 0.1200;
  output_a_cmp[1] = -0.7600;
  output_a_cmp[2] = 2.9700;
  output_a_cmp[3] = -8.8500;
  Vector<Real> output_a(input_a.size());
  filter_b.Filter(input_a, output_a);
  ASSERT(IsApproximatelyEqual(output_a, output_a_cmp, VERY_SMALL));

  DigitalFilter<Real> filter_c(impulse_resp);
  ASSERT(IsApproximatelyEqual(filter_c.Filter(0.6), 0.1200));
  ASSERT(IsApproximatelyEqual(filter_c.Filter(-3.5), -0.7600));

//#ifdef MCL_APPLE_ACCELERATE
//  DigitalFilter<Real> filter_ca(impulse_resp);
//  ASSERT(IsApproximatelyEqual(filter_ca.FilterAppleDsp(0.6), 0.1200));
//  ASSERT(IsApproximatelyEqual(filter_ca.FilterAppleDsp(-3.5), -0.7600));
//#endif

  DigitalFilter<Real> filter_cb(impulse_resp);
  ASSERT(IsApproximatelyEqual(filter_cb.Filter(0.6), 0.1200));
  ASSERT(IsApproximatelyEqual(filter_cb.Filter(-3.5), -0.7600));

//#ifdef MCL_APPLE_ACCELERATE
//  DigitalFilter<Real> filter_cd(impulse_resp);
//  ASSERT(IsApproximatelyEqual(filter_cd.FilterAppleDsp(0.6), 0.1200));
//  ASSERT(IsApproximatelyEqual(filter_cd.FilterStraight(-3.5), -0.7600));
//#endif

  DigitalFilter<Real> filter_ce(impulse_resp);
  ASSERT(IsApproximatelyEqual(filter_ce.Filter(0.6), 0.1200));
//#ifdef MCL_APPLE_ACCELERATE
//  ASSERT(IsApproximatelyEqual(filter_ce.FilterAppleDsp(-3.5), -0.7600));
//#endif

  // Testing copy constructor
  DigitalFilter<Real> filter_d = filter_c;
  ASSERT(IsApproximatelyEqual(filter_d.Filter(5.6), 2.9700));
  ASSERT(IsApproximatelyEqual(filter_d.Filter(2.3), -8.8500));

  ASSERT(IsApproximatelyEqual(filter_c.Filter(5.6), 2.9700));
  ASSERT(IsApproximatelyEqual(filter_c.Filter(2.3), -8.8500));


  DigitalFilter<Real> filter_i;
  ASSERT(IsApproximatelyEqual(filter_i.Filter(1.2), 1.2));
  ASSERT(IsApproximatelyEqual(filter_i.Filter(-0.2), -0.2));

  Vector<Real> impulse_resp_b = {0.1, 0.3, -0.2, 1.2, -4.5, 0.0, -2.1, -1.2};
  DigitalFilter<Real> filter_l(impulse_resp_b);
  Vector<Real> input_b = {0.3377, 0.9001, 0.3692, 0.1112, 0.7803, 0.3897, 0.2417, 0.4039, 0.0965, 0.1320, 0.9421, 0.9561};
  Vector<Real> output_b_cmp = {0.033770000000000, 0.191320000000000, 0.239410000000000, 0.347100000000000, -0.401980000000000, -3.356590000000000, -2.252110000000000, -1.824530000000000, -4.816670000000000, -2.178800000000000, -2.260530000000000, -3.104640000000000};
  Vector<Real> output_b(input_b.size());
  filter_l.Filter(input_b, output_b);
  ASSERT(IsApproximatelyEqual(output_b_cmp, output_b, VERY_SMALL));

  filter_l.SetStateToZero();
  for (Int i=0; i<(Int)input_b.size(); ++i) {
    ASSERT(IsApproximatelyEqual(filter_l.Filter(input_b[i]), output_b_cmp[i]));
  }

//#ifdef MCL_APPLE_ACCELERATE
//  DigitalFilter<Real> filter_la(impulse_resp_b);
//  Real cmp_la[input_b.size()];
//  filter_la.FilterAppleDsp(input_b.data(), input_b.size(), cmp_la);
//  ASSERT(IsEqual(output_b_cmp, cmp_la));
//
//  DigitalFilter<Real> filter_lasplit(impulse_resp_b);
//  Real cmp_lasplit_a[4];
//  Real cmp_lasplit_b[8];
//  filter_lasplit.FilterAppleDsp(&input_b.data()[0], 4, cmp_lasplit_a);
//  ASSERT(IsEqual(&output_b_cmp.data()[0], cmp_lasplit_a, 4));
//  filter_lasplit.FilterAppleDsp(&input_b.data()[4], 8, cmp_lasplit_b);
//  ASSERT(IsEqual(&output_b_cmp.data()[4], cmp_lasplit_b, 8));
//#endif

  DigitalFilter<Real> filter_lb(impulse_resp_b);
  Vector<Real> cmp_lb(input_b.size(), 0.0);
  filter_lb.Filter(input_b, cmp_lb);
  ASSERT(IsApproximatelyEqual(output_b_cmp, cmp_lb, VERY_SMALL));

  Vector<Real> input_c = {0.8147, 0.9058, 0.1270, 0.9134, 0.6324, 0.0975, 0.2785, 0.5469, 0.9575, 0.9649, 0.1576, 0.9706, 0.9572, 0.4854, 0.8003, 0.1419, 0.4218, 0.9157, 0.7922, 0.9595};
  Vector<Real> impulse_resp_c = {0.6948, 0.3171, 0.9502, 0.0344, 0.4387, 0.3816, 0.7655, 0.7952, 0.1869, 0.4898, 0.4456, 0.6463, 0.7094, 0.7547, 0.2760, 0.6797};
  DigitalFilter<Real> filter_m(impulse_resp_c);
  Vector<Real> output_c_cmp = {0.566053560000000, 0.887691210000000, 1.149596720000000, 1.563618860000000, 1.238274470000000, 1.848822500000000, 1.881767519999999, 2.373108650000000, 2.702443100000000, 3.155909820000000, 3.544349419999999, 3.760939330000000, 3.860796740000000, 5.071760400000001, 5.228588220000000, 5.070855620000001, 5.216075850000000, 4.336750739999999, 5.636061180000000,5.665156830000000};
  Vector<Real> output_c(input_c.size());
  filter_m.Filter(input_c, output_c);
  ASSERT(IsApproximatelyEqual(output_c_cmp, output_c, VERY_SMALL));


  // Various attempt to check that the batch processing does not mess up
  filter_m.SetStateToZero();
  Vector<Real> input_c_sub_a(input_c.begin(), input_c.begin()+16);
  Vector<Real> output_c_cmp_sub_a(output_c_cmp.begin(), output_c_cmp.begin()+16);
  Vector<Real> output_cc(input_c_sub_a.size());
  filter_m.Filter(input_c_sub_a, output_cc);
  ASSERT(IsApproximatelyEqual(output_cc, output_c_cmp_sub_a, VERY_SMALL));
  ASSERT(IsApproximatelyEqual(filter_m.Filter(input_c[16]), output_c_cmp[16], VERY_SMALL));

  Vector<Real> input_c_sub_b(input_c.begin()+17, input_c.end());
  Vector<Real> output_c_cmp_sub_b(output_c_cmp.begin()+17, output_c_cmp.end());
  Vector<Real> output_c_sub_b(input_c_sub_b.size());
  filter_m.Filter(input_c_sub_b, output_c_sub_b);
  ASSERT(IsApproximatelyEqual(output_c_sub_b, output_c_cmp_sub_b, VERY_SMALL));

  //
  filter_m.SetStateToZero();
  ASSERT(IsEqual(filter_m.Filter(input_c[0]), output_c_cmp[0]));
  ASSERT(IsApproximatelyEqual(filter_m.Filter(input_c[1]), output_c_cmp[1]));
  ASSERT(IsApproximatelyEqual(filter_m.Filter(input_c[2]), output_c_cmp[2]));
  Vector<Real> input_c_sub_ab(input_c.begin()+3, input_c.begin()+19);
  Vector<Real> output_c_cmp_sub_ab(output_c_cmp.begin()+3, output_c_cmp.begin()+19);
  Vector<Real> output_ccc(input_c_sub_ab.size());
  filter_m.Filter(input_c_sub_ab, output_ccc);
  ASSERT(IsApproximatelyEqual(output_ccc, output_c_cmp_sub_ab, VERY_SMALL));
  ASSERT(IsApproximatelyEqual(filter_m.Filter(input_c[19]), output_c_cmp[19]));


  //
  Vector<Real> impulse_response_k = {0.8147, 0.9058, 0.1270, 0.9134, 0.6324};
  DigitalFilter<Real> filter_k(impulse_response_k);
  Vector<Real> input_k = input_c;
  Vector<Real> output_k_cmp = {0.663736090000000, 1.475910520000000, 1.027407440000000, 1.718367160000000, 2.701277000000000, 1.457092690000000, 1.310138610000000, 1.865475550000000, 1.799813030000000, 2.038904730000000, 1.799667500000000, 2.276484260000000, 3.165878180000000, 2.139907940000000, 2.199456410000000, 2.390277390000000, 1.622509220000000, 2.184069510000000, 2.164136180000000, 2.090582990000000};
  Vector<Real> output_k(input_k.size());
  filter_k.Filter(input_k, output_k);
  ASSERT(IsApproximatelyEqual(output_k, output_k_cmp, VERY_SMALL));

  //
  filter_k.SetStateToZero();
  for (Int i=0; i<(Int)input_c.size()-1; ++i) {
    ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[i]), output_k_cmp[i]));
  }

  //
  filter_k.SetStateToZero();
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[0]), output_k_cmp[0]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[1]), output_k_cmp[1]));
  Vector<Real> input_k_sub_a = Vector<Real>(input_k.begin()+2,
                                                      input_k.begin()+7);
  Vector<Real> output_k_cmp_sub_a = Vector<Real>(output_k_cmp.begin()+2,
                                                           output_k_cmp.begin()+7);
  
  Vector<Real> output_k_sub_a(input_k_sub_a.size());
  filter_k.Filter(input_k_sub_a, output_k_sub_a);
  ASSERT(IsApproximatelyEqual(output_k_sub_a, output_k_cmp_sub_a, VERY_SMALL));

  Vector<Real> output_k_sub_b(2);
  filter_k.Filter(Vector<Real>(input_k.begin()+7, input_k.begin()+9), output_k_sub_b);
  ASSERT(IsApproximatelyEqual(output_k_sub_b,
                 Vector<Real>(output_k_cmp.begin()+7, output_k_cmp.begin()+9), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[9]), output_k_cmp[9]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[10]), output_k_cmp[10]));
  
  Vector<Real> output_k_sub_c(9);
  filter_k.Filter(Vector<Real>(input_k.begin()+11, input_k.begin()+20), output_k_sub_c);
  ASSERT(IsApproximatelyEqual(output_k_sub_c, Vector<Real>(output_k_cmp.begin()+11, output_k_cmp.begin()+20), VERY_SMALL));

  //
  filter_k.SetStateToZero();
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[0]), output_k_cmp[0]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[1]), output_k_cmp[1]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[2]), output_k_cmp[2]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[3]), output_k_cmp[3]));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[4]), output_k_cmp[4]));
  Vector<Real> output_k_sub_d(5);
  filter_k.Filter(Vector<Real>(input_k.begin()+5, input_k.begin()+10), output_k_sub_d);
  ASSERT(IsApproximatelyEqual(output_k_sub_d, Vector<Real>(output_k_cmp.begin()+5, output_k_cmp.begin()+10), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[10]), output_k_cmp[10]));
  Vector<Real> output_k_sub_e(8);
  filter_k.Filter(Vector<Real>(input_k.begin()+11, input_k.begin()+19), output_k_sub_e);
  ASSERT(IsApproximatelyEqual(output_k_sub_e, Vector<Real>(output_k_cmp.begin()+11, output_k_cmp.begin()+19), VERY_SMALL));
  ASSERT(IsApproximatelyEqual(filter_k.Filter(input_k[19]), output_k_cmp[19]));

  // Testing slow update of filter
  DigitalFilter<Real> filter_t(UnaryVector<Real>(1.0));
  ASSERT(IsApproximatelyEqual(filter_t.Filter(0.76), 0.76));
  ASSERT(IsApproximatelyEqual(filter_t.Filter(1.0), 1.0));
  filter_t.SetNumeratorCoeffs(UnaryVector<Real>(0.3), 2);
  ASSERT(IsApproximatelyEqual(filter_t.Filter(1.0), 0.5*1.0+0.5*0.3));
  ASSERT(IsApproximatelyEqual(filter_t.Filter(1.0), 0.3));

  
  DigitalFilter<Real> filter_u(BinaryVector<Real>(1.0, 0.5));
  ASSERT(IsApproximatelyEqual(filter_u.Filter(0.76), 0.76));
  ASSERT(IsApproximatelyEqual(filter_u.Filter(1.0), 1.0+0.76*0.5));
  filter_u.SetNumeratorCoeffs(BinaryVector<Real>(2.0, -0.5), 3);
  ASSERT(IsApproximatelyEqual(filter_u.Filter(0.76), 1.18));
  ASSERT(IsApproximatelyEqual(filter_u.Filter(1.0), 1.54));
  ASSERT(IsApproximatelyEqual(filter_u.Filter(0.76), 1.02));
  ASSERT(IsApproximatelyEqual(filter_u.Filter(1.0), 1.62));

//  MaxGradientFilter filter_y(1.0);
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(0.0), 0.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.0), 1.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(0.0), 0.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-1.0), -1.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-3.0), -2.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-3.0), -3.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.5), -2.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.5), -1.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.5), 0.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.5), 1.0));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(1.5), 1.5));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-1.5), 0.5));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-1.5), -0.5));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-2.5), -1.5));
//  ASSERT(IsApproximatelyEqual(filter_y.Filter(-2.5), -2.5));


  return true;
}

inline void DigitalFilterSpeedTests() {

  constexpr size_t kernel_length_pow2 = 1024;
  constexpr size_t kernel_length_non_pow2 = 1023;
  constexpr size_t batch_size = 128;
  constexpr Real input_seconds = 10.0;
  constexpr Real input_single_seconds = input_seconds*10.0;
  constexpr Real sampling_frequency = 44100;
  
  const size_t num_input_samples = round(input_seconds*sampling_frequency);
  const size_t num_input_single_samples = round(input_single_seconds*sampling_frequency);
  
  RandomGenerator random_generator;

  Vector<Real> impulse_response = random_generator.Rand(kernel_length_pow2);
  DigitalFilter<Real> fir_filter(impulse_response);
  Vector<Real> input = random_generator.Rand(num_input_samples);
  Vector<Real> output(num_input_samples);
  Vector<Real> batch_input = random_generator.Rand(batch_size);
  Vector<Real> batch_output(batch_size);
  
  clock_t launch=clock();
  for (size_t i = 0; i<floor(num_input_samples/batch_size); i++) {
    fir_filter.Filter(batch_input, batch_output);
  }
  clock_t done=clock();

  std::cout<<"Fir filter speed (batch; filter length is not power of 2): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_seconds*100<<"% \n";

  launch=clock();
  for (size_t i=0; i<input.size(); ++i) {
    output[i] = fir_filter.Filter(input[i]);
  }
  done=clock();

  std::cout<<"Fir filter speed (sequential; filter length is not power of 2): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_seconds*100<<"% \n";



  DigitalFilter<Real> fir_filter_b(random_generator.Rand(kernel_length_non_pow2));

  launch=clock();
  for (size_t i = 0; i<floor(num_input_samples/batch_size); i++) {
    fir_filter_b.Filter(batch_input, batch_output);
  }
  done=clock();

  std::cout<<"Fir filter speed (batch; filter length is power of 2): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_seconds*100<<"% \n";

  launch=clock();
  for (size_t i=0; i<input.size(); ++i) {
    output[i] = fir_filter_b.Filter(input[i]);
  }
  done=clock();

  std::cout<<"Fir filter speed (sequential; filter length is power of 2): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_seconds*100<<"% \n";
  
  
  
  Vector<Real> input_single = random_generator.Rand(num_input_single_samples);
  Vector<Real> output_single(num_input_single_samples, 0.0);
  
  DigitalFilter<Real> fir_filter_c(UnaryVector<Real>(1.32));

  launch=clock();
  for (size_t i = 0; i<floor(num_input_single_samples/batch_size); i++) {
    fir_filter_c.Filter(batch_input, batch_output);
  }
  done=clock();

  std::cout<<"Fir filter speed (batch; filter length is 1): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_single_seconds*100<<"% \n";

  launch=clock();
  for (size_t i=0; i<input_single.size(); ++i) {
    output_single[i] = fir_filter_c.Filter(input_single[i]);
  }
  done=clock();

  std::cout<<"Fir filter speed (sequential; filter length is 1): "<<
  (done - launch) / ((Real) CLOCKS_PER_SEC) / input_single_seconds*100<<"% \n";
}

} // namespace mcl
