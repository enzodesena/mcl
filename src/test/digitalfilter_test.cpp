/*
 iirfilter_test.cpp
 Matlab Cpp Library (MCL)
 Copyright (c) 2012, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@me.com
 
 Last committed:     $Revision: 93 $
 Last changed date:  $Date: 2012-06-07 09:58:54 +0100 (Thu, 07 Jun 2012) $
 */

#include "iirfilter.h"
#include "firfilter.h"
#include <vector>
#include "mcltypes.h"
#include "mcl.h"

namespace mcl {

bool IirFilter::Test() {
  IirFilter filter_a = IirFilter::IdenticalFilter();
  assert(IsEqual(filter_a.Filter(1.2), 1.2));
  
  IirFilter filter_b = IirFilter::GainFilter(0.76);
  assert(IsEqual(filter_b.Filter(1.2), 0.912));
  
  std::vector<Real> B;
  std::vector<Real> A;
  B.push_back(0.76);
  B.push_back(-3.06);
  B.push_back(1.76);
  B.push_back(1.76);
  A.push_back(1.0);
  A.push_back(-0.5);
  A.push_back(0.23);
  A.push_back(0.75);
  
  IirFilter filter_c(B,A);
  
  assert(IsEqual(filter_c.Filter(1.2), 0.912));
  assert(IsEqual(filter_c.Filter(0.2), -3.064));
  assert(IsEqual(filter_c.Filter(0.5), 0.13824));
  assert(IsEqual(filter_c.Filter(-1.2), 0.11184));
  assert(IsEqual(filter_c.Filter(0.0), 7.226124799999999999999999999999));
  assert(IsEqual(filter_c.Filter(0.0), 2.251659199999999999999999999999));
  
  IirFilter filter_d = IirFilter::WallFilter(carpet_pile, 44100);
  assert(IsEqual(filter_d.Filter(1.0), 0.562666833756030));
  assert(IsEqual(filter_d.Filter(0.5), 0.315580130841020));
  
  IirFilter filter_e = IirFilter::WallFilter(carpet_cotton, 44100);
  assert(IsEqual(filter_e.Filter(1.0), 0.687580695329600));
  assert(IsEqual(filter_e.Filter(0.5), 0.322032066558733));
  
  IirFilter filter_f = IirFilter::WallFilter(wall_bricks, 44100);
  assert(IsEqual(filter_f.Filter(1.0), 0.978495798553620));
  assert(IsEqual(filter_f.Filter(0.5), 0.490594444043047));
  
  IirFilter filter_g = IirFilter::WallFilter(ceiling_tile, 44100);
  assert(IsEqual(filter_g.Filter(1.0), 0.168413736374283));
  assert(IsEqual(filter_g.Filter(0.5), 0.151668254946940));
  
  return true;
}
  
  
bool FirFilter::Test() {
  using mcl::IsEqual;
  
  std::vector<Real> impulse_resp(3);
  impulse_resp[0] = 0.2;
  impulse_resp[1] = -0.1;
  impulse_resp[2] = 2.5;
  
  std::vector<Real> impulse(3);
  impulse[0] = 1.0;
  impulse[1] = 0.0;
  impulse[2] = 0.0;
  
  FirFilter filter_a(impulse_resp);
  assert(IsEqual(filter_a.Filter(impulse), impulse_resp));
  
  FirFilter filter_b(impulse_resp);
  std::vector<Real> input_a(4);
  input_a[0] = 0.6;
  input_a[1] = -3.5;
  input_a[2] = 5.6;
  input_a[3] = 2.3;
  
  std::vector<Real> output_a_cmp(4);
  output_a_cmp[0] = 0.1200;
  output_a_cmp[1] = -0.7600;
  output_a_cmp[2] = 2.9700;
  output_a_cmp[3] = -8.8500;
  
  assert(IsEqual(filter_b.Filter(input_a), output_a_cmp));
  
  FirFilter filter_c(impulse_resp);
  assert(IsEqual(filter_c.Filter(0.6), 0.1200));
  assert(IsEqual(filter_c.Filter(-3.5), -0.7600));
  
  // Testing copy constructor
  FirFilter filter_d = filter_c;
  assert(IsEqual(filter_d.Filter(5.6), 2.9700));
  assert(IsEqual(filter_d.Filter(2.3), -8.8500));
  
  assert(IsEqual(filter_c.Filter(5.6), 2.9700));
  assert(IsEqual(filter_c.Filter(2.3), -8.8500));
  
  
  return true;
}


} // namespace mcl
