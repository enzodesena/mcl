/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#include "iirfilter.h"
#include "comparisonop.h"
#include "vectorop.h"
#include "butter.h"
#include <limits>

namespace mcl {
  
Real IirFilter::GetNumeratorCoefficient(const Int coeff_id) const noexcept {
  ASSERT(coeff_id>=0 && coeff_id<=order());
  return B_[coeff_id]*A0_;
}
  
void IirFilter::SetNumeratorCoefficient(const Int coeff_id,
                                          const Real value) noexcept {
  ASSERT(coeff_id >= 0 && coeff_id<=order());
  B_[coeff_id] = value/A0_;
}
  
Vector<Real> IirFilter::B() const {
  // Return the non-normalised version
  return Multiply(B_, A0_);
}

Vector<Real> IirFilter::A() const {
  // Return the non-normalised version
  return Multiply(A_, A0_);
}

// Constructor
IirFilter::IirFilter() :
  B_(mcl::UnaryVector<Real>(1.0)), A_(mcl::UnaryVector<Real>(1.0)) {
  Int size = B_.length();
  state_ = new Real[size];
  for (Int i=0; i<size; ++i) { state_[i] = 0.0; }
}
  
IirFilter::IirFilter(Vector<Real> B, Vector<Real> A) : 
          B_(B), A_(A) {
  // TODO: implement also for B.size != A.size
  ASSERT(B.length() == A.length());
  ASSERT(B.length() >= 1);
  
  A0_ = A[0];
  if (! IsEqual(A[0], 1.0, std::numeric_limits<Real>::epsilon())) {
    B_ = Multiply(B, (Real) 1.0 / A[0]);
    A_ = Multiply(A, (Real) 1.0 / A[0]);
  }
  
  
  Int size = B.length();
  state_ = new Real[size];
  for (Int i=0; i<size; ++i) { state_[i] = 0.0; }
}

// Copy constructor
IirFilter::IirFilter(const IirFilter& copy) :
          B_(copy.B_), A_(copy.A_), A0_(copy.A0_) {
  Int size = B_.length();
  state_ = new Real[size];
  for (Int i=0; i<size; ++i) {
    state_[i] = copy.state_[i];
  }
}

// Assignment constructor
IirFilter& IirFilter::operator= (const IirFilter& other) {
  if (this != &other) {
    Int size = other.B_.length();
    delete [] state_;
    state_ = new Real[size];
    for (Int i=0; i<size; ++i) { state_[i] = other.state_[i]; }
    B_ = other.B_;
    A_ = other.A_;
    A0_ = other.A0_;
  }
  return *this;
}

IirFilter::~IirFilter() { delete [] state_; }

Int IirFilter::order() const noexcept {
  return Max(B_.length(), A_.length())-1;
}

Real IirFilter::GetDenominatorCoefficient(const Int coeff_id) const noexcept {
  return A_.at(coeff_id);
}
  
void IirFilter::SetCoefficients(const Vector<Real>& B,
                                 const Vector<Real>& A) noexcept {
  // TODO: implement case where length changes.
  ASSERT(B_.length() == B.length());
  ASSERT(A_.length() == A.length());
  
  B_ = B;
  A_ = A;
  A0_ = A[0];
}
  
void IirFilter::SetCoefficients(const IirFilter& other_filter) noexcept {
  const Int filter_order = order();
  assert(filter_order == other_filter.order());
  
  for (Int i=0; i<=filter_order; ++i) {
    SetNumeratorCoefficient(i, other_filter.GetNumeratorCoefficient(i));
    SetDenominatorCoefficient(i, other_filter.GetDenominatorCoefficient(i));
  }
}
  
  
void IirFilter::SetDenominatorCoefficient(const Int coeff_id,
                                            const Real value) noexcept {
  ASSERT(coeff_id >= 0 &&coeff_id<(Int)A_.length());
  A_[coeff_id] = value;
  if (coeff_id == 0) { A0_ = value; }
}

void IirFilter::Filter(const Real* input_data, const Int num_samples,
                       Real* output_data) noexcept {
  if (B_.length() == 1) {
    Multiply(input_data, num_samples, B_[0], output_data);
  } else {
    FilterSerial(input_data, num_samples, output_data);
  }
}
  
  
Real IirFilter::Filter(Real input) noexcept {
  Real v = input; // The temporary value in the recursive branch.
  Real output(0.0);
  
  Int size = B_.length();
  
  // Speed up return for simple gain filters
  if (size == 1) { return v*B_[0]; }
  
  // The index i in both loops refers to the branch in the classic plot of a 
  // direct form II, with the highest branch (the one multiplied by b(0) only)
  // being i=0.
  for (Int i=1; i<size; ++i) { 
    v += state_[i-1]*(-A_[i]);
    output += state_[i-1]*B_[i];
  }
  
  for (Int i=(size-1); i>=1; --i) {
    state_[i] = state_[i-1];
  }
  
  state_[0] = v;
  
  output += v*B_[0];
  
  return output;
}
  
void IirFilter::Reset() {
  const Int size = B_.length();
  for (Int i=0; i<size; ++i) { state_[i] = 0.0; }
}



IirFilter GainFilter(Real gain) {
  Vector<Real> B(1);
  Vector<Real> A(1);
  B[0] = 1.0*gain;
  A[0] = 1.0;
  
  return IirFilter(B,A);
}

IirFilter IdenticalFilter() { return GainFilter(1.0); }

IirFilter WallFilter(WallType wall_type, Real sampling_frequency) {
  // TODO: implement for frequencies other than 44100
  if (! IsEqual(sampling_frequency, 44100)) {
    mcl::Logger::GetInstance().LogError("Attempting to use a wall filter "
        "designed for 44100 Hz sampling frequency with a sampling frequency "
        "of %f Hz. The filter response will be inaccurate.", sampling_frequency);
  }
  
  Vector<Real> B;
  Vector<Real> A;
  
  switch (wall_type) {
    case kRigid: {
      B.push_back(1.0);
      A.push_back(1.0);
      break;
    }
    case kCarpetPile: {
      // B_carpet_pile=[0.562666833756030  -1.032627191365576   0.469961155406544];
      // A_carpet_pile=[1.000000000000000  -1.896102349247713   0.896352947528892];
      B.push_back(0.562666833756030);
      B.push_back(-1.032627191365576);
      B.push_back(0.469961155406544);
      A.push_back(1.000000000000000);
      A.push_back(-1.896102349247713);
      A.push_back(0.896352947528892);
      break;
    }
    case kCarpetCotton: {
      // B_carpet_cotton = [0.687580695329600  -1.920746652319969   1.789915765926473  -0.556749690855965];
      // A_carpet_cotton = [1.000000000000000  -2.761840732459190   2.536820778736938  -0.774942833868750];
      B.push_back(0.687580695329600);
      B.push_back(-1.920746652319969);
      B.push_back(1.789915765926473);
      B.push_back(-0.556749690855965);
      A.push_back(1.000000000000000);
      A.push_back(-2.761840732459190);
      A.push_back(2.536820778736938);
      A.push_back(-0.774942833868750);
      break;
    }
    case kWallBricks: {
      // B_wall=[0.978495798553620  -1.817487798457697   0.839209660516074];
      // A_wall=[1.000000000000000  -1.858806492488240   0.859035906864860];
      B.push_back(0.978495798553620);
      B.push_back(-1.817487798457697);
      B.push_back(0.839209660516074);
      A.push_back(1.000000000000000);
      A.push_back(-1.858806492488240);
      A.push_back(0.859035906864860);
      break;
    }
    case kCeilingTile: {
      // B_ceiling=[0.168413736374283  -0.243270224986791   0.074863520490536];
      // A_ceiling=[1.000000000000000  -1.845049094190385   0.845565720138466];
      B.push_back(0.168413736374283);
      B.push_back(-0.243270224986791);
      B.push_back(0.074863520490536);
      A.push_back(1.000000000000000);
      A.push_back(-1.845049094190385);
      A.push_back(0.845565720138466);
      break;
    }
    default: {
      mcl::Logger::GetInstance().LogError("Unrecognised type of wall filter. Reverting to a completely absorptive one.");
      B.push_back(0.0);
      A.push_back(1.0);
      break;
    }
  }
  
  return IirFilter(B,A);
}

IirFilter PinkifierFilter() {
  Vector<Real> poles(5);
  poles[0] = 0.9986823;
  poles[1] = 0.9914651;
  poles[2] = 0.9580812;
  poles[3] = 0.8090598;
  poles[4] = 0.2896591;
  Vector<Real> zeros(5);
  zeros[0] = 0.9963594;
  zeros[1] = 0.9808756;
  zeros[2] = 0.9097290;
  zeros[3] = 0.6128445;
  zeros[4] = -0.0324723;
  
  std::vector<Complex> num = Poly(zeros);
  std::vector<Complex> den = Poly(poles);
  
  return IirFilter(RealPart(num),RealPart(den));
}

  
Vector<Real> IirFilterBank::Filter(const Real input) {
  const Int N = filters_.length();
  Vector<Real> outputs(N);
  for (Int i=0; i<N; ++i) {
    outputs[i] = filters_[i].Filter(input);
  }
  return outputs;
}

  
std::vector<Vector<Real> >
IirFilterBank::Filter(const Vector<Real>& input) {
  const Int N = filters_.length();
  std::vector<Vector<Real> > outputs(N);
  for (Int i=0; i<N; ++i) {
    outputs[i] = filters_[i].Filter(input);
  }
  return outputs;
}

  
void IirFilterBank::Reset() {
  const Int N = filters_.length();
  for (Int i=0; i<N; ++i) {
    filters_[i].Reset();
  }
}

  
} // namespace mcl




