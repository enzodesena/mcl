/*
 firfilter.h
 MCL
 Copyright (c) 2012, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@me.com
 
 */

#ifndef MCL_FIRFILTER_H
#define MCL_FIRFILTER_H

#include <vector>

#include "mcltypes.h"
#include "digitalfilter.h"

namespace mcl {
/** FIR Filter */
class FirFilter : public DigitalFilter {
public:
  /** Constructs a default FIR filter, i.e. identical filter */
  FirFilter();
  
  /** Constructs an FIR filter with impulse response B. */
  FirFilter(std::vector<Real> B);
  
  /** 
   Returns the output of the filter for an input equal to `input`.
   For example, if B=1, A=1, output will be equal to input. 
   As a second example, if B=[0,1], A=[1], you will have 
   (1) Filter(0.5)==0 and then
   (2) Filter(0.0)==0.5
   */
  virtual Real Filter(const Real input);
  
  /** Returns the output of the filter for an input signal equal to `input`. */
  virtual std::vector<Real> Filter(const std::vector<Real>& input);
  
  /** 
   Updates the filter coefficients. May cause articafts if the coefficients are
   updated too rapidly. If the `impulse_response` changes length than the
   filter is resetted to zero state.
   */
  void UpdateFilter(std::vector<Real> impulse_response, Int update_length = 1);
  
  /** Resets the state of the filter */
  void Reset();
  
  /** Returns the impulse response of the filter */
  std::vector<Real> impulse_response();
  
  /** Constructs a filter for which output==gain*input always. */
  static FirFilter GainFilter(const Real gain);
  
  /** Tests */
  static bool Test();
  
  static void SpeedTests();
  
  virtual ~FirFilter() {}
  
private:
  std::vector<Real> FilterSequential(const std::vector<Real>& input);
  /** Method called to slowly update the filter coefficients. It is called
   every time one of the Filter method is called and is activated only
   if updating_ = true. TODO: uniformise action between sequential and
   batch. */
  void UpdateCoefficients();
  
  std::vector<float> impulse_response_;
  std::vector<float> impulse_response_old_;
  Int update_index_;
  Int update_length_;
  
  std::vector<float> coefficients_;
  std::vector<float> delay_line_;
  Int counter_;
  Int length_;
};
  
} // namespace mcl


#endif
