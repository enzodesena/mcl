/*
 transformop.h
 MCL
 
 This contains transforming operations (fft, ifft, hilbert...)
 
 Authors: Enzo De Sena, enzodesena@me.com
 
 */

#ifndef MCL_TRANSFORMOP_H
#define MCL_TRANSFORMOP_H

#include "mcltypes.h"
#include <vector>


namespace mcl {

  

/**
 Performs the fft of the input signal.
 Equivalent to Matlab's fft(input, n_point)
 */
std::vector<Complex> Fft(const std::vector<Complex>& input, UInt n_point);

/** 
 Performs the ifft of the input signal.
 Equivalent to Matlab's ifft(input, n_point)
 */
std::vector<Complex> Ifft(const std::vector<Complex>& input, UInt n_point);

/** 
 Performs the equivalent of Matlab's Hilbert (i.e. computes the so-called
 discrete-time analytic signal).
 */
std::vector<Complex> Hilbert(const std::vector<Real>& input);

/** 
 Returns the real cepstrum of the real sequence X.
 Equivalent to Matlab's rceps(vector)
 */
std::vector<Real> RCeps(const std::vector<Real>& vector);

/** 
 Returns the (unique) minimum-phase sequence that has the same real
 cepstrum as vector. Equivalent to Matlab's [~, out] = rceps(vector).
 */
std::vector<Real> MinPhase(const std::vector<Real>& vector);


bool TransformOpTest();

} /**< namespace mcl */

#endif
