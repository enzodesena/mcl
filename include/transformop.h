/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
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
std::vector<Complex> Fft(const std::vector<Complex>& input,
                                 Int n_point) noexcept;

/**
 Performs the fft of the real input signal.
 Equivalent to Voice Box's rfft(input, n_point)
 */
std::vector<Complex> Rfft(const Vector<Real>& input,
                                  Int n_point) noexcept;

  
/**
 Performs the fft of real vectors.
 Equivalent to Voice Box's rfft(input, n_point)
 */
std::vector<std::vector<Complex> >
Rfft(const std::vector<Vector<Real> >& input, Int n_point) noexcept;
  
  
/**
 Performs the inverse fft of conjugate symmetric spectrum.
 Equivalent to Voice Box's rfft(input, n_point)
 */
Vector<Real> Irfft(const std::vector<Complex>& input,
                                Int n_point) noexcept;


/**
 Performs the inverse fft of conjugate symmetric spectra.
 Equivalent to Voice Box's rfft(input, n_point)
 */
std::vector<Vector<Real> >
Irfft(const std::vector<std::vector<Complex> >& input, Int n_point) noexcept;
  
/** 
 Performs the ifft of the input signal.
 Equivalent to Matlab's ifft(input, n_point)
 */
std::vector<Complex> Ifft(
  const std::vector<Complex>& input,
  Int n_point) noexcept;

/** 
 Performs the equivalent of Matlab's Hilbert (i.e. computes the so-called
 discrete-time analytic signal).
 */
std::vector<Complex> Hilbert(const Vector<Real>& input) noexcept;

/** 
 Returns the real cepstrum of the real sequence X.
 Equivalent to Matlab's rceps(vector)
 */
Vector<Real> RCeps(const Vector<Real>& vector) noexcept;

/** 
 Returns the (unique) minimum-phase sequence that has the same real
 cepstrum as vector. Equivalent to Matlab's [~, out] = rceps(vector).
 */
Vector<Real> MinPhase(const Vector<Real>& vector) noexcept;

  
/** Equivalent to Matlab's xcorr(vect_a, vect_b) */
Vector<Real> XCorr(const Vector<Real>& vector_a,
                                const Vector<Real>& vector_b);
// The method XCorr naturally is placed in vectorop, but since it depends
// on the Fft method, I place it here, so someone who doesn't want to
// compile with KissFFT won't get compilation errors.

bool TransformOpTest();

} /**< namespace mcl */

#endif
