/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include <vector>

namespace mcl {

  

/**
 Performs the fft of the input signal.
 Equivalent to Matlab's fft(input, n_point)
 */
Vector<Complex> Fft(const Vector<Complex>& input,
                                 Int n_point) noexcept;

/**
 Performs the fft of the real input signal.
 Equivalent to Voice Box's rfft(input, n_point)
 */
Vector<Complex> Rfft(const Vector<Real>& input,
                                  Int n_point) noexcept;

  
/**
 Performs the fft of real vectors.
 Equivalent to Voice Box's rfft(input, n_point)
 */
std::vector<Vector<Complex> >
Rfft(const Vector<Vector<Real> >& input, Int n_point) noexcept;
  
  
/**
 Performs the inverse fft of conjugate symmetric spectrum.
 Equivalent to Voice Box's rfft(input, n_point)
 */
Vector<Real> Irfft(const Vector<Complex>& input,
                                Int n_point) noexcept;


/**
 Performs the inverse fft of conjugate symmetric spectra.
 Equivalent to Voice Box's rfft(input, n_point)
 */
Vector<Vector<Real> >
Irfft(const std::vector<Vector<Complex> >& input, Int n_point) noexcept;
  
/** 
 Performs the ifft of the input signal.
 Equivalent to Matlab's ifft(input, n_point)
 */
Vector<Complex> Ifft(
  const Vector<Complex>& input,
  Int n_point) noexcept;

/** 
 Performs the equivalent of Matlab's Hilbert (i.e. computes the so-called
 discrete-time analytic signal).
 */
Vector<Complex> Hilbert(const Vector<Real>& input) noexcept;

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

