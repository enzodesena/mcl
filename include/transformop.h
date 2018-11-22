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
template<typename T>
inline Vector<Complex<T>> Fft(
  const Vector<Complex<T>>& input,
  size_t n_point) noexcept;

/**
 Performs the fft of the real input signal.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Complex<T>> Rfft(
  const Vector<T>& input,
  size_t n_point) noexcept;

  
/**
 Performs the fft of real vectors.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Vector<Complex<T>>> Rfft(
  const Vector<Vector<T> >& input,
  size_t n_point) noexcept;
  
  
/**
 Performs the inverse fft of conjugate symmetric spectrum.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<T> Irfft(
  const Vector<Complex<T>>& input,
  size_t n_point) noexcept;


/**
 Performs the inverse fft of conjugate symmetric spectra.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Vector<T>> Irfft(
  const std::vector<Vector<Complex<T>> >& input,
  size_t n_point) noexcept;
  
/** 
 Performs the ifft of the input signal.
 Equivalent to Matlab's ifft(input, n_point)
 */
template<typename T>
inline Vector<Complex<T>> Ifft(
  const Vector<Complex<T>>& input,
  Int n_point) noexcept;

/** 
 Performs the equivalent of Matlab's Hilbert (i.e. computes the so-called
 discrete-time analytic signal).
 */
template<typename T>
inline Vector<Complex<T>> Hilbert(
  const Vector<T>& input) noexcept;

/** 
 Returns the real cepstrum of the real sequence X.
 Equivalent to Matlab's rceps(vector)
 */
template<typename T>
inline Vector<T> RCeps(
  const Vector<T>& vector) noexcept;

/** 
 Returns the (unique) minimum-phase sequence that has the same real
 cepstrum as vector. Equivalent to Matlab's [~, out] = rceps(vector).
 */
template<typename T>
inline Vector<T> MinPhase(
  const Vector<T>& vector) noexcept;

  
/** Equivalent to Matlab's xcorr(vect_a, vect_b) */
template<typename T>
inline Vector<T> XCorr(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b);
// The method XCorr naturally is placed in vectorop, but since it depends
// on the Fft method, I place it here, so someone who doesn't want to
// compile with KissFFT won't get compilation errors.

bool TransformOpTest();

} /**< namespace mcl */

