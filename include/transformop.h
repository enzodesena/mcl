/*
 MCL
 Copyright (c) 2012-18, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */

#pragma once

#include "mcltypes.h"
#include "vector.h"
#include "elementaryop.h"
#include "kissfft.hh"

namespace mcl {

// Forward declarations
template<typename T>
bool IsApproximatelyReal(
  const Vector<T>& vector,
  const T precision) noexcept;
// End of forward declarations

/**
 Performs the fft of the input signal.
 Equivalent to Matlab's fft(input, n_point)
 */
template<typename T>
inline Vector<Complex<T>> Fft(
  const Vector<Complex<T>>& input,
  size_t n_point) noexcept
{
  Vector<Complex<T>> padded(n_point);
  ZeroPad(input, padded);
  kissfft<T> fft((int) n_point, false);
  Vector<Complex<T>> outbuf(n_point);
  fft.transform(&padded[0], &outbuf[0]);
  return outbuf;
}

/**
 Performs the ifft of the input signal.
 Equivalent to Matlab's ifft(input, n_point)
 */
template<typename T>
inline Vector<Complex<T>> Ifft(
  const Vector<Complex<T>>& input,
  Int n_point) noexcept
{
  Vector<Complex<T>> padded(n_point);
  ZeroPad(input, padded);
  kissfft<T> fft((int) n_point, true);
  Vector<Complex<T>> outbuf(n_point);
  
  fft.transform(&padded[0], &outbuf[0]);
  
  return Multiply(outbuf, Complex<T>(T(1.0)/n_point, T(0.0)));
}

/**
 Performs the fft of the real input signal.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Complex<T>> Rfft(
  const Vector<T>& input,
  size_t n_point) noexcept
{
  return Elements(
    Fft(CastToComplex(input), n_point), 0,
    (Int) floor(1.0+((double)n_point)/2.0)-1);
}

  
/**
 Performs the fft of real vectors.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Vector<Complex<T>>> Rfft(
  const Vector<Vector<T> >& input,
  size_t n_point) noexcept
{
  Vector<Vector<Complex<T>>> outputs(input.length());
  for (size_t i=0; i<input.length(); ++i) {
    outputs[i] = Rfft(input[i], n_point);
  }
  return outputs;
}
  
  
  
/**
 Performs the inverse fft of conjugate symmetric spectrum.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<T> Irfft(
  const Vector<Complex<T>>& input,
  size_t n_point) noexcept
{
// If n_point is even, then it includes the Nyquist term (the only
  // non-repeated term, together with DC) and is of dimension M=N/2 + 1
  // whereas if N is odd then there is no Nyquist term
  // and the input is of dimension M=(N+1)/2.
  Int M = (Rem(n_point, (size_t) 2) == 0) ? (n_point/2 + 1) : (n_point+1)/2;
  Vector<Complex<T>> zero_padded(M);
  ZeroPad(input, zero_padded);
  Vector<Complex<T>> spectrum;
  if (Rem(n_point, (size_t) 2) == 0)
  { // If n_point is even
    spectrum = Concatenate(zero_padded, Flip(Conj(
      Elements(zero_padded, 1, n_point/2-1))));
  }
  else
  { // If n_point is odd
    spectrum = Concatenate(zero_padded,
      Flip(Conj(Elements(zero_padded, 1, (n_point+1)/2-1))));
  }
  ASSERT(spectrum.length() == n_point);
  Vector<Complex<T>> output = Ifft(spectrum, n_point);
  ASSERT(IsApproximatelyReal(output));
  return RealPart(output);
}


/**
 Performs the inverse fft of conjugate symmetric spectra.
 Equivalent to Voice Box's rfft(input, n_point)
 */
template<typename T>
inline Vector<Vector<T>> Irfft(
  const std::vector<Vector<Complex<T>> >& input,
  size_t n_point) noexcept
{
  Vector<Vector<T>> outputs;
  for (Int i=0; i<(Int)input.length(); ++i) {
    outputs.PushBack(Irfft(input[i], n_point));
  }
  return outputs;
}
  

/** 
 Performs the equivalent of Matlab's Hilbert (i.e. computes the so-called
 discrete-time analytic signal).
 */
template<typename T>
inline Vector<Complex<T>> Hilbert(
  const Vector<T>& input) noexcept
{
  Int n = input.length();
  Vector<Complex<T>> x = Fft(ComplexVector(input), n);
  Vector<Complex<T>> h = Zeros<Complex<T>>(n);
  
  if (n > 0 && 2*floor(n/2) == n)
  {
    // even and nonempty
    h[0] = 1.0;
    h[n/2] = 1.0;
    for (Int i=1; i<(n/2); ++i) {
      h[i] = 2.0;
    }
  }
  else if (n>0)
  {
    // odd and nonempty
    h[0] = 1.0;
    for (Int i=1; i<((n+1)/2); ++i) {
      h[i] = 2.0;
    }
  }
  
  // x = ifft(x.*h(:,ones(1,size(x,2))));
  return Ifft(Multiply(x, h), n);
}

/** 
 Returns the real cepstrum of the real sequence X.
 Equivalent to Matlab's rceps(vector)
 */
template<typename T>
inline Vector<T> RCeps(
  const Vector<T>& x) noexcept
{
  Int n = Length(x);
  Vector<T> fftxabs = Abs(Fft(CastToComplex<T>(x), n));
  //  TODO: implement this check.
  //  if any(fftxabs == 0)
  //    error(generatemsgid('ZeroInFFT'),...
  //          'The Fourier transform of X contains zeros. Therefore, the T cepstrum of X does not exist.');
  //  end
  return RealPart(Ifft(ComplexVector(Log(fftxabs)), n));
}

/** 
 Returns the (unique) minimum-phase sequence that has the same real
 cepstrum as vector. Equivalent to Matlab's [~, out] = rceps(vector).
 */
template<typename T>
inline Vector<T> MinPhase(
  const Vector<T>& x) noexcept
{
  Int n = Length(x);
  Int odd = Fix(Rem((Int) n, (Int) 2));
  Vector<T> wn_1 = Concatenate(
    UnaryVector((T) 1.0),
    Multiply(Ones<T>((n+odd)/2-1), (T) 2.0));
  Vector<T> wn_2 = Concatenate(wn_1, Ones<T>(1-(Int)Rem((Int)n,(Int)2)));
  Vector<T> wn = Concatenate(wn_2, Zeros<T>((n+odd)/2-1));
  return RealPart(Ifft(Exp(Fft(ComplexVector(Multiply(wn, RCeps(x))),n)),n));
}

  
// The method XCorr naturally is placed in vectorop, but since it depends
// on the Fft method, I place it here, so someone who doesn't want to
// compile with KissFFT won't get compilation errors.
/** Equivalent to Matlab's xcorr(vect_a, vect_b) */
template<typename T>
inline Vector<T> XCorr(
  const Vector<T>& vector_a,
  const Vector<T>& vector_b)
{
  // TODO: implement for different sizes
  ASSERT(vector_a.length() == vector_b.length());
  Int M = (Int)vector_a.length();
  Int n_fft = (UInt) pow(2.0, NextPow2(2*M-1));
  Vector<Complex<T>> x = Fft(ComplexVector(vector_a), n_fft);
  Vector<Complex<T>> y = Fft(ComplexVector(vector_b), n_fft);
  Vector<Complex<T>> c = Ifft(Multiply(x, Conj(y)), n_fft);
  
  // Ignore residual imaginary part
  Vector<T> c_T = RealPart(c);
  Vector<T> output = Zeros<T>(2*M-1);
  Int end = c_T.length();
  Int maxlag = M-1;
  Int k = 0; // running index
  for (Int i=end-maxlag+1-1; i<=(end-1); ++i)
  {
    output[k++] = c_T[i];
  }
  for (Int i=1-1; i<=(maxlag+1-1); ++i)
  {
    output[k++] = c_T[i];
  }
  return output;
}

bool TransformOpTest();

} /**< namespace mcl */

