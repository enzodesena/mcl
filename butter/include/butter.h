/*
 MCL
 
 Authors: user261002 @ stackoverflow.com
 Authors: Enzo De Sena, enzodesena@gmail.com
 
 Adapted from "http://stackoverflow.com/questions/10373184/bandpass-butterworth-filter-implementation-in-c"
 */


#ifndef MCL_BUTTER_H
#define MCL_BUTTER_H

#include <vector>
#include "iirfilter.h"

namespace mcl {
  
template<typename T>
inline Vector<T> ComputeLP(
  int FilterOrder) noexcept
{
  int m;
  int i;
  Vector<T> NumCoeffs(FilterOrder+1);
  NumCoeffs[0] = 1;
  NumCoeffs[1] = FilterOrder;
  m = FilterOrder/2;
  for( i=2; i <= m; ++i) {
    NumCoeffs[i] =(T) (FilterOrder-i+1)*NumCoeffs[i-1]/i;
    NumCoeffs[FilterOrder-i]= NumCoeffs[i];
  }
  NumCoeffs[FilterOrder-1] = FilterOrder;
  NumCoeffs[FilterOrder] = 1;
  return NumCoeffs;
}

template<typename T>
inline Vector<T> ComputeHP(
  int FilterOrder) noexcept
{
  int i;
  Vector<T> NumCoeffs = ComputeLP<T>(FilterOrder);
  //if(NumCoeffs == NULL ) return( NULL );
  for( i = 0; i <= FilterOrder; ++i)
  {
    if( i % 2 )
    {
      NumCoeffs[i] = -NumCoeffs[i];
    }
  }
  return NumCoeffs;
}
  
template<typename T>
inline Vector<T> ComputeNumCoeffs(
  int FilterOrder,
  T Lcutoff,
  T Ucutoff,
  Vector<T> DenC) noexcept
{
  ASSERT(std::isgreaterequal(Lcutoff, 0.0) && std::isless(Ucutoff, 1.0));
  int i;
  Vector<T> NumCoeffs(2*FilterOrder+1);
  Vector<Complex<T> > NormalizedKernel(2*FilterOrder+1);
  // if( NormalizedKernel == NULL ) return( NULL );
  Vector<T> TCoeffs = ComputeHP<T>(FilterOrder);
  // if( TCoeffs == NULL ) return( NULL );
  for( i = 0; i < FilterOrder; ++i)
  {
    NumCoeffs[2*i] = TCoeffs[i];
    NumCoeffs[2*i+1] = 0.0;
  }
  NumCoeffs[2*FilterOrder] = TCoeffs[FilterOrder];
  Vector<T> cp(2);
  T Bw, Wn;
  cp[0] = 2*2.0*tan(PI * Lcutoff/ 2.0);
  cp[1] = 2*2.0*tan(PI * Ucutoff / 2.0);

  Bw = cp[1] - cp[0];
  //center frequency
  Wn = sqrt(cp[0]*cp[1]);
  Wn = 2*atan2(Wn,4);

  const std::complex<T> result = std::complex<T>(-1,0);

  for(int k = 0; k<(2*FilterOrder+1); k++)
  {
    NormalizedKernel[k] = std::exp(-sqrt(result)*Wn*((T) k));
  }
  T b=0;
  T den=0;
  for(int d = 0; d<(2*FilterOrder+1); d++)
  {
    b+=real(NormalizedKernel[d]*NumCoeffs[d]);
    den+=real(NormalizedKernel[d]*DenC[d]);
  }
  for(int c = 0; c<(2*FilterOrder+1); c++)
  {
    NumCoeffs[c]=(NumCoeffs[c]*den)/b;
  }
  return NumCoeffs;
}

template<typename T>
inline Vector<T> ComputeDenCoeffs(
  int FilterOrder,
  T Lcutoff,
  T Ucutoff) noexcept
{
  if ((Lcutoff < 0.0) | (Ucutoff > 1.0)) { ASSERT(false); }
  
  int k;            // loop variables
  T theta;     // PI * (Ucutoff - Lcutoff) / 2.0
  T cp;        // cosine of phi
  T st;        // sine of theta
  T ct;        // cosine of theta
  T s2t;       // sine of 2*theta
  T c2t;       // cosine 0f 2*theta
  T PoleAngle;      // pole angle
  T SinPoleAngle;     // sine of pole angle
  T CosPoleAngle;     // cosine of pole angle
  T a;         // workspace variables
  
  cp = cos(PI * (Ucutoff + Lcutoff) / 2.0);
  theta = PI * (Ucutoff - Lcutoff) / 2.0;
  st = sin(theta);
  ct = cos(theta);
  s2t = 2.0*st*ct;        // sine of 2*theta
  c2t = 2.0*ct*ct - 1.0;  // cosine of 2*theta
  
  Vector<T> RCoeffs(2 * FilterOrder); // z^-2 coefficients
  Vector<T> TCoeffs(2 * FilterOrder); // z^-1 coefficients
  
  for( k = 0; k < FilterOrder; ++k )
  {
    PoleAngle = PI * (T)(2*k+1)/(T)(2*FilterOrder);
    SinPoleAngle = sin(PoleAngle);
    CosPoleAngle = cos(PoleAngle);
    a = 1.0 + s2t*SinPoleAngle;
    RCoeffs[2*k] = c2t/a;
    RCoeffs[2*k+1] = s2t*CosPoleAngle/a;
    TCoeffs[2*k] = -2.0*cp*(ct+st*SinPoleAngle)/a;
    TCoeffs[2*k+1] = -2.0*cp*st*CosPoleAngle/a;
  }
  
  Vector<T> DenomCoeffsTemp = TrinomialMultiply(FilterOrder, TCoeffs, RCoeffs);
  
  Vector<T> DenomCoeffs(2*FilterOrder+1);
  
  DenomCoeffs[0] = 1.0;
  DenomCoeffs[1] = DenomCoeffsTemp[0];
  DenomCoeffs[2] = DenomCoeffsTemp[2];
  for( k = 3; k <= 2*FilterOrder; ++k )
  {
    DenomCoeffs[k] = DenomCoeffsTemp[2*k-2];
  }
  
  return DenomCoeffs;
}
  
/**
 Construncts a bandpass butterworth filter. Equivalent to Matlab's
 butter(order, [w_low, w_high])
 */
template<typename T>
inline IirFilter<T> Butter(
  const Int order,
  const Real w_low,
  const Real w_high) noexcept
{
  Vector<T> DenC = ComputeDenCoeffs((int) order, w_low, w_high);
  Vector<T> NumC = ComputeNumCoeffs((int) order, w_low, w_high, DenC);
  Vector<T> denominator(DenC.begin(), DenC.end());
  Vector<T> numerator(NumC.begin(), NumC.end());
  return IirFilter<T>(numerator, denominator);
}

/** Constructs a digital octave filter with given center frequency
 and sampling frequency.
 */
template<typename T>
inline IirFilter<T> OctaveFilter(
  const Int order,
  const T center_frequency,
  const T sampling_frequency) noexcept
{
  // beta = pi/2/N/sin(pi/2/N);
  T beta = PI/2.0/((T) order)/sin(PI/2.0/((T) order));
  
  // alpha = (1+sqrt(1+8*beta^2))/4/beta;
  T alpha = (1.0+sqrt(1.0+8.0*pow(beta,2.0)))/4.0/beta;
  
  // W1 = Fc/(Fs/2)*sqrt(1/2)/alpha;
  T W1 = center_frequency/(sampling_frequency/2.0)*sqrt(1.0/2.0)/alpha;
  
  // W2 = Fc/(Fs/2)*sqrt(2)*alpha;
  T W2 = center_frequency/(sampling_frequency/2.0)*sqrt(2.0)*alpha;
  
  // [B,A] = butter(N,[W1,W2]);
  return Butter<T>(order, W1, W2);
}

template<typename T>
inline IirFilterBank<T> OctaveFilterBank(
  const Int order,
  const Int num_bands,
  const T starting_frequency,
  const T sampling_frequency) noexcept
{
  T current_frequency = starting_frequency;
  Vector<IirFilter<T>> filters;
  for (Int i=0; i<num_bands; ++i)
  {
    mcl::IirFilter
    octave_filter = mcl::OctaveFilter<T>(
      order,
      current_frequency,
      sampling_frequency);
    filters.PushBack(octave_filter);
    current_frequency = current_frequency * 2.0;
  }
  return IirFilterBank<T>(filters);
}

template<typename T>
inline Vector<T> TrinomialMultiply(
  int FilterOrder,
  Vector<T> b,
  Vector<T> c) noexcept
{
  int i, j;
  
  Vector<double> RetVal(4 * FilterOrder);
  //if( RetVal == NULL ) return( NULL );
  
  RetVal[2] = c[0];
  RetVal[3] = c[1];
  RetVal[0] = b[0];
  RetVal[1] = b[1];
  
  for( i = 1; i < FilterOrder; ++i )
  {
    RetVal[2*(2*i+1)]   += c[2*i] * RetVal[2*(2*i-1)]  - c[2*i+1] * RetVal[2*(2*i-1)+1];
    RetVal[2*(2*i+1)+1] += c[2*i] * RetVal[2*(2*i-1)+1] + c[2*i+1] * RetVal[2*(2*i-1)];
    
    for( j = 2*i; j > 1; --j )
    {
      RetVal[2*j]   += b[2*i] * RetVal[2*(j-1)]   - b[2*i+1] * RetVal[2*(j-1)+1] +
        c[2*i] * RetVal[2*(j-2)]   - c[2*i+1] * RetVal[2*(j-2)+1];
      RetVal[2*j+1] += b[2*i] * RetVal[2*(j-1)+1] + b[2*i+1] * RetVal[2*(j-1)] +
        c[2*i] * RetVal[2*(j-2)+1] + c[2*i+1] * RetVal[2*(j-2)];
    }
    
    RetVal[2] += b[2*i] * RetVal[0] - b[2*i+1] * RetVal[1] + c[2*i];
    RetVal[3] += b[2*i] * RetVal[1] + b[2*i+1] * RetVal[0] + c[2*i+1];
    RetVal[0] += b[2*i];
    RetVal[1] += b[2*i+1];
  }
  
  return RetVal;
}

  
} // namespace mcl
  
#endif
