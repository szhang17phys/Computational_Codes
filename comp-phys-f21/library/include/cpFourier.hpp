#ifndef CPFOURIER_H
#define CPFOURIER_H

//************************************************************
// cpFourier: Fourier Transform functions
//
// o cpFFT: simple implementation of a Fast Fourier Transform
// o cpDFT: non-optimized discrete Fourier Transform
//
// Modifications
// 06-Jan-11  add Ind2T, Ind2F, cpDFT
// 21-Dec-10  created
//************************************************************

#include <vector>

//------------------------------------------------------------
// Array ordering functions
//   Array indexing conventions
//   index   t =          f =
//   ---------------------------------
//   0       0            0
//   1       dt           +1/Ndt
//   ...     ...          ...
//   N/2-1   (N/2-1)dt    +(N/2-1)/Ndt
//   N/2     (N/2)dt      +-1/2dt
//   N/2+1   (N/2+1)dt    -(N/2-1)/Ndt
//   ...     ...          ...
//   N-1     (N-1)dt      -1/Ndt
//-----------------------------------------------------------
double Ind2T( int ind, int N, double dt );   // time from array index
int T2Ind( double t, int N, double dt );
double Ind2F( int ind, int N, double dt );   // frequency from array index
int F2Ind( double f, int N, double dt );


//------------------------------------------------------------
// cpFFT: 1D Fast Fourier Transform
// Inputs
// o real/imag: arrays containing real and imaginary values h_k(H_n)
//                array ordering (real and imag) as above
//                these are overwritten on output with H_n(h_k)
// o N          dimension of data & imag (must be an integer power of 2 !!!)
// o dir        = +1 ==> forward transform: H_n = Sum_k h_k e^{2pi i kn/N}
//              = -1 ==> reverse transform: h_k = 1/N Sum_n H_n e^{-2pi i kn/N}
// Outputs
// o real/imag: arrays containing transformed H_n
//
// Returns: true if FFT succeeded
//------------------------------------------------------------
bool cpFFT( double* real, double* imag, const int N, const int dir=1 );
inline bool cpFFT( std::vector<double>& real, std::vector<double>& imag, const int dir=1 ){
  return( cpFFT( &real[0],&imag[0],(int)real.size(),dir ) );
}

//------------------------------------------------------------
// cpDFT: 1D Discrete Fourier Transform
// Inputs/Outputs: as in cpFFT
//------------------------------------------------------------
bool cpDFT( std::vector<double>& real, std::vector<double>& imag, const int dir=1 );


#endif
