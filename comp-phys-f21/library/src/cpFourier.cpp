#include <cmath>

#include "cpFourier.hpp"
#include "cpUtils.hpp"
using namespace std;

//===================================================================
// Array Indexing
//   assumes t_0 = 0
//===================================================================
double Ind2T( int ind, int N, double dt ){ return( ind * dt ); }

int T2Ind( double t, int N, double dt )
{
  if( t < 0 ) return -1;              // t < 0
  for( int j = 0; j < N; ++j )        // t in region
    if( t < (j+1)*dt ) return j;      //  return j if t < t_{j+1}
  return -1;                          // t > upper bound
} // T2Ind

double Ind2F( int ind, int N, double dt )
{
  double df = 1. / (N*dt);
  if( ind <= N/2 ) return( ind * df );     // f >= 0 & f = +=f_c
  else             return( (ind-N) * df ); // f <  0
} // Ind2F

int F2Ind( double f, int N, double dt )
{
  double df = 1. / (N*dt);
  if( f >= 0 ){                            // positive frequencies (and 0)
    for( int j = 0; j <= N/2; ++j )
      if( f < (j+1)*df ) return j;         //  return j if f < f_{j+1}
    return -1;                             // f outside of bound
  } // f >= 0

  else{                                    // negative frequencies
    for( int j = N/2; j < N; --j )
      if( f < (j-N+1)*df ) return j;
    return -1;
  } // f < 0
} // F2Ind

//===================================================================
// cpFFT: 1D Fast Fourier Transform
//   based on four1 in Numerical Recipes 12.2.1
//   and SimpleFFT in CSM app. 9B 
//===================================================================
bool cpFFT( double* real, double* imag, const int N, const int dir )
{
  // forward or reverse tranform
  int sgn = dir / abs(dir);

  //----------------------
  // check that N = 2^m & determine m
  //----------------------
  int pow = 1;                              // pow will hold m
  int N2  = N/2;
  int tmp = N2;
  if( N < 2 || N&(N-1) ) return( false );   // N & (N-1) only zero if N=2^m
  while( (tmp >>= 1) > 0 ) pow++;           // divide by 2

  //----------------------
  // Bit-Reverse input arrays
  //----------------------
  int mm;
  int jj = N2;
  for( int i = 1; i < N-1; ++i ){    // [0] and [N-1] comp's never reversed
    if( jj > i ){
      cpSwap( real[jj],real[i] );
      cpSwap( imag[jj],imag[i] );
    } // i < jj

    // update jj to be the bit-reversal of next i
    // (not sure why this algorithm works generally)
    mm = N2;
    while( (mm >= 1) && (mm <= jj) ){
      jj -= mm;        // j = j-m
      mm >>= 1;        // m = m/2
    } // while
    jj += mm;
  } // i loop

  //-----------------------------
  // Danielson-Lanczos Section
  //-----------------------------
  int inc,ip;
  double theta,cosTh,sinTh,wr,wi,wtmp,tmpRe,tmpIm;

  // Loop over stages of the FFT
  // in each stage there will be N/2 2-point Transforms
  // jj picks the array elements used for the 2-point transforms (butterfly)
  // - corresponds to the difference in array index between the elements used 
  //   for the current butterfly, i.e. 1,2,4,...
  jj = 1;
  for( int p = 1; p <= pow; ++p ){
    inc   = jj << 1;             // increment to next butterfly
    wr    = 1.0;                 // Re( W^0 )
    wi    = 0.0;                 // Im( W^0 )
    theta = sgn * cpPI / jj;     // base argument for calc of W_N^j
    cosTh = cos( theta );        //   each subsequent power of W_N will be
    sinTh = sin( theta );        //   multiplied by e^{i theta}

    // j = groups of butterflies w/ different powers of W_N multiplying
    // 2nd element
    for( int j = 0; j < jj; ++j ){
      for( int i = j; i < N; i += inc ){ // i  = index of 1st elem in butterfly
	ip    = i + jj;                  // ip = index of 2nd elem in butterfly

	// multiply 2nd butterfly element by appropriate power of W_N
	tmpRe = wr*real[ip] - wi*imag[ip]; // R(X_1)*R(W_N^j) - I(X_1)*I(W_N^j)
	tmpIm = wi*real[ip] + wr*imag[ip]; // R(X_1)*I(W_N^j) + I(X_1)*R(W_N^j)

	// the butterfly
	real[ip] = real[i] - tmpRe;        // x_1 = x_0 - x_1
	imag[ip] = imag[i] - tmpIm;
	real[i] += tmpRe;                  // x_0 = x_0 + x_1
	imag[i] += tmpIm;
      } // i loop

      // update power of W_N for 2nd element of next butterfly group
      //   W = W * e^{i theta}
      wtmp = wr;
      wr   = wr*cosTh - wi*sinTh;    // Re(W)
      wi   = wi*cosTh + wtmp*sinTh;  // Im(W)
    } // j loop

    jj = inc;  // update difference in array elements for next stage
  } // p loop

  //-------------------------
  // multiply by 1/N for reverse FFT
  //-------------------------
  if( dir < 0 ){
    for( int n = 0; n < N; ++n ){
      real[n] = real[n] / N;
      imag[n] = imag[n] / N;
    } // n loop
  } // dir < 0

  return( true );
} // cpFFT

//==================================================================
// cpDFT: 1D Discrete Fourier Transform (not fast)
//==================================================================
bool cpDFT( vector<double>& real, vector<double>& imag, const int dir )
{
  int N   = (int)real.size();
  int sgn = dir / abs( dir );
  vector<double> rTmp( real );
  vector<double> iTmp( imag );

  double theta,sinTh,cosTh;

  // loop over elements of output vectors
  for( int j = 0; j < N; ++j ){
    real[j] = 0;
    imag[j] = 0;

    // FT for element j
    for( int k = 0; k < N; ++k ){
      theta    = sgn * 2. * cpPI * j * k / N;
      cosTh    = cos( theta );
      sinTh    = sin( theta );
      real[j] += cosTh*rTmp[k] - sinTh*iTmp[k];  // e^{2pi i jk/N} h_k
      imag[j] += cosTh*iTmp[k] + sinTh*rTmp[k];
    } // k loop

  } // j loop

  //-------------------------
  // multiply by 1/N for reverse FFT
  //-------------------------
  if( dir < 0 ){
    for( int n = 0; n < N; ++n ){
      real[n] = real[n] / N;
      imag[n] = imag[n] / N;
    } // n loop
  } // dir < 0

  return( true );
} // cpDFT
