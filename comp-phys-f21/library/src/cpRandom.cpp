//**************************************************************
// cpRandom: non-inline formulas
//**************************************************************
#include <cmath>

#include "cpRandom.hpp"
#include "cpFunctions.hpp"
using namespace std;

//============================================================
// randDevNrm: normal distribution w/ mean=0, sigma=1
//   p(x) dx = exp( -x^2 /2 ) / sqrt( 2pi )
//   Interval: -infinity - +infinity
//   Uses: Leva's Ratio-of-Uniforms method from Numerical Recipes
//============================================================
double cpRandom::randDevNrm()
{
  const double a = 0.19600;
  const double b = 0.25472;
  const double s = 0.449871;
  const double t = 0.386595;
  const double bnd = 0.85776;   // sqrt(2/e)
  const double Qlo = 0.27597;
  const double Qhi = 0.27846;

  double u,v,x,y,Q;

  // loop until we get a point in region
  do{
    // generate uniform deviates u,v
    u = randDouble( 0.0,1.0 );    // u: [0,1]
    v = randDouble( -bnd,bnd );   // v: [-sqrt(2/e),sqrt(2/e)]

    // variables used in squeezes
    x = u - s;
    y = abs(v) + t;
    Q = x*x + y*( a*y - b*x );
  } while( Q > Qlo &&
	   ( Q > Qhi || v*v > -4.0 * log(u) * u*u ) );
  // the above checks whether
  //   1. Q is inside the inner squeeze (Qlo)  --> done 
  //   2. Q is outside the outer squeeze (Qhi) --> generate another point
  //   3. Q is between the squeezes            --> check the real boundary

  return( v/u );  // return x = v/u
} // randDevNrm

//============================================================
// randDevNrmBM: normal distribution w/ mean=0, sigma=1
//   p(x) dx = exp( -x^2 /2 ) / sqrt( 2pi )
//   Interval: -infinity - +infinity
//   Uses: Box-Muller method from Numerical Recipes
//============================================================
double cpRandom::randDevNrmBM()
{
  double v1,v2,rsq,y1,fac;

  //-----------------------
  // no extra deviate available: generate new pair
  //-----------------------
  if( _y2 == 0. ){
    do{
      v1  = randDouble( -1.0,1.0 );      // pick two uniform numbers [-1,1]
      v2  = randDouble( -1.0,1.0 );
      rsq = v1*v1 + v2*v2;               // in unit circle?
    } while( rsq >= 1.0 || rsq == 0. );  // bad radius generated

    // new values:
    // y1 = sqrt( -2 lnx1 ) cos( 2pi x2)
    // y2 = sqrt( -2 lnx1 ) sin( 2pi x2)
    fac = sqrt( -2.0 * log( rsq ) / rsq );
    y1  = v2 * fac;
    _y2 = v1 * fac;  // store y2 for next time
  } // _y2 == 0

  //-------------------
  // extra deviate available: use previously stored y2
  //-------------------
  else{
    y1  = _y2;
    _y2 = 0.;   // will need new values next time
  } // _y2 != 0

  return( y1 ); // deviate centered at 0 w/ width 1
} // randDevNrmBM


//============================================================
// randDevExp: exponential decay distribution from Transformation Method
//   p(x) dx = lambda exp( -lambda x ) dx
//   x       = -lambda ln( r )
//============================================================
double cpRandom::randDevExp( double lambda )
{
  double r = randDouble( 0.0,1.0 );  // uniform deviate on [0,1]
  return( -log(r) / lambda );        // transformation
} // randDevExp


//============================================================
// randDevPoisson: Poisson Distribution
//   q(k) dk = lambda^([k]) exp( -lambda ) / [k]!  dk
//   Interval: k = 0,1,2,.....,infinity
//   Uses: various methods from Numerical Recipes
//============================================================
int cpRandom::randDevPoisson( double lambda )
{
  double u,u2,v,v2,p,t;
  double explam,sqlam,loglam;
  int k;

  //-----------------------
  // Small lambda case: use product-of-uniforms method
  //-----------------------
  if( lambda < 5. ){
    explam = exp( -lambda );
    k = -1;
    t = 1;
    do{  // sum exponential time intervals until Sum > lambda
      ++k;
      t *= randDouble();
    } while( t > explam );
  } // lambda < 5

  //-----------------------
  // Large lambda case: use ratio-of-uniforms method
  //-----------------------
  else{
    sqlam = sqrt( lambda );
    loglam = log( lambda );

    for( ;; ){
      u = 0.64 * randDouble();
      v = -0.68 + 1.28*randDouble();

      //---------------
      // fast rejection: outside of outer squeeze
      //---------------
      if( lambda > 13.5 ){
	v2 = v*v;
	if( v >= 0. ){        // above upper part of outer squeeze
	  if( v2 > 6.5*u*( 0.64-u )*( u+0.2 ) ) continue; }
	else{                 // below lower part of outer squeeze
	  if( v2 > 9.6*u*( 0.66-u )*( u+0.07 ) ) continue; }
      } // lambda > 13.5

      k = (int)floor( sqlam*(v/u) + lambda + 0.5 ); // floor for negatives
      if( k < 0 ) continue;
      u2 = u*u;

      //--------------
      // fast accept: inside inner squeeze
      //--------------
      if( lambda > 13.5 ){
	if( v >= 0. ){        // below upper part of inner squeeze
	  if( v2 < 15.2*u2*( 0.61-u )*( 0.8-u ) ) break; }
	else{                 // above lower part of inner squeeze
	  if( v2 < 6.76*u2*( 0.62-u )*( 1.4-u ) ) break; }
      } // lambda > 13.5

      //--------------
      // between squeezes: must do full calculation of p(k)
      // note: don't store ln(Gamma) as in NR ==> slower execution
      //--------------
      p = sqlam * exp( -lambda + k*loglam - gammln(k+1.) );
      if( u2 < p ) break;
    } // for
  } // lambda > 5

  return k;
} // randDevPoisson

//=====================================================================================
// cpRandMT: Mersenne Twister functions
//=====================================================================================
// Set the random generator sequence
void cpRandMT::init( Ullong seed )
{
  _fCount = 312;
  _fMt[0] = _seed;

  for ( int i = 1; i < 312; i++ ) 
    _fMt[i] =  (6364136223846793005ULL * (_fMt[i-1] ^ (_fMt[i-1] >> 62)) + i);
} // cpRandMT::init

// The 64-bit generator
Ullong cpRandMT::randInt64()
{
  Ullong y;

  const int kN = 312;
  const int kM = 156;
  const Ullong kMatrixA   = 0xb5026f5aa96619e9;
  const Ullong kUpperMask = 0xffffffff80000000;
  const Ullong kLowerMask = 0x7fffffff;

  static Ullong mag01[2] = { 0ULL, kMatrixA };

  if ( _fCount >= kN ){
    register int i;

    for ( i = 0; i < kN-kM; i++ ) {
      y = (_fMt[i] & kUpperMask) | (_fMt[i+1] & kLowerMask);
      _fMt[i] = _fMt[i+kM] ^ (y >> 1) ^ mag01[(int)(y & 1ULL)];
    }

    for (   ; i < kN-1    ; i++) {
      y = (_fMt[i] & kUpperMask) | (_fMt[i+1] & kLowerMask);
      _fMt[i] = _fMt[i+kM-kN] ^ (y >> 1) ^ mag01[(int)(y & 1ULL)];
    }

    y = (_fMt[kN-1] & kUpperMask) | (_fMt[0] & kLowerMask);
    _fMt[kN-1] = _fMt[kM-1] ^ (y >> 1) ^ mag01[(int)(y & 1ULL)];

    _fCount = 0;
  } // _fCount >= kN

  y  = _fMt[_fCount++];
  y ^= (y >> 29) & 0x5555555555555555ULL;
  y ^= (y << 17) & 0x71d67fffeda60000ULL;
  y ^= (y << 37) & 0xfff7eee000000000ULL;
  y ^= (y >> 43);

  return y;

} // cpRandMT::randInt64
