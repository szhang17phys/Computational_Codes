#include "cpFunctions.hpp"
#include "cpUtils.hpp"
using namespace std;

//============================================================
// Returns (pointer to) the function specified by id
// Output: name = name of function
//============================================================
cp1Dfunction get1Dfunction( const int id, string& name )
{
  cp1Dfunction func;

  switch( id ){
  case idPolyN:
    func = fPolyN;    name = sPolyN;    break;
  case idSqrt1mX:
    func = fSqrt1mX;  name = sSqrt1mX;  break;
  case idSin_2pix:
    func = fSin_2pix; name = sSin_2pix; break;
  case idS2_2pix:
    func = fS2_2pix;  name = sS2_2pix;  break;
  case idGauss:
    func = fGauss;  name = sGauss;  break;
  case idExp:
    func = fExp;  name = sExp;  break;
  case idEcG:
    func = fEcG;  name = sEcG;  break;
  case idExp2:
    func = fExp2;  name = sExp2;  break;
  case idExp2x2:
    func = fExp2x2;  name = sExp2x2;  break;
  case idEcGalt:
    func = fEcGalt;  name = sEcGalt;  break;
  case idPoisson:
    func = fPoisson;  name = sPoisson;  break;
  } // switch

  return func;
} // get1Dfunction

//=============================================================
// fPolyN: O(N) polynomial
//   f = Sum_{i=0,N} p_i x^i (bad form according to NR)
//=============================================================
double fPolyN( const double x, const vector<double>& p )
{
  double sum   = 0.0;
  double xpowi = 1.0;

  for( int i = 0; i < (int)p.size(); ++i ){
    sum   += p[i] * xpowi;
    xpowi *= x;
  } // i loop

  return sum;
} // fPolyN

//=============================================================
// fPoisson: Poisson distribution
//   f = nu^{k} e^{-nu} / k!   (for k = integer)
//=============================================================
double fPoisson( const double x, const vector<double>& p )
{
  int k = (int)x;
  if( k < 0 ) return -1;   // k must be >= 0

  // construct log of distribution
  double nu = p[0];
  double lnP = k*log(nu) - nu - gammln( k + 1.0 );
  return( exp(lnP) );
} // fPoisson

//============================================================
// fLegendre: Legendre polynomials
// Inputs
//   n  = order of polynomial desired
//   x  = point at which to evaluate polynomial and derviative
// Outputs
//   dPdx = derivative dP_{n}(x)/dx
// Returns: P_n(x)
//=============================================================
double fLegendre( const int n, const double x, double& dPdx )
{
  double Pm1,Pm2,P,dP;

  // n = -1 case: (only used for recursion)
  if( n == -1 ){
    dPdx = 0.0;
    return 0.0;
  } // n = -1

  // n = 0 case: P_{0} = 1
  if( n == 0 ){
    dPdx = 0.0;
    return 1.0;
  } // n = 0

  // n > 0 - use recurrence
  // P_n = [ (2n-1) x P_{n-1} - (n-1) P_{n-2} ] / n
  // dP_n/dx = n [ x P_{n} - P_{n-1} ] / (x^2 -1)
  else{
    Pm1  = fLegendre( n-1,x,dP );      // P_{n-1}
    Pm2  = fLegendre( n-2,x,dP );      // P_{n-2}
    P    = ( (2.0*n - 1.0)*x*Pm1 - (n - 1.0)*Pm2 ) / n;  // P_{n}
    dPdx = n * (x*P - Pm1) / (x*x - 1.0); // the derivative from P_{n}, P_{n-1}
    return P;
  } // n > 0
} // fLegendre

//============================================================
// fHermite: orthonormal Hermite polynomials
// Inputs
//   n  = order of polynomial desired
//   x  = point at which to evaluate polynomial and derviative
// Outputs
//   dHdx = derivative dH_{n}(x)/dx
// Returns: H^tilde_n(x)
//=============================================================
double fHermite( const int n, const double x, double& dHdx )
{
  double Hm1,dH;

  // n = -1 case (only used for recursion)
  if( n == -1 ){
    dHdx = 0.0;
    return 0.0;
  } // n = -1

  // n = 0 case: H_{0} = 1/pi^(1/4)
  if( n == 0 ){
    dHdx = 0.0;
    return 0.7511255444649425;
  } // n = 0

  // n > 0 - use recurrence
  // H_{j+1} = x sqrt( 2/(j+1) ) H_{j} - sqrt( j/(j+1) ) H_{j-1}
  else{
    Hm1  = fHermite( n-1,x,dH );       // H_{n-1}
    dHdx = sqrt( 2.0*n ) * Hm1;        // the derivative from H_{n-1}
    return( x * sqrt( 2.0/n ) * Hm1    // H_{n}
	    - sqrt( (double)(n-1)/n ) * fHermite( n-2,x,dH ) );
  } // n > 0
} // fHermite

//============================================================
// fLaguerre: Laguerre polynomials
// Inputs
//   n  = order of polynomial desired
//   x  = point at which to evaluate polynomial and derviative
// Outputs
//   dPdx = derivative dP_{n}(x)/dx
// Returns: P_n(x)
//=============================================================
double fLaguerre( const int n, const double x, double& dPdx )
{
  double Pm1,Pm2,P,dP;

  // n = -1 case: (only used for recursion)
  if( n == -1 ){
    dPdx = 0.0;
    return 0.0;
  } // n = -1

  // n = 0 case: P_{0} = 1
  if( n == 0 ){
    dPdx = 0.0;
    return 1.0;
  } // n = 0

  // n > 0 - use recurrence
  // L_n = [ (2n - 1- x) L_{n-1} - (n-1) L_{n-2} ] / n
  // dL_n/dx = n [ L_{n} - P_{n-1} ] / x
  else{
    Pm1  = fLaguerre( n-1,x,dP );      // P_{n-1}
    Pm2  = fLaguerre( n-2,x,dP );      // P_{n-2}
    P    = ( (2.0*n - 1.0 - x)*Pm1 - (n - 1.0)*Pm2 ) / n;  // L_{n}
    dPdx = n*(P - Pm1) / x; // the derivative from L_{n}, L_{n-1}
    return P;
  } // n > 0
} // fLaguerre


//=============================================================
// Gamma Functions: from NR sect 6.2
//=============================================================

//------------------------------------------------------------
// cpGammaP: user interface for Incomplete Gamma Function from NR gammp
//   P(a,x) = Int_0^x e^(-t) t^(a-1) dt / Gamma(a)
// Inputs
//   a,x  as in function above
// Returns: P(a,x) using best method
//
// Note: Guassian Quadrature method (gammpapprox) for large a
//       not yet implemented - have to check that NR method is correct
//------------------------------------------------------------
double cpGammaP( const double a, const double x )
{
  if( x < 0.0 || a <= 0.0 ) return -1;        // bad inputs
  if( x == 0.0 )            return 0.0;       // P(a,0) = 0
  else if( x < (a+1.0) )    return gser(a,x); // series rep for small x
  else              return( 1.0 - gcf(a,x) ); // continued fraction
} // cpGammaP

//------------------------------------------------------------
// cpGammaQ: user interface for Complementary Incomplete Gamma Function 
//           from NR gammq
//   Q(a,x) = 1 - P(a,x) = Int_x^infinity e^(-t) t^(a-1) dt / Gamma(a)
// Inputs
//   a,x  as in function above
// Returns: Q(a,x) using best method
//
// Note: Guassian Quadrature method (gammpapprox) for large a
//       not yet implemented - have to check that NR method is correct
//------------------------------------------------------------
double cpGammaQ( const double a, const double x )
{
  if( x < 0.0 || a <= 0.0 ) return -1;              // bad inputs
  if( x == 0.0 )            return 1.0;             // P(a,0) = 1
  else if( x < (a+1.0) ) return( 1.0 - gser(a,x) ); // series rep for small x
  else                      return gcf(a,x);        // continued fraction
} // cpGammaQ

//------------------------------------------------------------
// cpInvGammaP: user interface for Inverse Incomplete Gamma Function 
//   from NR invgammp
// uses Newton's root finding method (Halley's implementation)
//   to find x such that P(a,x) = p
// two cases are considered
// a > 1:  initial guess from Abramowitz and Stegun 26.2.22 and 26.4.17
// a <= 1: initial guess from P(a,1) ~ 0.253a + 0.12a^2
//   then P(a,x) ~ P(a,1) x^a                             x <  1
//               ~ P(a,1) + [1 - P(a,1)][1 - exp(1-x)]    x >= 1
// Inputs
//   p,a  as in function above
// Returns: x
//------------------------------------------------------------
double cpInvGammaP( const double p, const double a )
{
  const double EPS = 1.0e-8;       // accuracy is EPS^2
  double gln = gammln( a );
  double a1  = a - 1.0;

  double x,err;
  double t,u,pp,lna1,afac;

  if( p <= 0.0 ) return 0;
  if( p > 1.0 || a <= 0 ) return -1; // must have a>0, p<=1

  //-------------------
  // a > 1 initial guess from Abramowitz and Stegun
  //-------------------
  if( a > 1.0 ){
    lna1 = log( a1 );
    afac = exp( a1*(lna1 - 1.0) - gln );
    pp   = (p < 0.5) ? p : 1.0 - p;
    t    = sqrt( -2.0 * log(pp) );
    x    = ( 2.30753 + t*0.27061 ) / ( 1.0 + t*(0.99229 + t*0.4481) ) - t;
    if( p < 0.5 ) x = -x;
    x    = max( 1.0e-3, a*pow( 1.0 - 1.0/(9.0*a) - x/(3.0*sqrt(a)), 3 ) );
  } // a > 1

  //--------------------
  // a < 1: initial guess from above
  //--------------------
  else{
    t = 1.0 - a*( 0.253 + a*0.12 );
    if( p < t ) x = pow( p/t, 1.0/a );
    else        x = 1.0 - log( 1.0 - (p-t)/(1.0-t) );
  } // a < 1

  //--------------------
  // Now Halley's method root finding
  //--------------------
  for( int j = 0; j < 12; ++j ){
    if( x <= 0.0 ) return 0.0;
    err = cpGammaP(a,x) - p;
    if( a > 1.0 ) t = afac*exp( -(x-a1) + a1*(log(x) - lna1) );
    else          t = exp( -x + a1*log(x) - gln );
    u  = err/t;
    x -= ( t = u/(1.0 - 0.5*min(1.0,u*((a-1.0)/x - 1.0))) );
    if( x <= 0.0 ) x = 0.5*(x + t);
    if( abs(t) < EPS*x ) break;
  } // j loop

  return x;
} // cpInvGammaP

//-------------------------------------------------------------
// gammln: ln( Gamma(x) ) from N.R. gammln
// Gamma(x) = (x+gamma+1/2)^(x+1/2) exp( -(x+gamma+1/2) )
//          * sqrt(2pi) [ c0 + c1/(x+1) + c2/(x+2) + ... + cN/(x+N) ]
//-------------------------------------------------------------
double gammln( const double xx )
{
  double x,y,tmp,ser;
  static const double cof[14] = {
    57.1562356658629235,      -59.5979603554754912,     14.1360979747417471,
    -0.491913816097620199,    0.339946499848118887e-4,  0.465236289270485756e-4,
    -0.983744753048795646e-4, 0.158088703224912494e-3,  -0.210264441724104883e-3,
    0.217439618115212643e-3,  -0.164318106536763890e-3, 0.844182239838527433e-4,
    -0.261908384015814087e-4, 0.368991826595316234e-5
  };

  if( xx <= 0 ) return 0;

  y = x = xx;
  tmp = x + 671.0/128.0;
  tmp = ( x+0.5 )*log( tmp ) - tmp;

  ser = 0.999999999999997092;
  for( int i = 0; i < 14; ++i ) ser += cof[i]/++y;
  return( tmp + log( cpSqrt2PI * ser / x ) );
} // gammln

//-------------------------------------------------------------
// gser: Series representation of P(a,x) - internal use only
// gamma(a,x) = e^(-x) x^a Sum_n=0^infinity x^n Gamma(a) / Gamma(a+1+n)
//   from NR eqn 6.2.5
//-------------------------------------------------------------
double gser( const double a, const double x )
{
  double ap = a;
  double del = 1.0 / a;
  double sum = del;
  double gln = gammln(a);
  for( ;; ){
    ++ap;              
    del *= x/ap;
    sum += del;
    if( abs(del) < abs(sum)*Gamma_EPS )
      return( sum * exp( -x + a*log(x) - gln ) );
  } // for loop
} // gser

//-------------------------------------------------------------
// gcf: continued fraction rep of Q(a,x) - internal use only
// Gamma(a,x) ~ e^(-x) x^a [ 1/(x+1-a-) 1.(1-a)/(x+3-a-) 2.(2-a)/(x+5-a-)... ]
//   from NR eqn 6.2.7
// This evaluates the continued fraction by a modified Lenz's method
//   from NR sect 5.2
//-------------------------------------------------------------
double gcf( const double a, const double x )
{
  double gln = gammln(a);
  double b   = x + 1.0 - a;
  double c   = 1.0 / Gamma_FPMIN;
  double d   = 1.0 / b;
  double h   = d;
  double an,del;
  for( int i = 1; ; ++i ){
    an = -i * ( i - a );
    b += 2.0;
    d  = an*d + b;
    if( abs(d) < Gamma_FPMIN ) d = Gamma_FPMIN;
    c  = b + an/c;
    if( abs(c) < Gamma_FPMIN ) c = Gamma_FPMIN;
    d   = 1.0 / d;
    del = d * c;
    h  *= del;
    if( abs(del - 1.0) <= Gamma_EPS ) break;
  } // i loop
  return( exp( -x + a*log(x) - gln ) * h );
} // gcf

//===========================================================
// sortSI: straight insertion sort (efficient for N<20)
//   sorts array into ascending numerical order
//===========================================================
template<class T>
void sortSI( vector<T> &arr )
{
  int i,n = arr.size();
  T a;
  for( int j = 1; j < n; j++ ){
    a = arr[j];
    i = j;
    while( i > 0 && arr[i-1] > a ){
      arr[i] = arr[i-1];
      i--;
    } // while i
    arr[i] = a;
  } // j
} // sortSI
