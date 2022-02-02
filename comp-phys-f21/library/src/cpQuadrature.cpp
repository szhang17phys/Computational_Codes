#include <cmath>

#include "cpQuadrature.hpp"
using namespace std;

//======================================================================
// cpTrapezoid: based on Numerical Recipes Trapzd
//======================================================================
double cpTrapezoid::next(){
  double x,sum,hx2;

  _nIter++;  // update iteration number to value for this step

  //-----------------
  // n=1: evaluate integral using two endpoints (a,b)
  //      integral = (1/2)(f(b)+f(a)) * (b-a)
  //-----------------
  if( _nIter == 1 ){
    _to_add = 1;   // number of points to add in next iteration
    return( _s = 0.5*(_b - _a)*(func(_b) + func(_a)) ); // the base integral
  } // _nIter = 1

  //----------------
  // n > 1: add 2^(n-2) new points to existing integral
  //        I_n = (1/2)I_{n-1} + (b-a)/2^(n-1) Sum_{i=new-pts} f_{i}
  //----------------
  else{
    // spacing of pts to be added (2^{n-2} new points)
    // this is twice the step size (h) for this iteration
    hx2 = ( _b - _a ) / (double)_to_add;

    // construct Sum_{i=new-pts} f_{i}
    x   = _a + 0.5*hx2; // 1st new point
    sum = 0.0;
    for( int j = 0; j < _to_add; ++j, x+=hx2 ) sum += func(x);

    // update in preparation for next iteration
    _to_add *= 2;      // number of points to add in next iteration
    return( _s = 0.5*(_s + hx2*sum) );  // new value for integral
  } // _nIter > 1

} // next

//======================================================================
// qRect: integrates function f from a to b using rectangular approx
//   note - only implemented for a fixed number of steps
//======================================================================
double qRect( int& j, cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const int n )
{
  if( n < 0 ) return -1;

  j = n;
  int nInterv = (int)pow( 2.0,n );
  double dx = ( b - a ) / (double)nInterv;

  // accumulate f(x) at low edges of intervals
  double x=a,s=0.0;
  for( int i = 0; i < nInterv; ++i, x+=dx ) s += f(x,p);

  return( s*dx );
} // qRect

double qRect( cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const int n ){
  int j; return( qRect(j,f,p,a,b,n) ); }

//======================================================================
// qTrap: integrates function f from a to b using trapezoid approx.
//======================================================================
double qTrap( int& j, cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps, const int n )
{
  cpTrapezoid t( f,p,a,b );
  int nStep = (n < 0) ? JMAX : n;  // number of steps to take

  // Add points to integral
  double s=-1.0,olds=0.0;  // current and prev value of integral
  for( j = 0; j <= nStep; ++j ){
    s = t.next();              // add next set of points to integral
    if( n < 0 && j > JMIN ){   // test for convergence if requested (n<0)
      if( converge( s,olds,eps ) ) return s;
      olds = s;
    } // n<0 && j>JMIN
  } // j = steps

  if( n < 0 )                           // integral did not converge
    throw( "Too many steps in qTrap" );
  return s;                             // fixed number of steps
} // qTrap

double qTrap( cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps ){
  int j; return( qTrap(j,f,p,a,b,eps) ); }

//======================================================================
// qSimp: integrates function f from a to b using Simpson's rule
//======================================================================
double qSimp( int& j, cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps, const int n )
{
  if( n == 0 ) return 0;           // need at least 2 intervals for Simpson

  cpTrapezoid t( f,p,a,b );
  int nStep = (n < 0) ? JMAX : n;  // number of steps to take

  // Add points to integral
  double s,os=0.0;       // current and prev value of integral
  double st,ost=0.0;     // current and prev value of trapezoid approx
  for( j = 0; j <= nStep; ++j ){
    st = t.next();               // add next set of points to trapezoid integral
    s  = ( 4.0*st - ost ) / 3.0; // Simpson's rule

    if( n < 0 && j > JMIN )      // test for convergence if requested (n<0)
      if( converge( s,os,eps ) ) return s;

    os  = s;    // save this step to be used as "old" in next step
    ost = st;
  } // j = steps

  if( n < 0 )                           // integral did not converge
    throw( "Too many steps in qSimp" );
  return s;                             // fixed number of steps (n>0)
} // qSimp

double qSimp( cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps ){
  int j; return( qSimp(j,f,p,a,b,eps) ); }

//======================================================================
// qRomb: integrates function f from a to b using Romberg method
//======================================================================
double qRomb( int& j, cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps )
{
  cpTrapezoid t( f,p,a,b );

  // initialize
  double ss;
  vector<double> s( JMAXP,0 );     // succesive trapezoid approximations
  vector<double> h( JMAXP,0 );     // relative step sizes (compared to b-a)
  h[0] = 1.0;

  //-------------------
  // Add points to integral
  //-------------------
  double dy;
  for( j = 1; j <= JMAX; ++j ){
    s[j-1] = t.next();       // integral estimate for current step (j-1)

    // extrapolate to zero step size using previous K points
    if( j >= K ){
      ss = PolyInterp( h,s,j-K,K,0.0,dy );

      // convergence based on error in extrapolation
      // dy O(K) - O(K-1) extrapolation
      if( abs(dy) <= eps*abs(ss) || abs(dy) <= eps ) return ss;
    } // j >= JMIN

    // next point (relative step size) for extrapolation
    // this is 1/4 of previous step size to make extrapolation a polynomial
    //   in h^2
    h[j] = 0.25 * h[j-1];
  } // j = steps

  throw( "Too many steps in qRomb" );
} // qRomb

double qRomb( cp1Dfunction f, const vector<double>& p,
	      const double a, const double b, const double eps ){
  int j; return( qRomb(j,f,p,a,b,eps) ); }


//=======================================================================
// Gaussian Quadrature Functions
//=======================================================================

//------------------------------------------------------
// qGauss: does the Gaussian Quadrature integration using method specified
//------------------------------------------------------
double qGauss( int& j, cp1Dfunction f, const vector<double>& p,
	       const double a, const double b,
	       const int idMethod, const double eps, const int n )
{
  if( n == 0 ) return 0;           // need at least 1 point

  vector<double> x;
  vector<double> w;
  double s = 0,os = -cpHUGE;

  //-----------------
  // Select Method
  //-----------------
  cpGaussWtRt fWtRt;
  cpWeightFunc fWeight;
  switch( idMethod ){
  case idGaHe:   // Gauss-Hermite
    fWtRt = GaussHermite;  fWeight = fWeightGaHe;  break;
  case idGaLa:   // Gauss-Laguerre
    fWtRt = GaussLaguerre;  fWeight = fWeightGaLa;  break;
  default:       // Gauss-Legendre
    fWtRt = GaussLegendre;  fWeight = fWeightGaLe;  break;
  } // switch

  //----------------------
  // Gaussian Quadrature for number of points specified
  // or until convergence
  //----------------------
  int jbeg = ( n < 0 ) ? 1    : n;
  int jend = ( n < 0 ) ? JMAXGAUSS : n+1;
  for( j = jbeg; j < jend; ++j ){
    x.resize( j,0.0 );       // set up root and weight vectors
    w.resize( j,0.0 );       //   for this number of points
    fWtRt( a,b,x,w );        // extract roots & weights using specified method

    s = 0.0;
    for( int i = 0; i < j; ++i )            // construct Sum w_i f(x_i)
      s += w[i] * f( x[i],p ) / fWeight( x[i] );

    if( n < 0 )
      if( converge( s,os,eps ) ) return s;  // test for conv (if n<0)

    os = s;    // save this step to be used as "old" in next step
  } // j loop

  if( n < 0 )                           // integral did not converge
    throw( "Too many steps in qGauss" );

  j = n;        // fixed number of steps
  return s;
} // qGauss


//-------------------------------------------------------
// GaussLegendre: adapted from NR gauher
//-------------------------------------------------------
void GaussLegendre( const double x1, const double x2, 
		    vector<double>& x, vector<double>& w )
{
  // order of Legendre polynomial to be used from size of vectors
  int n = x.size();

  //---------------------
  // Loop over roots on the standard [-1,1] interval
  // note: roots are symmetric about midpoint
  // - so we only have to do half of them
  //---------------------
  double z,z1;            // estimate and previous estimate of root
  double Pj,dPj;          // Legendre poly of order j and derivative
  int m = ( n + 1 ) / 2;  // n=even --> no root at 0; n=odd --> root at zero
  int its = 0;

  double xm = ( x2 + x1 ) / 2.0;  // midpoint
  double xl = ( x2 - x1 ) / 2.0;  // difference

  for( int irt = 0; irt < m; ++irt ){
    //----------------
    // construct initial guess for roots
    //   couldn't find from where these are obtained
    //----------------
    z = cos( cpPI * (irt+0.75) / (n+0.5) );

    //---------------
    // refine the root estimate using Newton's method
    // use absolute convergence criterion since we will have roots at 0
    //---------------
    for( its = 0; its < MAXITGAUSS; ++its ){
      Pj = fLegendre( n,z,dPj );          // value of H_{n} at current root est
      z1 = z;                            // store current estimate of root
      z  = z1 - Pj/dPj;                  // new root est: Newton-Raphson formula
      if( abs( z - z1 ) <= EPSGAUSS ) break;  // convergence achieved
    } // its (refinement) loop
    if( its >= MAXITGAUSS ) throw( "GaussLegendre: too many iterations" );

    //--------------
    // Scale the roots to the desired interval
    //   x_remap = x (b-a)/2 + (b+a)/2
    //   w_remap = w (b-a)/2
    // - note that roots are symmetric about the midpoint
    //--------------
    x[irt]     = xm - z*xl;                         // scaled root
    w[irt]     = xl * 2.0/((1.0 - z*z) * dPj*dPj);  // scaled weight
    if( n-1-irt != irt ){        // symmetric roots
      x[n-1-irt] = xm + z*xl;    //   don't need to redo the middle one
      w[n-1-irt] = w[irt];       //   for n = odd
    } // n-1-irt != irt
  } // irt (root) loop

} // GaussLegendre


//-------------------------------------------------------
// GaussHermite: adapted from NR gauher
//-------------------------------------------------------
void GaussHermite( const double x1, const double x2, 
		   vector<double>& x, vector<double>& w )
{
  // order of Hermite polynomial to be used from size of vectors
  int n = x.size();

  //---------------------
  // Loop over roots - starting from largest root
  // note: roots are symmetric about 0 - so we only have to do half of them
  //---------------------
  double z,z1;            // estimate and previous estimate of root
  double Pj,dPj;          // Hermite poly of order j and derivative
  int m = ( n + 1 ) / 2;  // n=even --> no root at 0; n=odd --> root at zero
  int its = 0;

  for( int irt = 0; irt < m; ++irt ){
    //----------------
    // construct initial guess for roots
    //   couldn't find from where these are obtained
    //----------------
    if( irt == 0 )                         // largest root
      z = sqrt( (double)(2*n+1) ) - 1.85575*pow( (double)(2*n+1), -0.16667 );
    else if( irt == 1 )                    // 2nd largest root
      z -= 1.14 * pow( (double)n,0.426 ) / z;
    else if( irt == 2 )                    // 3rd largest root
      z = 1.86*z - 0.86*x[0];
    else if( irt == 3 )                    // 4th largest root
      z = 1.91*z - 0.91*x[1];
    else                                 // all others
      z = 2.0*z - x[irt-2];

    //---------------
    // refine the root estimate using Newton's method
    // use absolute convergence criterion since we will have roots at 0
    //---------------
    for( its = 0; its < MAXITGAUSS; ++its ){
      Pj = fHermite( n,z,dPj );          // value of H_{n} at current root est
      z1 = z;                            // store current estimate of root
      z  = z1 - Pj/dPj;                  // new root est: Newton-Raphson formula
      if( abs( z - z1 ) <= EPSGAUSS ) break;  // convergence achieved
    } // its (refinement) loop
    if( its >= MAXITGAUSS ) throw( "GaussHermite: too many iterations" );

    //--------------
    // Store the root, its symmetric (negative) counterpart, and the weights
    //--------------
    x[irt]     = z;                     // the roots
    x[n-1-irt] = -z;
    w[irt]     = 2.0 / ( dPj*dPj );     // the weights
    w[n-1-irt] = w[irt];
  } // irt (root) loop

} // GaussHermite


//-------------------------------------------------------
// GaussLaguerre: adapted from NR gauher
//-------------------------------------------------------
void GaussLaguerre( const double x1, const double x2, 
		   vector<double>& x, vector<double>& w )
{
  // order of Laguerre polynomial to be used from size of vectors
  int n = x.size();

  //---------------------
  // Loop over roots - starting from smallest root
  // note: roots do not exhibit a nice symmetry 
  // - so we have to do them all
  //---------------------
  double z,z1;            // estimate and previous estimate of root
  double Pj,dPj;          // Laguerre poly of order j and derivative
  int its = 0;

  for( int irt = 0; irt < n; ++irt ){
    //----------------
    // construct initial guess for roots
    //   couldn't find from where these are obtained
    //----------------
    if( irt == 0 )                         // smallest root
      z = 3.0 / ( 1.0 + 2.4*n );
    else if( irt == 1 )                    // 2nd smallest root
      z += 15.0 / ( 1.0 + 2.5*n );
    else                                   // all others
      z += (z - x[irt-2]) * (1.0 + 2.55*(irt-1)) / (1.9*(irt-1));

    //---------------
    // refine the root estimate using Newton's method
    // use absolute convergence criterion since we will have roots at 0
    //---------------
    for( its = 0; its < MAXITGAUSS; ++its ){
      Pj = fLaguerre( n,z,dPj );          // value of H_{n} at current root est
      z1 = z;                            // store current estimate of root
      z  = z1 - Pj/dPj;                  // new root est: Newton-Raphson formula
      if( abs( z - z1 ) <= EPSGAUSS ) break;  // convergence achieved
    } // its (refinement) loop
    if( its >= MAXITGAUSS ) throw( "GaussLaguerre: too many iterations" );

    //--------------
    // Store the root & weight
    //--------------
    x[irt] = z;                     // the root
    w[irt] = 1.0 / ( z*dPj*dPj );   // the weight
  } // irt (root) loop

} // GaussLaguerre
