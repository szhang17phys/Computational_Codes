#ifndef CPQUADRATURE_H
#define CPQUADRATURE_H

//************************************************************************
// cpQuadrature: numerical integration functions
//
// Modifications
// 28-Aug-18  includes all in this file
// 18-Feb-10  add GaussLaguerre
// 16-Feb-10  add cpGaussWtFunc, qGauss
//            add idXXX (method IDs), JMAXGAUSS
// 09-Feb-10  add base class cpQuadrature
//            change variables _it to _to_add; _n to _nIter
// 08-Feb-10  add absolute test to converge() and to qRomb()
// 14-Oct-09  created
//************************************************************************

#include <cmath>
#include <limits>

#include "cpUtils.hpp"
#include "cpFunctions.hpp"
#include "cpInterpolate.hpp"

const int JMAX  = 20;       // max steps = 2^JMAX
const int JMAXP = JMAX + 1;
const int JMIN  = 5;        // min # of steps (2^5): avoids spurious converg
const int K     = 5;        // for use in Romberg: O(4) interpolating function
const int JMAXGAUSS  = 25;        // for Gaussian Quadrature
const int MAXITGAUSS = 10;        //   max iterations for root-finding
const double EPSGAUSS = 1.0e-14;  //   convergence condition for root-finding
const double CONV = 1.0e-10;  // convergence condition

// Integration Method IDs
const int idRect = 0;   // Rectangle
const int idTrap = 1;   // Trapezoid
const int idSimp = 2;   // Simpson
const int idRomb = 3;   // Romberg
const int idGaLe = 4;   // Gauss-Legendre
const int idGaHe = 5;   // Gauss-Hermite
const int idGaLa = 6;   // Gauss-Laguerre
const std::string sMethod[10] = { "Rectangle","Trapezoid","Simpson","Romberg",
			     "Gauss-Legendre","Gauss-Hermite","Gauss-Laguerre" };

//=====================================================================
// typedef specifying form of functions returning roots and weights
//   associated with orthogonal polynomials used in Gaussian Quadrature
// Input
//   arg 1 = lower bound of region (if used)
//   arg 2 = upper bound of region (if used)
// Output
//   arg 3 = vector of roots in region
//   arg 4 = vector of weights for Gaussian quadrature
//=====================================================================
typedef void (*cpGaussWtRt)( const double, const double, 
			     std::vector<double>&, std::vector<double>& );

// typedef for the weight functions for different Gauss-Quadrature methods
typedef double (*cpWeightFunc)( const double );


//=======================================================================
// cpQuadrature: base class for numerical integration
// Inputs
//   f = function to integrate
//   p = parameters used by function
//   a,b = lower and upper bounds of integration
//=======================================================================
class cpQuadrature {
public:
  cpQuadrature() {};
  cpQuadrature( cp1Dfunction f, const std::vector<double>& p,
	       const double a, const double b ) :
    _f(f), _p(p), _a(a), _b(b) { _nIter = 0; _to_add = 0; }

  // interface to function used in evaluating the integral
  // this gets changed in the daughter class if we have to transform the func
  virtual double func( const double x ){ return( fBase(x) ); }

  // interface to function that hides the (const) parameters
  double fBase( const double x ){ return( _f(x,_p) ); }

  // updates integral to include new points (defined in daughters)
  virtual double next() { return _s; }

  // this allows to reset limits for use with a transformed function
  void resetLimits( const double a, const double b ) { _a = a; _b = b; }

// protected (rather than private) allows daughter classes to access these vars
protected:
  cp1Dfunction _f;   // function to integrate
  std::vector<double> _p; // parameters of function
  double _a,_b;      // lower and upper bounds
  double _s;         // current estimate of integral
  // current refinement level & points to add
  int _nIter;        // for standard Trapezoid: N(points) = 2^{_n-1} + 1
  int _to_add;       // for standard Trapezoid: 2^{n-2})
}; // cpQuadrature

//=======================================================================
// cpTrapezoid: base actions for the standard Trapezoid rule
// Inputs: same as cpQuadrature
//=======================================================================
class cpTrapezoid : public cpQuadrature {
public:
  cpTrapezoid() {};
  cpTrapezoid( cp1Dfunction f, const std::vector<double>& p,
	       const double a, const double b ) : cpQuadrature(f,p,a,b) {};

  // use just the base function for normal Trapezoid rule
  double func( const double x ){ return( fBase(x) ); }

  // updates integral to include new points
  double next();

}; // cpTrapezoid


//======================================================================
// Numerical Integration algorithms
// Come in two versions: normal use & test
//   j   = (OUTPUT) number of steps taken to convergence (test only)
//   f   = function to integrate
//   p   = vector of function parameters
//   a,b = limits of integration
//   eps = convergence if new/old < eps
//   n   = number of intervals = 2^n (test only)
//         if n=-1 (default) integrate until convergence
//         note: Romberg method does not use this
//   idMethod = ID for (Gaussian Quadrature) integration method
//======================================================================
double qRect( int& j, cp1Dfunction f, const std::vector<double>& p,  // Rectangular
	      const double a, const double b, const int n );
double qRect( cp1Dfunction f, const std::vector<double>& p,
	      const double a, const double b, const int n );

double qTrap( int& j, cp1Dfunction f, const std::vector<double>& p,  // Trapezoid
	      const double a, const double b,
	      const double eps=CONV, const int n=-1 );
double qTrap( cp1Dfunction f, const std::vector<double>& p,
	      const double a, const double b, const double eps=CONV );

double qSimp( int& j, cp1Dfunction f, const std::vector<double>& p,  // Simpson
	      const double a, const double b,
	      const double eps=CONV, const int n=-1 );
double qSimp( cp1Dfunction f, const std::vector<double>& p,
	      const double a, const double b, const double eps=CONV );

double qRomb( int& j, cp1Dfunction f, const std::vector<double>& p,  // Romberg
	      const double a, const double b, const double eps=CONV );
double qRomb( cp1Dfunction f, const std::vector<double>& p,
	      const double a, const double b, const double eps=CONV );

double qGauss( int& j, cp1Dfunction f, const std::vector<double>& p,  // Gauss
	       const double a, const double b,
	       const int idMethod=idGaLe, 
	       const double eps=CONV, const int n=-1 );

//======================================================================
// converge: common test for convergence by comparing change between 
//             new(s1) and old(s0) result to a fraction of old result
//======================================================================
inline bool converge( const double s1, const double s0, const double eps ){
  return( ( (std::abs(s1-s0) < eps*std::abs(s0)) || (std::abs(s1-s0) < eps)
	    || (s1==0.0 && s0==0.0)  ) ); }

//======================================================================
// Functions related to Gaussian Quadrature
//======================================================================

//-------------------------------
// GaussLegendre: roots and weights for Gauss-Legendre quadrature
// Input: x1,x2 - mapped to  Gauss-Hermite assumes -1 - +1
// Gauss-Legendre W(x) = 1
//-------------------------------
inline double fWeightGaLe( const double x ) { return 1.0; }
void GaussLegendre( const double x1, const double x2, 
		    std::vector<double>& x, std::vector<double>& w );

//-------------------------------
// GaussHermite: roots and weights for Gauss-Hermite quadrature
//   using orthonormal Hermite polynomials
// Input: x1,x2 (unused) - Gauss-Hermite assumes -infinity - +infinity
// Gauss-Hermite W(x) = exp( -x^2 )
//-------------------------------
inline double fWeightGaHe( const double x ) { return( exp( -x*x ) ); }
void GaussHermite( const double x1, const double x2, 
		   std::vector<double>& x, std::vector<double>& w );

//-------------------------------
// GaussLaguerre: roots and weights for Gauss-Laguerre quadrature
//   using Laguerre polynomials
// Input: x1,x2 (unused) - Gauss-Laguerre assumes 0 - +infinity
// Gauss-Laguerre W(x) = exp(-x)
//-------------------------------
inline double fWeightGaLa( const double x ) { return( exp(-x) ); }
void GaussLaguerre( const double x1, const double x2, 
		   std::vector<double>& x, std::vector<double>& w );


#endif
