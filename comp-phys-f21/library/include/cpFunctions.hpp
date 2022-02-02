#ifndef CPFUNCTIONS_H
#define CPFUNCTIONS_H

//**********************************************************************
// cpFunctions: compilation of testing, mathematical and other functions
//   mainly from Numerical Recipes
//
// Modifications
// 12-Apr-10  add cpInvGammaP
// 30-Mar-10  add fPoisson
//            fix incorrect normalization in fExp
// 15-Mar-10  add Gamma functions + <limits>
// 19-Feb-10  add fLaguerre
// 17-Feb-10  add fLegendre
// 15-Feb-10  add fHermite, fExp2,fExp2x2
// 06-Feb-10  add fGauss, fExp
// 01-Feb-10  add function IDs, get1Dfunction
//            add fSin_2pix, fS2_2pix
// 15-Dec-09  move constants cpPI and cpSqrt2PI to cpUtils
//**********************************************************************

#include <limits>
#include <cmath>
#include <vector>
#include <string>

#include "cpUtils.hpp"

//-------------------------------------------------------------------
// Testing Functions
//   note: each function has an ID and a descriptor string
//     they can be accessed by id using get1Dfunction
//   also note: all functions with code in the .hpp file must be
//     declared inline (for the Mac compiler)
//------------------------------------------------------------------

//-------------------------------
// typedef specifying form of 1D functions
//   arg 1 = variable (x)
//   arg 2 = vector of function parameters
//-------------------------------
typedef double (*cp1Dfunction)( double, const std::vector<double>& );

const int n1Dfuncs = 11;     // change this every time a new function is added
cp1Dfunction get1Dfunction( const int id, std::string& name );

// Degree N polynomial - params: p_0, p_1, p_2,... (degree from size of vector)
const int idPolyN = 0;
const std::string sPolyN( "O(N) Polynomial" );
double fPolyN( const double x, const std::vector<double>& p );

// sqrt( 1 - x^2 )
const int idSqrt1mX = 1;
const std::string sSqrt1mX( "sqrt(1 - x^2)" );
inline double fSqrt1mX( const double x, const std::vector<double>& p ){
  return( sqrt( 1.0 - x*x ) ); }

// sin(2pi x) - no params
const int idSin_2pix = 2;
const std::string sSin_2pix( "2#pi sin(2#pi x)" );
inline double fSin_2pix( const double x, const std::vector<double>& p ){
  return( 2.0*cpPI * sin( 2.0*cpPI*x ) ); }

// sin^2(2pi x) - no params
const int idS2_2pix = 3;
const std::string sS2_2pix( "sin^2(2#pi x)" );
inline double fS2_2pix( const double x, const std::vector<double>& p ){
  return( sin( 2.0*cpPI*x ) * sin( 2.0*cpPI*x ) ); }

// Gaussian: p[0] = mu, p[1] = sigma
const int idGauss = 4;
const std::string sGauss( "G( x ; #mu,#sigma )" );
inline double fGauss( const double x, const std::vector<double>& p ){
  return( exp( -(x-p[0])*(x-p[0]) / (2.0*p[1]*p[1]) ) 
	  / ( sqrt(2.0*cpPI)*p[1] ) ); }

// Exponential Decay: p[0] = lambda
const int idExp = 5;
const std::string sExp( "#lambda exp( -#lambda x )" );
inline double fExp( const double x, const std::vector<double>& p ){
  return( (x >= 0.0) ? ( p[0] * exp( -x * p[0] ) ) : 0.0 ); }

// Exponential Decay conv. w/ Gaussian
// p[0] = lambda ; p[1] = x^{meas} (mu in Gauss) ; p[2] = sigma of Gauss
// x = x^{true} (the integration variable)
const int idEcG = 6;
const std::string sEcG( "Exp conv. Gauss" );
inline double fEcG( const double x, const std::vector<double>& p ){
  std::vector<double> p1( p.begin()+1, p.end() );
  return( fExp( x,p ) * fGauss( x,p1  ) ); }

// exp[ -x^2 ]
const int idExp2 = 7;
const std::string sExp2( "exp(-x^{2}) / #sqrt{#pi}" );
inline double fExp2( const double x, const std::vector<double>& p ){
  return( exp(-x*x) / sqrt(cpPI) ); }

// x^2 exp[ -x^2 ]
const int idExp2x2 = 8;
const std::string sExp2x2( "2 x^{2} exp(-x^{2}) / #sqrt(#pi)" );
inline double fExp2x2( const double x, const std::vector<double>& p ){
  return( 2.0 * x*x * exp(-x*x) / sqrt(cpPI) ); }

// Exponential Decay conv. w/ Gaussian (alt version)
// fEcGalt: alternate Exp conv Gauss
//   f    = exp[ -x^2/2 - lambda sigma x ]
//   p[0] = lambda * sigma
// note: need to multiply this by the normalization
//   exp(-lambda x^m) / sqrt(2pi) lambda
const int idEcGalt = 9;
const std::string sEcGalt( "Exp conv. Gauss (alt)" );
inline double fEcGalt( const double x, const std::vector<double>& p ){
  return( exp(-x*x/2 - p[0]*x) ); }

// Poisson Distribution
// f(k;nu) = nu^{k} exp(-nu) / Gamma(k+1)
// p[0] = nu
const int idPoisson = 10;
const std::string sPoisson( "Poisson: nu^{k} e^{-k} / k!" );
double fPoisson( const double x, const std::vector<double>& p );


//------------------------------------------------------------------
// Mathematical Functions
//------------------------------------------------------------------
// Orthogonal Polynomials
double fLegendre( const int n, const double x, double& dPdx );
double fHermite( const int n, const double x, double& dHdx );
double fLaguerre( const int n, const double x, double& dPdx );

//------------------------------------------------------------------
// Gamma Functions
//   the user should call cpGammaP, cpGammaQ, cpInvGammaP
//------------------------------------------------------------------
// constants for Gamma Functions
const double Gamma_EPS   = std::numeric_limits<double>:: epsilon();
const double Gamma_FPMIN =  std::numeric_limits<double>:: min() / Gamma_EPS;

// Incomplete Gamma Function
// P(a,x) = Int_0^x e^(-t) t^(a-1) dt / Gamma(a)
double cpGammaP( const double a, const double x );

// Complemenntary Incomplete Gamma Function
// Q(a,x) = 1 - P(a,x) - Int_x^infinity e^(-t) t^(a-1) dt / Gamma(a)
double cpGammaQ( const double a, const double x );

// Inverse Incomplete Gamma Function
// solves P(a,x) = p for x
double cpInvGammaP( const double p, const double a );

// internal functions
double gammln( const double a );    // ln[Gamma(a)]
double gser( const double a, const double x );  // series rep. of P(a,x)
double gcf( const double a, const double x );   // continued fraction rep P(a,x)

//------------------------------------------------------------------
// Sorting Functions
//------------------------------------------------------------------
template<class T>
void sortSI( std::vector<T>& arr );   // Straight Insertion sort (piksrt in NR)

#endif
