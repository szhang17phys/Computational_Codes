//***********************************************************************
// Run tests of numerical integration
//
// $ ./tests <test> <function-number> <xlo> <xhi> [params...]
//
// function-number: see function list below
// xlo, xhi:        limits of integration
//
// Test = 1: Convergence vs N(steps) for Rect. Trap. Simpson Romberg
// Test = 2: Resolution Study
// Test = 3: Gaussian Quadrature roots & weights
// Test = 4: Integral with specific method
//
// Modifications
// 13-Aug-18  update for Fall 2018
// 19-Feb-10  add Gauss-Laguerre, testResAlt
// 16-Feb-10  add testRootWeight(), testIntegral()
//            add Gauss-Hermite (idGaHe)
// 08-Feb-10  add testResolution
//***********************************************************************
using namespace std;

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>

#include "cpQuadrature.hpp"
#include "cpUtils.hpp"
#include "cpFunctions.hpp"

// Function Prototypes
inline double dev( double x1, double x2 ){ 
  return( (x2 == 0.0) ? x1 : x1/x2 - 1.0 ); }
int testConvvsN( cp1Dfunction f, double a, double b, int argc, char* argv[] );
int testResolution( cp1Dfunction f, double a, double b, int argc, char* argv[] );
int testResAlt( cp1Dfunction f, double a, double b, int argc, char* argv[] );
int testRootWeight( cp1Dfunction f, double a, double b, int argc, char* argv[] );
int testIntegral( cp1Dfunction f, double a, double b, int argc, char* argv[] );

//================================================================
int main( int argc, char* argv[] )
{
  //---------------------
  // Set up an array of functions from cpFunctions
  //---------------------
  cp1Dfunction Func;
  string sFunc;

  //---------------------
  // Read command line inputs: test number
  //---------------------
  if( argc < 4 ){
    cout << "Quadrature: ./testQuadrature <test> <function-id> <xlo> <xhi> [params]" 
	 << endl;
    cout << endl;
    cout << "Test=1 : Convergence vs N(steps)" << endl;
    cout << "  param 1     exact integral" << endl;
    cout << "  param 2     log_2(nSteps-max)" << endl;
    cout << "  param 3...  function parameters" << endl;
    cout << endl;
    cout << "Test=2 : Resolution" << endl;
    cout << "  param 1     method: "
	 << idTrap << "=Trap " << idSimp << "=Simp " << idRomb << "=Romb "
	 << idGaLa << "=Gauss-Laguerre "
	 << endl;
    cout << "  param 2     Convergence [def=" << CONV << "]" << endl;
    cout << "  param 3     N(bins) [def=50]" << endl;
    cout << "  param 4     sigma [def=0.25/lambda]" << endl;
    cout << "  param 5     Plot tlo [def=-1/lambda]" << endl;
    cout << "  param 6     Plot thi [def=4/lambda]" << endl;
    cout << endl;
    cout << "Test=3 : Gaussian Quadrature Roots & Weights" << endl;
    cout << "  param 1     method: "
	 << idGaLe << "=Gauss-Legendre " << idGaHe << "=Gauss-Hermite "
	 << idGaLa << "=Gauss-Laguerre "
	 << endl;
    cout << "  param 2     Number of points" << endl;
    cout << endl;
    cout << "Test=4 : Integrate by Method" << endl;
    cout << "  param 1     method: "
	 << idTrap << "=Trap " << idSimp << "=Simp " << idRomb << "=Romb "
	 << idGaLe << "=Gauss-Legendre " << idGaHe << "=Gauss-Hermite "
	 << idGaLa << "=Gauss-Laguerre "
	 << endl;
    cout << "  param 2     exact integral" << endl;
    cout << "  param 3     Convergence [def=" << CONV << "]" << endl;
    cout << "  param 4...  function parameters" << endl;
    cout << endl;
    cout << "Test=5 : Resolution (Gauss Quadrature)" << endl;
    cout << "  param 1     Convergence [def=" << CONV << "]" << endl;
    cout << "  param 2     N(bins) [def=50]" << endl;
    cout << "  param 3     sigma [def=0.25/lambda]" << endl;
    cout << "  param 4     Plot tlo [def=-1/lambda]" << endl;
    cout << "  param 5     Plot thi [def=4/lambda]" << endl;
    cout << endl;
    cout << "Functions:  " << endl;
    for( int i = 0; i < n1Dfuncs; ++i ){
      Func = get1Dfunction( i, sFunc );
      cout << "  " << i << " : " << sFunc << endl;
    } // i
    return -1;
  } // argc < 2
  int test_no = atoi( argv[1] );
  int func_id = atoi( argv[2] );
  double xlo  = atod( argv[3] );
  double xhi  = atod( argv[4] );

  //--------------------------
  // Get Function from list in cpFunctions
  //---------------------------
  Func = get1Dfunction( func_id, sFunc );

  //---------------------
  // The Tests (each is a separate function)
  //---------------------
  int stat;
  switch( test_no ){

  case 1:                       // convergence vs N(steps)
    cout << "CONVERGENCE COMPARISON: " << sFunc << endl;
    stat = testConvvsN( Func,xlo,xhi,argc,argv );
    break;

  case 2:                       // experimental Resolution
    cout << "EXPERIMENTAL RESOLUTION: " << sFunc << endl;
    stat = testResolution( Func,xlo,xhi,argc,argv );
    break;

  case 3:                       // Gaussian Quadrature
    cout << "ROOTS & WEIGHTS: " << endl;
    stat = testRootWeight( Func,xlo,xhi,argc,argv );
    break;

  case 4:                       // General Integral
    cout << "INTEGRATION: " << sFunc << endl;
    stat = testIntegral( Func,xlo,xhi,argc,argv );
    break;

  case 5:                       // experimental Resolution (alt)
    cout << "EXP RES (GAUSS QUAD): " << sFunc << endl;
    stat = testResAlt( Func,xlo,xhi,argc,argv );
    break;

  default:                      // bad test number
    cout << "Invalid test number: " << test_no << endl;
    break;

  } // switch test_no

  return 0;
} // main

//================================================================
// testConvvsN: integral vs N(step) for various methods
// Parameters
//   argv[5]    = analytic integral value
//   argv[6]    = N(step) ==> 2^N(step) intervals
//   argv[7...] = function params
// Output
//   for i = 0, N
//   <2^i> <Rectangular integral> <Trapezoid integral> <Simpson integral>
//   <2^M> <Romberg integral>
//================================================================
int testConvvsN( cp1Dfunction f, double a, double b, int argc, char* argv[] )
{
  //-------------------
  // Read parameters from command line
  //-------------------
  double exactInt = 0.0;     // exact integral
  if( argc > 5 ) exactInt = atod( argv[5] );

  int log2nStep = 10;     // log_2(Number of steps)
  if( argc > 6 ) log2nStep = atoi( argv[6] );

  vector<double> pars;             // function params
  if( argc > 7 )
    for( int i = 7; i < argc; ++i ) pars.push_back( atod( argv[i] ) );

  //---------------------
  // Loop over steps, printing out integral estimate each time
  //---------------------
  int nR,nT,nS;
  double qR,qT,qS;
  double dR,dT,dS;
  int prec = 9;
  cout << left << setw(4) << log2nStep 
       << "intervals "
       << "Rectangle                         | "
       << "Trapezoid                         | "
       << "Simpson" << endl;
  for( int i = 0; i < log2nStep; ++i ){
    qR = qRect( nR,f,pars,a,b,i );    dR = dev( qR,exactInt );
    qT = qTrap( nT,f,pars,a,b,0,i );  dT = dev( qT,exactInt );
    qS = qSimp( nS,f,pars,a,b,0,i );  dS = dev( qS,exactInt );

    cout.unsetf( ios_base::scientific );
    cout.unsetf( ios_base::showpos );
    cout << left << setw(4) << i << " " << setw( 8 ) << (int)pow(2.0,i) 
	 << showpos << scientific << setprecision(prec) 
	 << " " << qR << " " << dR
	 << " | " << qT << " " << dT
	 << " | " << qS << " " << dS << endl;
  } // i = log2nStep loop

  // Romberg integration
  int nRo;
  double qRo = qRomb( nRo,f,pars,a,b );
  double dRo = dev( qRo,exactInt );
  cout.unsetf( ios_base::scientific );
  cout.unsetf( ios_base::showpos );
  cout << left << setw(4) << nRo << " " << setw( 8 ) << (int)pow(2.0,nRo) 
       << scientific << setprecision(prec) 
       << " " << qRo << " " << dRo << " Romberg" << endl;

  return 0;
} // testConvvsN

//================================================================
// testResolution: integral vs N(step) for various methods
// Parameters
//   argv[5]  = integration method
//   argv[6]  = convergence
//   argv[7]  = N(bins)
//   argv[8]  = sigma
//   argv[9]  = plot tlo
//   argv[10] = plot thi
// Output
//   convolution integral for bins of x^{meas}
//================================================================
int testResolution( cp1Dfunction f, double a, double b, int argc, char* argv[] )
{
  //-------------------
  // Read parameters from command line
  //-------------------
  double lambda = 1.0;    // plot everything per unit lifetime (tau=1/lambda=1)

  int idMethod = 0;       // integration method
  if( argc > 5 ) idMethod = atoi( argv[5] );

  double converge = CONV;                     // convergence condition
  if( argc > 6 ) converge = atod( argv[6] );

  int nBins = 50;         // number of output bins
  if( argc > 7 ) nBins = atoi( argv[7] );

  double sigma = 0.25 / lambda;    // Gaussian sigma
  if( argc > 8 ) sigma = atod( argv[8] ) / lambda;

  double tlo = -1.0/lambda, thi = 4.0/lambda;  // limits of output plot
  if( argc > 9 )  tlo = atod( argv[9] )  / lambda;
  if( argc > 10 ) thi = atod( argv[10] ) / lambda;

  double dt = ( thi - tlo ) / (double)nBins;  // bin size

  //---------------------
  // Parameters for convolution function
  //---------------------
  vector<double> pars(3);
  pars[0] = lambda;  // lambda
  pars[1] = 0.0;     // t^{true} will be set in loop over bins
  pars[2] = sigma;   // sigma

  //---------------------
  // Loop over each t^{meas} bin doing integral over t^{true}
  //---------------------
  double q;
  int nSteps;
  cout << "Plot-Region: " << tlo << " " << thi << " : " << nBins << " bins" 
       << endl;
  for( int i = 0; i < nBins; ++i ){
    pars[1] = tlo + (i+0.5)*dt;     // t^{meas} at center of bin

    // do the integral for this bin using selected method
    switch( idMethod ){
    case idTrap:
      q = qTrap( nSteps,f,pars,a,b,converge ); break;
    case idSimp:
      q = qSimp( nSteps,f,pars,a,b,converge ); break;
    case idRomb:
      q = qRomb( nSteps,f,pars,a,b,converge ); break;
    case idGaLa:
      q = qGauss( nSteps,f,pars,a,b,idMethod,converge ); break;
    } // switch

    cout << scientific << setprecision(9) 
	 << pars[1] << " \t" << q << " "
	 << fixed << setw(10) << nSteps << endl;
  } // i = nBins loop

  return 0;
} // testResolution

//================================================================
// testResolution: integral vs N(step) for various methods
//   should use f = fEcGalt
// Parameters
//   argv[5] = convergence
//   argv[6]  = N(bins)
//   argv[7]  = sigma
//   argv[8]  = plot tlo
//   argv[9]  = plot thi
// Output
//   convolution integral for bins of x^{meas}
//================================================================
int testResAlt( cp1Dfunction f, double a, double b, int argc, char* argv[] )
{
  if( f != fEcGalt ) throw( "testResAlt: use function fEcGalt" );

  //-------------------
  // Read parameters from command line
  //-------------------
  double lambda = 1.0;    // plot everything per unit lifetime (tau=1/lambda=1)

  double converge = CONV;                     // convergence condition
  if( argc > 5 ) converge = atod( argv[5] );

  int nBins = 50;         // number of output bins
  if( argc > 6 ) nBins = atoi( argv[6] );

  double sigma = 0.25 / lambda;    // Gaussian sigma
  if( argc > 7 ) sigma = atod( argv[7] ) / lambda;

  double tlo = -1.0/lambda, thi = 4.0/lambda;  // limits of output plot
  if( argc > 8 ) tlo = atod( argv[8] ) / lambda;
  if( argc > 9 ) thi = atod( argv[9] ) / lambda;

  double dt = ( thi - tlo ) / (double)nBins;  // bin size

  //---------------------
  // Parameters for convolution function (fEcG)
  //  [1/(sqrt(2pi)lambda sigma)] exp[ -x^2/2 - lambda sigma x ]
  //---------------------
  vector<double> pars( 1, lambda*sigma );

  //---------------------
  // Loop over each t^{meas} bin doing integral over t^{true}
  //---------------------
  double tm,qInf,qCor,Nrm;
  int nSInf,nSCor;
  cout << "Plot-Region: " << tlo << " " << thi << " : " << nBins << " bins" 
       << endl;
  for( int i = 0; i < nBins; ++i ){
    tm  = tlo + (i+0.5)*dt;             // t^{meas} at center of bin

    // pre-factor in front of integrals
    Nrm = exp( -lambda * tm ) / ( sqrt(2.0*cpPI) * lambda );

    // do [0,infinity] integral using Gauss-Laguerre
    qInf = qGauss( nSInf,f,pars,0,0,idGaLa,converge );

    // correction integral using Gauss-Legendre
    // subtract off [0, -t^m/sigma] if t^m < 0
    if( tm < 0.0 )
      qCor = -qGauss( nSCor,f,pars,0.0,-tm/sigma,idGaLe,converge );
    // add [-t^m/sigma, 0] part if t^m >= 0
    else
      qCor = qGauss( nSCor,f,pars,-tm/sigma,0.0,idGaLe,converge );

    // print out results
    cout << scientific << setprecision(9) 
	 << tm << " \t" << Nrm * (qInf + qCor) << " "
	 << fixed << setw(10) << nSInf + nSCor << endl;
    //	 << " " << qInf << " " << qCor << endl;
  } // i = nBins loop

  return 0;
} // testResAlt

//================================================================
// testRootWeight: prints out Gaussian Roots and Weights for a given method
// Parameters
//   argv[5]    = Gaussian Quadrature method
//   argv[6]    = N(points)
// Output
//   root and weight for each point
//================================================================
int testRootWeight( cp1Dfunction f, double a, double b, int argc, char* argv[] )
{
  //-------------------
  // Read parameters from command line
  //-------------------
  int idMethod = idGaHe;       // integration method
  if( argc > 5 ) idMethod = atoi( argv[5] );

  int nPts = 1;                // number of points
  if( argc > 6 ) nPts = atoi( argv[6] );

  // vectors that will contain the roots and weights
  vector<double> x( nPts );
  vector<double> w( nPts );

  // choose the method
  switch( idMethod ){
  case idGaLe:
    GaussLegendre( a,b,x,w );   break;
  case idGaHe:
    GaussHermite( a,b,x,w );   break;
  case idGaLa:
    GaussLaguerre( a,b,x,w );   break;
  } // switch

  // print out
  cout << sMethod[idMethod] << " : " << nPts << " points" << endl;
  for( int i = 0; i < nPts; ++i ){
    cout.unsetf( ios_base::showpos );
    cout << setw( 4 ) << i << " "
	 << showpos << fixed << setprecision( 10 ) << x[i] << "  " << w[i]
	 << endl;
  } // i loop

  return 0;
} // testRootWeight

//================================================================
// testIntegral: integrates a function using a specified method
// Parameters
//   argv[5]    = integration method
//   argv[6]    = analytic integral value
//   argv[7]    = convergence parameter
//   argv[8...] = function params
// Output
//================================================================
int testIntegral( cp1Dfunction f, double a, double b, int argc, char* argv[] )
{
  //-------------------
  // Read parameters from command line
  //-------------------
  int idMethod = idGaLe;       // integration method
  if( argc > 5 ) idMethod = atoi( argv[5] );
  
  double exactInt = 0.0;     // exact integral
  if( argc > 6 ) exactInt = atod( argv[6] );

  double converge = CONV;    // convergence condition
  if( argc > 7 ) converge = atod( argv[7] );

  vector<double> pars;             // function params
  if( argc > 8 )
    for( int i = 8; i < argc; ++i ) pars.push_back( atod( argv[i] ) );

  // do the integral depending on method chosen
  int nSteps;
  double q;
  switch( idMethod ){
  case idTrap:
    q = qTrap( nSteps,f,pars,a,b,converge ); break;
  case idSimp:
    q = qSimp( nSteps,f,pars,a,b,converge ); break;
  case idRomb:
    q = qRomb( nSteps,f,pars,a,b,converge ); break;
  case idGaLe: case idGaHe: case idGaLa:
    q = qGauss( nSteps,f,pars,a,b,idMethod,converge ); break;
  } // switch

  double dev = (exactInt != 0.0) ? q/exactInt - 1.0 : q;
  cout << sMethod[idMethod] << ": "
       << scientific << setprecision(9) 
       << q << "  " <<  dev << "  " << nSteps << endl;

  return 0;
} // testIntegral
