//****************************************************************
// Random/tests: test random number generation
// $ ./test_rand <test> <generator> <seed> [test-dependent-params]
//
// Modifications
// 22-Oct-18  fixed bug in testFilling usage printout
// 25-Sep-18  updated for Fall 2018
// 04-Apr-10  add timing to testNonUnifDev
//            add Gauss and GassBM
// 30-Mar-10  add testNonUnifDev
// 16-Mar-10  add chi^2 test to testUniformity
// 07-Mar-10  (re)created
//****************************************************************
using namespace std;

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cpRandom.hpp"
#include "cpStatistics.hpp"
#include "cpFunctions.hpp"
#include "cpQuadrature.hpp"
#include "cpUtils.hpp"

// Function prototypes
cpRandom* makePRNG( int gen_no, int seed, int argc, char* argv[], int& cur_par );
int FindBin( vector<double>& xbins, double x );
int testPeriod( cpRandom* Gen );
int testUniformity( cpRandom* Gen, char* test_par[] );
int testFilling2D( cpRandom* Gen, char* test_par[] );
int testAutoCorr( cpRandom* Gen, char* test_par[] );
int testNonUnifDev( cpRandom* Gen, char* test_par[] );

//================================================================
int main( int argc, char *argv[] )
{
  //-------------------
  // Read command line inputs
  //-------------------
  if( argc < 4 ){
    cout << "Usage: test_rand <test> <generator> <seed> [prng-pars] [test-pars]" 
	 << endl;

    cout << endl;
    cout << "Test=1 : Period of PRNG" << endl;
    cout << "  test-pars: none" << endl;

    cout << endl;
    cout << "Test=2 : Uniformity" << endl;
    cout << "  test-par 1: N(trials)" << endl;
    cout << "  test-par 2: N(bins) for histogram" << endl;
    cout << "  test-par 3: low edge of histogram" << endl;
    cout << "  test-par 4: high edge of histogram" << endl;

    cout << endl;
    cout << "Test=3 : 2D Filling" << endl;
    cout << "  test-par 0: L   - generate an LxL lattice" << endl;
    cout << "  test-par 1: tLo - ave # of pts per site (lower bnd)" << endl;
    cout << "  test-par 2: tHi - ave # of pts per site (upper bnd)" << endl;
    cout << "  test-par 3: N   - steps between tLo and tHi" << endl;

    cout << endl;
    cout << "Test=4 : Auto-correlation of PRNG" << endl;
    cout << "  test-par 1: N(exps)" << endl;
    cout << "  test-par 2: N(meas)   - measurements per exp" << endl;
    cout << "  test-par 3: kLo       - starting lag" << endl;
    cout << "  test-par 4: kHi       - ending lag" << endl;
    cout << "  test-par 5: kPrint    - print exp-by-exp for this k: 0=no print" 
	 << endl;

    cout << endl;
    cout << "Test=5 : Non-Uniform Deviates" << endl;
    cout << "  test-par 1: N points to generate" << endl;
    cout << "  test-par 2: PDF: 0=Exp, 1=Gauss, 2=GaussBM, 3=Poisson" << endl;
    cout << "  test-par 3: PDF: param 1 (lambda,mu,nu)" << endl;
    cout << "  test-par 4: PDF: param 2 (Gauss sigma)" << endl;

    cout << endl;
    cout << "PRNGs:" << endl;
    cout << "  " << idD1A1r << " = " << sPrng[idD1A1r] 
	 << " - prng-pars: none" << endl; 
    cout << "  " << idLocal << " = " << sPrng[idLocal] 
	 << " - prng-pars: none" << endl; 
    cout << "  " << idLCG << " = " << sPrng[idLCG] 
	 << " - prng-pars: a, c, m" << endl; 
    cout << "  " << idDevNrm << " = " << sPrng[idDevNrm] 
	 << " - prng-pars: mean, sigma" << endl; 
    cout << "  " << idMT << " = " << sPrng[idMT] 
	 << " - prng-pars: none" << endl; 

    return -1;
  } // argc < 2

  // Test & PRNG IDs and Random number seed
  int test_no = atoi( argv[1] );   // test
  int gen_no  = atoi( argv[2] );   // PRNG
  int seed    = atoi( argv[3] );   // seed

  // instantiate the PRNG
  int cur_par;
  cpRandom* Gen = makePRNG( gen_no,seed,argc,argv, cur_par );

  //---------------------
  // The Tests (each is a separate function)
  //---------------------
  int stat;
  switch( test_no ){
  case 1:                       // period of PRNG
    stat = testPeriod( Gen );
    break;

  case 2:                       // uniformity
    stat = testUniformity( Gen, &(argv[cur_par]) );
    break;

  case 3:                       // 2D filling
    stat = testFilling2D( Gen, &(argv[cur_par]) );
    break;

  case 4:                       // auto-correlation
    stat = testAutoCorr( Gen, &(argv[cur_par]) );
    break;

  case 5:                       // non-uniform deviates
    stat = testNonUnifDev( Gen, &(argv[cur_par]) );
    break;

  } // switch test_no

  return 0;
} // main

//===============================================================
// makePRNG: instantiates the PRNG specified by gen_no
// Output
//   cur_par: the index of the next parameter in argv[] (after reading PRNG pars)
//===============================================================
cpRandom* makePRNG( int gen_no, int seed, int argc, char* argv[], int& cur_par )
{
  cur_par = 4;               // staring position in argv[] for PRNG parameters
  unsigned int a = 3, c = 4, m = 32;  // default LCG parameters

  cpRandom* Gen;

  // choose PRNG based on gen_no
  switch( gen_no ){
  case idLocal:                                 // Local PRNG
    Gen = new cpRandLocal( (Ullong)seed );
    break;

  case idLCG:                                   // LCG
    if( argc > cur_par ){                         //  read non-default pars
      a = atoT<unsigned int>( argv[cur_par++] );
      c = atoT<unsigned int>( argv[cur_par++] ); 
      m = atoT<unsigned int>( argv[cur_par++] );
    } // argc > cur_par
    Gen = new cpRandLCG( seed,a,c,m );
    break;

  case idMT:                                    // Mersenne Twister
    Gen = new cpRandMT( (Ullong)seed );
    break;

  default:                                      // default D1(A1_r)
    Gen = new cpRandom( (Ullong)seed );
    break;
  } // switch gen_no

  return Gen;
} // makePRNG


//===============================================================
// FindBin: find histogram bin given x-value
//===============================================================
int FindBin( vector<double>& xbins, double x )
{
  int N = (int)xbins.size();
  if( N <= 1 ) return( 0 );

  double dx = xbins[1] - xbins[0];  // assume uniform bin width
  if( x < xbins[0] || x > xbins[N-1]+dx ) // out of range
    return( -1 );

  for( int i = 0; i < N; ++i )  // x is in range
    if( xbins[i] > x ) return( i-1 );

  return( N-1 );
} // FindBin


//===============================================================
// testPeriod: finds period of (32-bit integer) generator
//   this probably only makes sense for LCG generators
//===============================================================
int testPeriod( cpRandom* Gen )
{
  // information about PRNG
  Gen->writeParams();

  // loop until generated number = seed again
  // note: need to trap for case where the PRNG always generates the same #
  Uint next = Gen->getSeed();
  Uint prev = Gen->getSeed();
  Uint period = 0;   // period counter
  do{
    period++;                   // update period counter
    prev = next;                // store previous value
    next = Gen->randInt32();    // generate next number
    cout << "  " << period << " : " << next << endl;
  } while( next != Gen->getSeed() && next != prev );  // do loop
  cout << "Period = " << period << endl;

  return 0;
} // testPeriod

//===============================================================
// testUniformity: constructs 1D uniformity histogram
//===============================================================
int testUniformity( cpRandom *Gen, char* test_par[] )
{
  //----------------------------
  // Test Parameters
  //----------------------------
  int nTrials = atoi( test_par[0] );       // nTrials
  int nBins   = atoi( test_par[1] );       // nBins in output histo
  double lo   = atod( test_par[2] );       // lo and hi edge of bins
  double hi   = atod( test_par[3] );
  double ave  = (double)nTrials / (double)nBins; // expected ave / bin

  //---------------------
  // Generate trials: accumulate histogram
  //---------------------
  double dx = ( hi - lo ) / nBins;
  vector<double> x(nBins);            // x-axis
  for( int i = 0; i < nBins; ++i )
    x[i] = lo + i*dx;                // lower bin edges

  vector<double> y(nBins,0);          // y-axis: fractions
  double next;
  int ind;
  for( int i = 0; i < nTrials; ++i ){
    next = Gen->randDouble();
    if( ( ind=FindBin(x,next) ) >= 0 ) y[ind]++;
  } // i = nTrials

  //------------------------
  // Calculate chi^2 and probability
  // note: Ndof = nBins since already know the average
  //------------------------
  vector<double> pred( nBins,ave );   // all bins should have ave # entries
  int Ndof;
  double chi2  = cpChi2( y,pred,Ndof,0 );
  double Pchi2 = cpProbChi2( chi2,(double)Ndof );

  //------------------------
  // Write out results:
  //   center of bin
  //------------------------
  Gen->writeParams();
  cout << "Bins: " << nBins << "\t Trials: " << nTrials << endl;
  for( int i = 0; i < nBins; ++i )
    cout << x[i] + dx/2 << "\t" << y[i] << endl;

  cout << "chi^2: " << chi2 << endl;
  cout << "P(chi^2;" << Ndof << "): " << scientific
       << Pchi2 << endl;

  return 0;
} // testUniformity


//===============================================================
// testFilling2D:
//   o generates t pairs {x_i x_(i+1)} per site on an LxL lattice
//===============================================================
int testFilling2D( cpRandom *Gen, char* test_par[] )
{
  //----------------------------
  // Test Parameters
  //----------------------------
  int L      = atoi( test_par[0] );       // lattice size (sites)
  double tLo = atod( test_par[1] );       // ave pts per site (lower bound)
  double tHi = atod( test_par[2] );       // ave pts per site (upper bound)
  int nSteps = atoi( test_par[3] );       // # of steps b/w tLo and tHi

  //---------------------
  // construct bin lower bin edges for LxL lattice on x,y = [0-1]
  //---------------------
  double xlo = 0, xhi = 1;
  vector<double> bins(L);
  for( int i = 0; i < L; ++i )
    bins[i] = xlo + i*(xhi-xlo)/L;

  //----------------------
  // construct the lattice and zero it out
  //----------------------
  vector<vector<int> > lattice;   // this is a 2D vector = matrix
  vector<int> row(L,0);
  for( int i = 0; i < L; ++i ) lattice.push_back( row );

  Gen->writeParams();
  cout << "2D-Lattice-sites: " << L << " Steps: " << nSteps << endl;

  double x,y;
  int ix,iy;
  int N;
  vector<double> nEmpty;  // empty sites vs t
  vector<double> pred;    // pred of empty sites vs t (from exp)
  int iEmpty = 0;
  double dt = ( tHi - tLo ) / nSteps;
  //----------------------
  // Loop over values of t requested
  //----------------------
  for( double t = tLo; t <= tHi; t += dt ){
    N = (int)(t*L*L);              // number of pts to generate

    for( int i = 0; i < L; ++i ){  // zero out lattice
      for( int j = 0; j < L; ++j ) lattice[i][j] = 0;
    }

    //----------------------
    // Generate the t*L^2 x,y pairs and fill lattice
    // o sites contain 1 if a pair generated there
    //                 0 otherwise
    //----------------------
    for( int i = 0; i < N; ++i ){
      x = Gen->randDouble();  ix = FindBin( bins,x );
      y = Gen->randDouble();  iy = FindBin( bins,y );
      if( ix >= 0 && iy >= 0 ) lattice[iy][ix] = 1;
    } // i = trials

    //--------------------
    // Print out filled lattice and count of empty sites
    //--------------------
    nEmpty.push_back( 0 );
    for( int i = 0; i < L; ++i ){
      for( int j = 0; j < L; ++j ){
	if( lattice[i][j] == 0 ) nEmpty[iEmpty] += 1.0;
	//cout << lattice[i][j] << " ";
      } // j = lattice sites
      //cout << endl;
    } // i = lattice sites

    // write out result
    cout << "t: " << t << " \t Empty-sites: " << nEmpty[iEmpty]
	 << " \t Empty-fraction: " << nEmpty[iEmpty]/(L*L) << endl;

    //------------------------------
    // Calculate prediction from exponential decay for bin from
    //   [t,t+1]
    //------------------------------
    pred.push_back( L*L * exp(-(double)(t)) );

    iEmpty++;
  } // t-loop

  //-----------------------------
  // Now calculate chi^2 b/w empty sites vs t: data and pred
  //-----------------------------
  int Ndof;
  double chi2  = cpChi2( nEmpty,pred,Ndof,0 );
  double Pchi2 = cpProbChi2( chi2,Ndof );
  cout << "Chi^2: " << chi2 << "  P(chi2;" << Ndof << "): " << Pchi2 << endl;

  return 0;
} // testFilling2D


//===============================================================
// testAutoCorr: tests auto-correlation
//===============================================================
int testAutoCorr( cpRandom *Gen, char* test_par[] )
{
  //----------------------------
  // Test Parameters
  //----------------------------
  int nExp   = atoi( test_par[0] );       // # of exp's
  int nMeas = atoi( test_par[1] );        // measurements per exp
  int kLo    = atoi( test_par[2] );       // starting lag
  int kHi    = atoi( test_par[3] );       // ending lag
  int kPrint = atoi( test_par[4] );       // lag to print

  Gen->writeParams();
  cout << "N(exp): " << nExp << "  N(trial): " << nMeas 
       << "  k(lo): " << kLo << "  k(hi): " << kHi 
       << "  kPrint: " << kPrint
       << endl;

  vector<double> seq;
  for( int i = 0; i < nExp; ++i ){
    for( int j = 0; j < nMeas; ++j )
      seq.push_back( (double)Gen->randInt32() );
  } // i-loop

  //------------------------
  // print out all values of C(k) for selected kPrint
  //------------------------
  double CkPrint;
  if( kPrint > 0 ){
    for( int i = 0; i < nExp; ++i ){
      CkPrint = cpAutoCorr( seq,kPrint,i*nMeas,nMeas );
      cout << "  Exp: " << i 
	   << "  C(" << kPrint << ") = " << CkPrint << endl;
    } // i loop
  } // kPrint > 0

  //-------------------------
  // calculate auto-correlation for various values of k
  //-------------------------
  vector<double> Ck( nExp );
  double aveCk,errCk;
  for( int k = kLo; k <= kHi; ++k ){

    // C(k) for each experiment
    for( int i = 0; i < nExp; ++i )
      Ck[i] = cpAutoCorr( seq,k,i*nMeas,nMeas );

    // now calculate average for this lag
    aveCk = cpAveVar( Ck,errCk );
    errCk = sqrt( errCk / nExp );  // error on mean
    cout << "Lag: " << noshowpos << fixed << setw(3) << k 
	 << showpos << scientific << setprecision(5) 
	 << " <C(k)>: " << aveCk << "  Sig[C(k)]: " << errCk << endl;

  } // k-loop

  return 0;
} // testAutoCorr


//===============================================================
// testNonUnifDev: generate non-uniform deviate from specific PDFs
//===============================================================
int testNonUnifDev( cpRandom *Gen, char* test_par[] )
{
  const int idExp = 0;
  const int idGau = 1;
  const int idGbm = 2;
  const int idPoi = 3;
  string sModel[4] = { "Exponential","Gauss:Ratio-of-Uniforms",
		       "Gauss:Box-Muller", "Poisson" };

  //----------------------------
  // Test Parameters
  //----------------------------
  int nPts    = atoi( test_par[0] );       // # of points to generate
  int model   = atoi( test_par[1] );       // PDF ID (1=Exp, 2=Gauss, 3=Poiss)
  double mean = atod( test_par[2] );       // lambda, mu, or nu
  double sig  = 0;
  if( model == idGau || model == idGbm )
    sig       = atod( test_par[3] );       // Gaussian sigma

  //----------------------------
  // Set up Histogram parameters
  //----------------------------
  int nBins = 100;
  double xlo,xhi;
  switch( model ){
  case idExp:
    xlo = 0;             xhi = 10*mean;       break;
  case idGau: case idGbm:
    xlo = mean - 5*sig;  xhi = mean + 5*sig;  break;
  case idPoi:
    nBins = 5 * mean;    // Poisson distrib takes integers
    xlo = 0;             xhi = (double)nBins; break;
  } // switch

  // create histogram bins
  double dx = ( xhi - xlo ) / nBins;
  vector<double> xbin( nBins );      // x axis
  vector<double> ybin( nBins,0 );    // N_i / nPts
  vector<double> dybin( nBins,0 );   // error on ybin
  vector<double> ypred( nBins,0 );   // prediction for this bin
  for( int i = 0; i < nBins; ++i )
    xbin[i] = xlo + i*dx;

  // parameters for prediction
  vector<double> par;
  par.push_back( mean );
  if( model == idGau || model == idGbm) par.push_back( sig );

  //---------------------
  // Generate points - fill histogram
  //---------------------
  double x;
  int ix;
  cpTimer* cpT = new cpTimer();            // time the generation
  for( int ipt = 0; ipt < nPts; ++ipt ){

    // generate point based on which model we've chosen
    switch( model ){
    case idExp:
      x = Gen->randDevExp( mean );  break;
    case idGau:
      x = Gen->randDevNrm( mean,sig );  break;
    case idGbm:
      x = Gen->randDevNrmBM( mean,sig );  break;
    case idPoi:
      x = Gen->randDevPoisson( mean );  break;
    } // switch

    // fill bin corresponding to generated x
    ix = FindBin( xbin,x );
    if( ix >= 0 ) ybin[ix] += 1;
    
  } // ipt
  double tGen = cpT->elapsed();   // elapsed time

  //----------------------
  // Normalize histogram to number of events
  // Calculate Poisson errors
  // Calculate Prediction = Gauss-Legendre integral over bin
  // Print results
  //----------------------
  cout << sModel[model] << "  N(points): " << nPts 
       << "  N(bins): " << nBins << " : " << xlo << " - " << xhi << endl;
  double xcen;
  int nStep;
  for( int i = 0; i < nBins; ++i ){
    ybin[i]  = ybin[i] / nPts;                  // normalize
    dybin[i] = cpBinomErr( ybin[i], nPts );     // Poisson errors
    xcen     = xbin[i] + dx/2;                  // bin center
    switch( model ){                            // pred at bin center
    case idExp:
      ypred[i] = qGauss( nStep,fExp, par, xbin[i], xbin[i]+dx );
      break;
    case idGau: case idGbm:
      ypred[i] = qGauss( nStep, fGauss, par, xbin[i], xbin[i]+dx );
      break;
    case idPoi:
      xcen = xbin[i];  // Poisson uses only integers - change center
      ypred[i] = fPoisson( xcen, par ); // prediction is for this integer
      break;
    } // switch

    cout << fixed << setprecision(5) << xcen << scientific 
	 << "  " << ybin[i] << "  " << dybin[i] 
	 << "  " << ypred[i] << endl;
  } // i = bins

  //-----------------------
  // Now calculate chi^2 b/w generated points and prediction
  //   use 1 constraint b/c normalizing to total number of events
  //-----------------------
  int Ndof;
  double chi2 = cpChi2( ybin,ypred,dybin,Ndof,1 );
  cout << "chi^2: " << fixed << setprecision(4) << chi2 << " / " << Ndof
       << "  Probability = " << scientific << cpProbChi2( chi2, (double)Ndof )
       << endl;

  //----------------------
  // Print out the time
  //----------------------
  cout << "Time: " << fixed << tGen << " us" << endl;

  return 0;
} // testNonUnifDev
