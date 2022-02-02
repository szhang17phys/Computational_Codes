#ifndef CPSTATISTICS_H
#define CPSTATISTICS_H

//**********************************************************************
// cpStatistics: Statistical functions
//   mainly from Numerical Recipes
//
// Functions
//   cpAveVar       average and variance
//   cpWgtAve       weighted average
//   cpChi2         chi^2 calculations
//   cpProbChi2     chi^2 probability
//   cpInvProbChi2  inverse chi^2 probability (x s.t. P(x,nu) = p)
//   cpAutoCorr     auto-correlation
//   cpBinomErr     binomial error calculation
//   cpLinReg       least squares fit to linear mode
//
// Modifications
// 27-Oct-10  add cpInvProbChi2
// 09-Apr-10  add cpLinReg
// 07-Apr-10  add cpWgtAve
// 30-Mar-10  add cpBinomErr
//            modify least-squares cpChi2 to ignore bins w/ sigma=0
// 23-Mar-10  fix bug in cpChi2(Pearson) - data=0 vs pred=0
// 13-Mar-10  created
//**********************************************************************

#include <cmath>
#include <vector>
#include <string>
#include <limits>   // have to include b/c of cpFunctions - but why???

#include "cpUtils.hpp"
#include "cpFunctions.hpp"
#include "cpMatrix.hpp"

//-------------------------------------------------------------------
// cpAveVar: calculates average and variance of a set of measurements
// Inputs:
//   data  vector of measurements
//   start average entries in data vector starting from index start
//   N     number of data points to average
// Outputs:
//   var   variance (sigma^2)
// Returns: average of data
//------------------------------------------------------------------
double cpAveVar( const std::vector<double>& data, double& var, int start=0, int N=0 );

//-------------------------------------------------------------------
// cpWgtAve: calculates weighted average of a set of measurements
// Inputs:
//   data  vector of measurements
//   sigma vector of errors (uncorrelated)
// Outputs:
//   err   uncertainty on average
//   chi2  chi-squared of average
//   Ndof  number of degrees of freedom in average
// Returns: weighted average of data
//------------------------------------------------------------------
double cpWgtAve( const std::vector<double>& data, const std::vector<double>& sigma,
		 double& err, double& chi2, int& Ndof );

//------------------------------------------------------------------
// cpChi2: calculates chi^2 in a variety of ways
//   Pearson's chi^2 statistic (for Poisson prediction)
//     chi^2 = Sum_i=0,N-1 (data_i - pred_i)^2 / pred_i
//   Least Squares (for data with associated errors)
//     chi^2 = Sum_i=0,N-1 (data_i - pred_i)^2 / sigma^2
// Inputs
//   meas    vector of measurements
//   pred    vector of predictions
//   sig     vector of errors on measurements (Least Squares version only)
//   Nconstr number of constraints, =1 if e.g. pred normalized to meas
// Outputs
//   Ndof   degrees of freedom
// Returns: chi-squared
//-------------------------------------------------------------------
// Pearson's chi^2 version
double cpChi2( const std::vector<double>& meas, const std::vector<double>& pred,
	       int& Ndof, const int Nconstr=1 );
// Least Squares version
double cpChi2( const std::vector<double>& meas, const std::vector<double>& pred,
	       const std::vector<double>& sig, int& Ndof, const int Nconstr=1 );

//------------------------------------------------------------------
// cpChi2Prob: calculates chi^2 probability
//   P(chi^2; nu) = Q( nu/2, chi^2/2 )
//     = probability of observing a chi^2 greater than chi^2 
//       given nu degrees of freedom
// Inputs:
//   chi2  chisquared
//   nu    number of degrees of freedom
//------------------------------------------------------------------
inline double cpProbChi2( const double chi2, const double nu ){
  return( (chi2>=0 && nu>0) ? cpGammaQ( 0.5*nu, 0.5*chi2 ) : -1 ); }

//------------------------------------------------------------------
// cpInvProbChi2: calculates inverse chi^2 probability
//   returns x such that P(x; nu) = Q( nu/2, x/2 ) = p
// Note: since this function uses:
//   p = P_chi(x,nu) = Q(nu/2,x/2) = 1 - P(nu/2,x/2)
// we need to invert
//   1-p = P(nu/2,x/2)  to get  x/2
// Inputs:
//   prob  probability
//   nu    number of degrees of freedom
//------------------------------------------------------------------
inline double cpInvProbChi2( const double prob, const double nu ){
  return( (prob>=0 && prob <=1 && nu>0) ? 
    2.*cpInvGammaP( (1.0-prob), nu/2.0 ) : -1 ); }

//------------------------------------------------------------------
// cpAutoCorr: calculates auto-correlation
//   C(k) = [<x_{i} x_{i+k}> - <x_i>^2] / [ <x_i^2> - <x_i>^2]
// Inputs
//   x      vector with data sequence
//   k      lag
//   start  starting position in vector for calculation
//   Ntot   number of points to use
// Returns: C(k)
//--------------------------------------------------------------
double cpAutoCorr( const std::vector<double>& x, const int k, const int start=0,
		   int Ntot=-1 );


//------------------------------------------------------------------
// cpBinomErr: binomial error calculations
// Inputs:
//   Na    number of accepted events
//   N     total number of events
//   eff   = Na / N
// Returns: sigma(eff) = sqrt( eff (1-eff) / N )
//------------------------------------------------------------------
inline double cpBinomErr( const double eff, const int N ){
  return( sqrt( eff * (1.0-eff) / N ) ); }

inline double cpBinomErr( const int Na, const int N ){
  return( cpBinomErr( (double)Na/(double)N, N ) ); }


//------------------------------------------------------------------
// cpLinReg: linear least squares fit
// Inputs:
//   x     vector of N x measurements
//   y     vector of N y measurements
//   sig   vector of N errors on y
// Outputs
//   a,b   best fit model parameters: f(x) = a + bx
//   V     covariance matrix for a,b (must already be created as 2x2)
//   chi2  chi^2 of fit
// Returns: Number of degrees of freedom in fit
//          -1 for failure
//------------------------------------------------------------------
int cpLinReg( const std::vector<double>& x, const std::vector<double>& y,
	      const std::vector<double>& sig,
	      double& a, double& b, cpMatrix& V, double& chi2 );

#endif
