#include "cpStatistics.hpp"
using namespace std;

//-----------------------------------------------------------------
// cpAveVar
//-----------------------------------------------------------------
double cpAveVar( const vector<double>& data, double& var, int start, int N )
{
  // number of data points to average
  // default is to ave from start to end of vector
  if( N == 0 ) N = data.size() - start;

  // calculate average: Sum_i=start,start+N-1 x_i
  double ave = 0.0;
  for( int i = start; i < start+N; ++i ) ave += data[i];
  ave = ave / N;

  // variance using two-pass formula NR 14.1.8
  // var = (1/(N-1)) [ Sum_i=0,N-1 (x_i - <x>)^2 
  //                   - 1/N(Sum_i=0,N-1 (x_i - <x>))^2) ]
  double s, ep = 0.0;
  var = 0.0;
  for( int i = start; i < start+N; ++i ){
    s    = data[i] - ave;
    ep  += s;               // accumulate  x_i - <x>
    var += s*s;             // accumulate (x_i - <x>)^2
  } // i loop
  var = ( var - ep*ep/N ) / (N-1); // two-pass formula

  return ave;
} // cpAveVar

//-----------------------------------------------------------------
// cpWgtAve: Uncorrelated measurements version
//   <X>     = [Sum X_i / sigma_i^2] err^2
//   1/err^2 = Sum 1/sigma_i^2
//-----------------------------------------------------------------
double cpWgtAve( const vector<double>& data, const vector<double>& sigma,
		 double& err, double& chi2, int& Ndof )
{
  int N = data.size(); // assume that both data & sigma have same dimension

  // # of Degrees of Freedom = N_meas - 1
  // one d.o.f. is used to create average from data points themselves
  Ndof  = N - 1;

  // Accumulate sums for <X> and sigma
  double Xave = 0;
  err  = 0;
  for( int i = 0; i < N; ++i ){
    Xave += data[i] / (sigma[i]*sigma[i]);       // X_i / sigma_i^2
    err  += 1.0 / ( sigma[i]*sigma[i] );         // 1 / sigma_i^2
  } // i loop
  Xave = Xave / err;     // average gets weigthed by Sum (1/sigma_i^2)
  err  = sqrt( 1/err );  // error

  // Now calculate the chi^2 between data and average
  chi2 = 0;
  for( int i = 0; i < N; ++i )
    chi2 += (data[i] - Xave)*(data[i] - Xave) / (sigma[i]*sigma[i]);

  return Xave;
} // cpWgtAve

//-----------------------------------------------------------------
// cpChi2
//-----------------------------------------------------------------
// Pearson's chi^2 version
double cpChi2( const vector<double>& meas, const vector<double>& pred,
	       int& Ndof, const int Nconstr )
{
  // number of bins and initial value of N degrees of freedom
  int N = meas.size();
  Ndof  = N - Nconstr;

  // loop over bins, accumulate sums
  double chi2 = 0;
  for( int i = 0; i < N; ++i ){

    // zero prediction and data means one less degree of freedom
    if( pred[i] == 0.0 && meas[i] == 0.0 ) Ndof--;

    // zero prediction and non-zero data is inconsistent
    else if( pred[i] == 0.0 && meas[i] != 0.0 ) return cpHUGE;

    // non-zero data and prediction, add to sum
    else chi2 += (meas[i] - pred[i])*(meas[i] - pred[i]) / pred[i];
  } // i loop

  return chi2;
} // cpChi2 (Pearson)

// Least Squares version
double cpChi2( const vector<double>& meas, const vector<double>& pred,
	       const vector<double>& sig, int& Ndof, const int Nconstr )
{
  // number of bins and initial value of N degrees of freedom
  int N = meas.size();
  Ndof  = N - Nconstr;

  // loop over bins, accumulate sums
  double chi2 = 0;
  for( int i = 0; i < N; ++i ){

    // don't count bins with sigma = 0 (this is an approximation)
    if( sig[i] == 0.0 ) Ndof--;

    // non-zero data and prediction, add to sum
    else chi2 += (meas[i] - pred[i])*(meas[i] - pred[i]) / (sig[i]*sig[i]);
  } // i loop

  return chi2;
} // cpChi2 (Least Squares)


//-----------------------------------------------------------------
// cpAutoCorr
//-----------------------------------------------------------------
double cpAutoCorr( const vector<double>& x, const int k, const int start,
		   int Ntot )
{
  // calculate ending index for auto-corr sums
  // ind_end = last data index + 1
  int ind_end = x.size() - k;    // use values from start to full size of vector
  int Nuse    = ind_end - start; // N - k
  if( Ntot > 0 ){
    ind_end = start + Ntot - k;  // use only part of full vector
    Nuse    = Ntot - k;
  } // Ntot > 0

  if( Nuse < 2 )
    throw( "cpAutoCorr: must use >= 2 data points" );

  // Loop over index range requested, accumulate sums
  int j;
  double sumij = 0.0;
  double sumi = 0.0, sumi2 = 0.0;
  double sumj = 0.0, sumj2 = 0.0;
  for( int i = start; i < ind_end; ++i ){
    j = i + k;
    sumij += x[i] * x[j];                   // x_{i} x_{i+k}
    sumi  += x[i];  sumi2 += x[i] * x[i];   // x_{i}   & x_{i}^2
    sumj  += x[j];  sumj2 += x[j] * x[j];   // x_{i+k} & x_{i+k}^2
  } // i loop
  sumi2 -= (sumi * sumi / Nuse);            // Var[x_{i}]
  sumj2 -= (sumj * sumj / Nuse);            // Var[x_{i}]

  // construct C(k)
  if( sumi2*sumj2 > 0.0 )
    return( sumij - sumi*sumj/Nuse ) / sqrt(sumi2*sumj2); 
  else
    return cpHUGE;
} // cpAutoCorr

//-----------------------------------------------------------------
// cpLinReg
//-----------------------------------------------------------------
int cpLinReg( const vector<double>& x, const vector<double>& y,
	      const vector<double>& sig,
	      double& a, double& b, cpMatrix& V, double& chi2 )
{
  int N = x.size();   // number of measurements
  int Ndof;
  a = b = 0;

  //-------------------
  // Accumulate sums
  //-------------------
  double w,t;
  double S = 0, Sx = 0, Sy = 0, Stt = 0;

  // first S, S_x, S_y
  for( int i = 0; i < N; ++i ){
    w   = 1.0 / (sig[i]*sig[i]);
    S  += w;
    Sx += x[i] * w;
    Sy += y[i] * w;
  } // i loop

  // now we can accumulate S_{tt} and b sums
  double SxoS = Sx / S;
  for( int i = 0; i < N; ++i ){
    t    = ( x[i] - SxoS ) / sig[i];
    Stt += t * t;
    b   += t * y[i] / sig[i];
  } // i loop

  //------------------
  // Best fit parameters and Covariance
  //------------------
  b = b / Stt;                            // parameters
  a = ( Sy - Sx*b ) / S;
  V[0][0] = ( 1.0 + Sx*Sx/(S*Stt) ) / S;  // covariance matrix
  V[1][1] = 1.0 / Stt;
  V[0][1] = V[1][0] = -Sx / (S*Stt);

  //-------------------
  // Calculate chi^2
  //-------------------
  vector<double> f;                        // prediction at each x
  for( int i = 0; i < N; ++i ) f.push_back( a + b*x[i] );
  chi2 = cpChi2( y,f,sig,Ndof,2 );         // chi^2 w/ 2 constraints

  return Ndof;
} // cpLinReg
