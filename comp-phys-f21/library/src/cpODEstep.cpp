#include <iostream>
#include <iomanip>
#include <cmath>
#include "cpODEstep.hpp"
using namespace std;

//======================================================================
// cpODEstep::update: move from t_n,s_n --> t_(n+1),s_(n+1)
//   t_n --> t_(n+1) = t_n + h
//   s_n --> s_(n+1)
// Inputs
//   s = s_(n+1)
//======================================================================
void cpODEstep::update( const vector<double>& s )
{
  double t_n1 =  _ode->getTime() + _h;           // time t_(n+1)
  _ode->setTime( t_n1 );
  _ode->setState( s );                           // state s_(n+1)
  _ode->setRate( _ode->getRate(t_n1,s) );        // rate ds/dt_(n+1)
  _ode->setAnalyticState( t_n1 );                // analytic solution at t_(n+1)
} // cpODEstep::update

//---------------------------------------
// << operator
//---------------------------------------
ostream& operator <<( ostream& os, const cpODEstep& cls ){
  int prec = 6;
  os << "h: " << scientific << setprecision(prec) << cls._h;  // step size
  os << " " << *(cls._ode);                                   // state
  return os;
}; // <<

//=======================================================================
// cpODEstepEuler
//   t_[n+1] = t_n + h
//   s_(n+1) = s_n  +  h (ds/dt)_n
//=======================================================================
double cpODEstepEuler::step()
{
  int N    = _ode->getN();                     // number of equations
  double t = _ode->getTime();                  // current time (t_n)
  vector<double> s( _ode->getState() );        // current state (s_n)
  vector<double> dsdt( _ode->getRate(t,s) );   // & rate (ds/dt_n)

  for( int i = 0; i < N; ++i )                 // update to new state
    s[i] += _h * dsdt[i];                      // using current rates

  update( s );                                 // change ode state
  return _h;
} // cpODEstepEuler::step

//=======================================================================
// cpODEstepRK2
//   t_[n+1] = t_n + h
//   k_1 = h ds/dt( t_n, s_n )
//   k_2 = h ds/dt( t_n + h/2, s_n + k_1/2 )
//   s^(i)( t_[n+1] ) = s^(i)( t_n )  +  k_2
//=======================================================================
double cpODEstepRK2::step()
{
  int N      = _ode->getN();                   // number of equations
  double t   = _ode->getTime();                // current time  t_n
  double th2 = t + _h/2.;                      // half-step     t_(n+1/2)
  vector<double> s( _ode->getState() );        // current state (s_n)
  vector<double> k1( _ode->getRate(t,s) );     // & rate (ds/dt_n)

  // Step to t+1/2 using Euler
  vector<double> s12( N );
  for( int i = 0; i < N; ++i ) s12[i] = s[i] + _h * k1[i]/2.;  // s_(n+1/2)
  vector<double> k2( _ode->getRate(th2,s12) );                 // ds/dt_(n+1/2)

  // Now the step to t_(n+1)
  for( int i = 0; i < N; ++i ) s[i] += _h * k2[i];             // s_(n+1)

  update( s );  // change ode state
  return _h;
} // cpODEstepRK2::step

//=======================================================================
// cpODEstepRK4
//   t_[n+1] = t_n + h
//   k_1 = h ds/dt( t_n, s_n )
//   k_2 = h ds/dt( t_n + h/2, s_n + k_1/2 )
//   k_3 = h ds/dt( t_n + h/2, s_n + k_2/2 )
//   k_4 = h ds/dt( t_n + h, s_n + k_3 )
//   s^(i)( t_[n+1] ) = s^(i)( t_n )  +  (1/6)( k_2 + 2k_2 + 2k_3 + k_4 )
//=======================================================================
double cpODEstepRK4::step()
{
  int N      = _ode->getN();              // number of equations
  double t   = _ode->getTime();           // current time  t_n
  double th2 = t + _h/2.;                 // half-step     t_(n+1/2)
  double th  = t + _h;                    // full-step     t_(n+1)
  double h6  = _h / 6.0;

  vector<double> s( _ode->getState() );  // current state (s_n)
  vector<double> sTmp( N );              // temp state for calculations

  // k_1 = ds/dt( t_n,s_n )
  vector<double> k1( _ode->getRate(t,s) );

  // k_2 = ds/dt( t_n + h/2, s_n + k_1/2 )
  for( int i = 0; i < N; ++i ) sTmp[i] = s[i] + _h * k1[i]/2.;
  vector<double> k2( _ode->getRate(th2,sTmp) );

  // k_3 = h ds/dt( t_n + h/2, s_n + k_2/2 )
  for( int i = 0; i < N; ++i ) sTmp[i] = s[i] + _h * k2[i]/2.;
  vector<double> k3( _ode->getRate(th2,sTmp) );

  // k_4 = h ds/dt( t_n + h, s_n + k_3/2 )
  for( int i = 0; i < N; ++i ) sTmp[i] = s[i] + _h * k3[i];
  vector<double> k4( _ode->getRate(th,sTmp) );

  // Now the step to t_(n+1), s_(n+1)
  for( int i = 0; i < N; ++i ) 
    s[i] += h6 * ( k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i] );

  update( s );  // change ode state
  return _h;
} // cpODEstepRK4::step

//=======================================================================
// cpODEstepDP5: Dormand-Price adaptive step algorithm
//   see NR 17.2 for details
//=======================================================================
double cpODEstepDP5::step()
{
  int N    = _ode->getN();                     // number of equations
  double t = _ode->getTime();                  // current time  t_n
  vector<double> s( _ode->getState() );        // current state (s_n)
  vector<double> dsdt( _ode->getRate(t,s) );   // current rate (ds/dt_n)
  vector<double> sTmp( N );                    // temp state for calculations
  vector<double> err( N );                     // 5th-4th order error est

  double h_try = _h;     // trial step size
  double h_next;         // step size estimate for next step

  //--------------------
  // Iterate to find optimal step size for this step and next
  //--------------------
  int trial;
  for( trial = 0; trial < MAXITER; ++trial ){
    trialStep( h_try, s, dsdt, sTmp, err );    // trial step

    if( needNewStep( h_try, err, s, sTmp, h_next ) ) // did not converge
      h_try = h_next;                          //   set h for next iteration
    else                                       // converged for this step size
      break;                                   //   jump out of loop
  } // trial loop

  // algorithm did not converge
  if( trial == MAXITER ){
    cout << "ERROR: DP5 algorithm did not converge" << endl;
    return -1;
  } // iter == MAXITER
  
  //--------------------
  // Update to new state at used step size
  //--------------------
  _h = h_try;      // set step to what was used for this step
  update( sTmp );  // change ode state
  _h = h_next;     // now set step size for next step
  return _h;
} // cpODEstepDP5::step

//======================================
// cpODEstepDP5::trialStep
//   takes a trial step with step size h_try
// Inputs
//   h_try     = step size to try
//   sN, dsdtN = state and derivative at step n
// Outputs
//   sN1       = state at step n+1 (using h_try)
//   err       = vector Delta  y(5)_{n+1} - y(4)_{n+1}
//======================================
void cpODEstepDP5::trialStep( const double h_try, 
			      const vector<double>& sN, const vector<double>& dsdtN,
			      vector<double>& sN1, vector<double>& err )
{
  //---------------------
  // Dormand-Prince 5(4) coefficients from table in NR 17.2
  //---------------------
  static const double c[7] = { 0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0 };
  static const double a[7][7] = {
    { 0.0 },
    { 0.2 },
    { 3.0/40.0,        9.0/40.0 },
    { 44.0/45.0,      -56.0/15.0,      32.0/9.0 },
    { 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0 },
    { 9017.0/3168.0,  -355.0/33.0,     46732.0/5247.0,  49.0/176.0,  -5103.0/18656.0 },
    { 35.0/384.0,      0.0,            500.0/1113.0,    125.0/192.0, -2187.0/6784.0,     11.0/84.0 }
  };
  static const double e[7] = {  // b_i^(5th) - b_i^(4th)
    71.0/57600.0, 0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0
  };

  int N    = _ode->getN();                     // number of equations
  double t = _ode->getTime();                  // current time  t_n
  vector<vector<double> > k;                   // rate vector at each internal step

  //----------------------
  // Construct the sum over the 6 internal steps
  //   k_i = ds/dt( t_n + c_i, s_n + Sum_{j<i} h a_[i,j] k_j )
  //   s_{n+1} = s_{n} + Sum_j a_[7,j] k_j
  // at the end of these loops sN1 contains the s_{n+1}
  //----------------------
  k.push_back( dsdtN );                    // k_1 = ds/dt_n
  for( int istp = 1; istp < 7; ++istp ){   // loop over internal steps

    for( int i = 0; i < N; ++i ){          // loop over equations
      // accumulate terms for this internal step
      sN1[i] = sN[i];  // start at s_n
      // the state vector for internal step istp
      // when istp = 6 we have s_{n+1} from sum of 6 k_istp's
      for( int jstp = 0; jstp < istp; ++jstp )
	sN1[i] += h_try * a[istp][jstp] * k[jstp][i];
    } // i = equations in state vector

    // next rate estimate
    // when istp = 6 we have the ds/dt_{n+1}
    // which we could (but won't) use to start the next step
    k.push_back( _ode->getRate( t + c[istp]*h_try, sN1 ) );
  } // istp = internal steps

  //----------------------
  // Calculate errors: 5th - 4th order truncation
  // note: we use 7 terms here (incl. the final ds/dt_{n+1})
  //----------------------
  for( int i = 0; i < N; ++i ){
    err[i]    = 0.0;
    for( int j = 0; j < 7; ++j )
      err[i] += e[j] * k[j][i];
    err[i] *= h_try;
  } // i = equations

} // cpODEstepDP5::trialStep

//======================================
// cpODEstepDP5::needNewStep
// returns TRUE if a new step is needed based on error
//   calculates size of next step (no PI control)
//     h_{n+1} = safety * h_{n} * ( 1/err )^alpha
// Inputs
//   h      = current step size
//   err    = vector of Delta = s(5)_{n+1} - s(4)_{n+1}
//   sN,sN1 = s_{n}, s_{n+1}
// Outputs
//   h_next = new step size
//======================================
bool cpODEstepDP5::needNewStep( const double h,
				const vector<double>& err, 
				const vector<double>& sN,
				const vector<double>& sN1,
				double& h_next )
{
  static const double alpha=0.2, safe=0.9;
  static const double minscale=0.2, maxscale=10.0;
  double scale;

  // find the total error normalized to the tolerance
  double del = errorCalc( err, sN, sN1 );

  //---------------
  // Successful step
  //   increase step (w/in limits)
  //---------------
  if( del <= 1.0 ){
    scale  = (del == 0.0) ? maxscale : safe*pow(del,-alpha);
    h_next = h * min( scale,maxscale );
    return false;  // new step not needed
  } // del <= 1

  //---------------
  // Failed step: error too large
  //   decrease step (w/in limits)
  //---------------
  else{
    scale  = safe * pow( del,-alpha ); 
    h_next = h * max( scale,minscale );
    return true;  // new step needed
  } // del > 1

} // cpODEstepDP5::needNewStep

//======================================
// cpODEstepDP5::errorCalc
// returns error [0.0-1.0] for this step 
//   as Euclidean sum of individual state 
//   equation error terms
//   normalized to the tolerance
//   del = 1/Neqn Sum_i=1,N ( Delta_i / tolerance_i )^2
// Inputs
//   err    = vector of Delta = s(5)_{n+1} - s(4)_{n+1}
//   sN,sN1 = s_{n}, s_{n+1}
//======================================
double cpODEstepDP5::errorCalc( const vector<double>& err, 
				const vector<double>& sN,
				const vector<double>& sN1 )
{
  double esum = 0.0;
  double scale;
  int N = _ode->getN();

  // calculate error relative to tolerance
  for( int i = 0; i < N; ++i ){
    scale = tolerance( sN[i],sN1[i] );
    esum += (err[i]*err[i]) / (scale*scale);
  } // i = equations

  return( sqrt( esum/N ) );
} // cpODEstepDP5::errorCalc

//======================================
// cpODEstepDP5::tolerance
// returns the error tolerance scale
//   scale = tol + max(|s_n|,|s_n+1|) tol
//   includes both relative tolerance (on state, s) and
//     absolute tolerance (mainly for oscillatory cases)
//   uses max(s_n,s_n+1) to avoid division by zero
// Inputs
//   sNi, sN1i: the state at steps n, n+1 for the ith
//              ODE equation
//======================================
double cpODEstepDP5::tolerance( const double sNi, const double sN1i ) 
{
  return( _tol + _tol * max( abs(sNi),abs(sN1i) ) );
} // cpODEstepDP5::tolerance
