#include <iostream>
#include <iomanip>
#include <cmath>
#include "cpODEmodel.hpp"
using namespace std;

//===========================================================
// cpODEmodel (base class)
//===========================================================

//---------------------------------------
// Constructor
//---------------------------------------
cpODEmodel::cpODEmodel( const vector<double>& s0, const vector<double>& p,
			const double t0 )
{
  _nEqn = (int)s0.size();
  _t0   = t0;
  _t    = t0;
  _s0   = s0;
  _s    = s0;
  _sAnalytic = s0;
  _dsdt.assign( s0.size(),0 );  // rate vector set in specific constructors
} // cpODEmodel

//---------------------------------------
// reset: sets everything to initial conditions
//---------------------------------------
void cpODEmodel::reset()
{
  vector<double> ds0dt( _s0.size(),0 );
  setState( _s0 );
  setTime( _t0 );
  setRate( ds0dt );
} // reset

//---------------------------------------
// Deviation: fractional deviation b/w analytic and numeric state
//   dev = ( Ana - Num ) / Ana             if Ana not too close to zero
//       = ( Ana - Num ) / max( Num,Ana )  otherwise
//---------------------------------------
double cpODEmodel::Deviation( const double ana, const double num ) const {
  if( num == 0.0 && ana == 0.0 ) return 0;
  return( ( abs(ana) > 1.0e-10 ) ? 1.0 - num/ana :           // |ana| >~ 0
  	  ( ana - num ) / ( max( abs(num),abs(ana) ) ) );    // otherwise
} // Deviation

//---------------------------------------
// << operator
//---------------------------------------
ostream& operator <<( ostream& os, const cpODEmodel& cls ){
  int prec = 6;
  string indent( "    " );
  os << "t: " << scientific << setprecision(prec) << cls._t;  // time
  os << " N: " << cls._nEqn << endl;                          // N(eqn's)

  for( int i = 0; i < cls._nEqn; ++i ){                       // state
    os << indent << " s[" << i << "]: " << showpos
       << cls._sAnalytic[i] << " " << cls._s[i]
       << " " << cls.Deviation( cls._sAnalytic[i],cls._s[i] ) << endl;
    os << noshowpos;
  } // i = nEqn

  os << indent << " E:    " << showpos                       // energy
     << cls.Energy(1) << " " << cls.Energy(0)
     << " " << cls.Deviation( cls.Energy(1),cls.Energy(0) ) << endl;
  os << noshowpos;

  return os;
}; // <<


//===========================================================
// cpODEmodelFall1D
//===========================================================

//------------------------------------
// cpODEmodelFall1D::getRate
// Implements free fall near earth surface
//   dx/dt = v
//   dv/dt = -g
// note: ds/dt does not depend on t here
//------------------------------------
vector<double> cpODEmodelFall1D::getRate( const double t, 
					  const vector<double>& s )
{ 
  vector<double> dsdt( _nEqn );
  dsdt[0] = s[1];
  dsdt[1] = -_g;
  return( dsdt );
} // cpODEmodelFall1D::getRate

//------------------------------------
// cpODEmodelFall1D::analyticState
// Solves 1D falling analytically at time t
//   x(t) = x_0 + v_0 [t-t0] - 1/2 g [t-t0]^2
//   v(t) = v_0 - g [t-t0]
//------------------------------------
void cpODEmodelFall1D::setAnalyticState( const double t ) {
  _sAnalytic[0] = _s0[0] + _s0[1]*(t-_t0) - 0.5*_g*(t-_t0)*(t-_t0);
  _sAnalytic[1] = _s0[1] - _g*(t-_t0);
} // cpODEmodelFall1D::setAnalyticState

//------------------------------------
// cpODEmodelFall1D::Energy
// Returns E = K + U for current (ODE-numeric(type=0) or analytic(type=1)) state
//   K = 1/2 v^2
//   U = g x
//------------------------------------
double cpODEmodelFall1D::Energy( const int type ) const {
  double K,U;
  if( type == 0 ){           // ODE numerical solution
    K = 0.5 * _s[1]*_s[1];
    U = _g * _s[0];
  } // type = 0 
  else{                      // analytic solution
    K = 0.5 * _sAnalytic[1]*_sAnalytic[1];
    U = _g * _sAnalytic[0];
  } // type != 0
  return( K+U );
} // cpODEmodelFall1D::Energy


//===========================================================
// cpODEmodelFallAir
//===========================================================

//------------------------------------
// cpODEmodelFallAir::getRate
// Implements 1D falling with air resistance
//   dx/dt = v
//   dv/dt = -g( 1 - v^2 / v_t^2 )
// note: ds/dt does not depend on t here
//------------------------------------
vector<double> cpODEmodelFallAir::getRate( const double t, 
					  const vector<double>& s )
{ 
  vector<double> dsdt( _nEqn );
  dsdt[0] = s[1];
  dsdt[1] = -_g * ( 1.0 - (s[1]*s[1])/(_vt*_vt) );
  return( dsdt );
} // cpODEmodelFallAir::getRate

//------------------------------------
// cpODEmodelFallAir::analyticState
// Returns analytic solution for no air resistance case
//   x(t) = x_0 + v_0 [t-t0] - 1/2 g [t-t0]^2
//   v(t) = v_0 - g [t-t0]
//------------------------------------
void cpODEmodelFallAir::setAnalyticState( const double t ) {
  _sAnalytic[0] = _s0[0] + _s0[1]*(t-_t0) - 0.5*_g*(t-_t0)*(t-_t0);
  _sAnalytic[1] = _s0[1] - _g*(t-_t0);
} // cpODEmodelFallAir::setAnalyticState


//===========================================================
// cpODEmodelSHO
//===========================================================

//------------------------------------
// cpODEmodelSHO::getRate
// Implements Simple Harmonic Oscillator
//   dx/dt = v
//   dv/dt = -omega^2 x
// note: ds/dt does not depend on t here
//------------------------------------
vector<double> cpODEmodelSHO::getRate( const double t,
				       const vector<double>& s )
{
  vector<double> dsdt( _nEqn );
  dsdt[0] = s[1];
  dsdt[1] = -_omega * _omega * s[0];
  return( dsdt );
} // cpODEmodelSHO::getRate

//------------------------------------
// cpODEmodelSHO::analyticState
// Solves SHO analytically at time t (not general - assumes x0=A, v0=0)
//   x(t) = x_0 cos( omega t )
//   v(t) = v_0 - g [t-t0]
//------------------------------------
void cpODEmodelSHO::setAnalyticState( const double t ) {
  _sAnalytic[0] = _s0[0] * cos( _omega*(t-_t0) );
  _sAnalytic[1] = -_s0[0] * _omega * sin( _omega*(t-_t0) );
} // cpODEmodelSHO::setAnalyticState

//------------------------------------
// cpODEmodelSHO::Energy (assumes unit mass)
// Returns E = K + U for current (ODE-numeric(type=0) or analytic(type=1)) state
//   K = 1/2 v^2
//   U = 1/2 k x^2
//------------------------------------
double cpODEmodelSHO::Energy( const int type ) const {
  double K,U;
  double k = _omega*_omega;
  if( type == 0 ){           // ODE numerical solution
    K = 0.5 * _s[1]*_s[1];
    U = 0.5 * k * _s[0]*_s[0];
  } // type = 0 
  else{                      // analytic solution
    K = 0.5 * _sAnalytic[1]*_sAnalytic[1];
    U = 0.5 * k * _sAnalytic[0]*_sAnalytic[0];
  } // type != 0
  return( K+U );
} // cpODEmodelSHO::Energy


//===========================================================
// cpODEmodelDecay
//===========================================================

//------------------------------------
// cpODEmodelDecay::getRate
// Implements exponential decay
//   dN/dt = -lambda N
// note: ds/dt does not depend on t here
//------------------------------------
vector<double> cpODEmodelDecay::getRate( const double t, 
					 const vector<double>& s )
{
  vector<double> dsdt( _nEqn );
  dsdt[0] = -_lambda * s[0];
  return( dsdt );
} // cpODEmodelDecay::getRate

//------------------------------------
// cpODEmodelDecay::analyticState
//   N(t) = N0 exp( -lambda t )
//------------------------------------
void cpODEmodelDecay::setAnalyticState( const double t ) {
  _sAnalytic[0] = _s0[0] * exp( -_lambda * (t-_t0) );
} // cpODEmodelDecay::setAnalyticState
