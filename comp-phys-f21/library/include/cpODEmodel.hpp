#ifndef CPODEMODEL_H
#define CPODEMODEL_H

//************************************************************************
// cpODEmodel.hpp : Ordinary Differential Equations for different models
//   Problem is broken into a set of coupled 1st order ODEs
//     ds^(i)/dt = f( t,s^(0),...,s^(N-1) )
//   These classes specify the state vector (s) and the functions ds/dt
//
// Classes
//   cpODEmodel : base class (all other inherit from this)
//   cpODEmodelFall1D  : falling motion in 1D
//   cpODEmodelFallAir : falling motion in 1D with air resistance
//   cpODEmodelSHO     : Simple Harmonic Oscillator
//   cpODEmodelDecay   : exponential decay
//
// Methods
//   cpODEmodel(s0,p,[t0] )
//     s0 = vector of initial values at t0
//     p  = vector with other params for this model
//     t0 = initial time (default=0)
//
// Modifications
// 22-Sep-10  modify << operator (uses multiple lines)
// 13-Jan-10  modify << operator
//            add deviation function
//            remove getAnalyticState, change getState to return Num or Ana
// 21-Dec-09  add cpODEmodelFallAir
// 18-Dec-09  add t to getRate
// 17-Dec-09  add reset() and getAnalyticState
// 15-Dec-09  add cpODEmodelDecay
//************************************************************************

#include <iostream>
#include <vector>

#include "cpUtils.hpp"

//========================================================================
// cpODEmodel: base class
//========================================================================
class cpODEmodel {
  // output method
  friend std::ostream& operator <<( std::ostream&, const cpODEmodel& );

public:
  // Constructor: sets initial state
  // need to instantiate specific models to set initial rate
  cpODEmodel( const std::vector<double>& s0, const std::vector<double>& p,
	      const double t0=0.0 );

  // the function that implements the ds/dt calculations (dummy for base class)
  // most of the time this does not depend explicitly on t
  // but for some ODEs it does
  virtual std::vector<double> getRate(const double t, const std::vector<double>& s ){ 
    return _s; }

  // the analytic solution (if it exists, if not use virtual method)
  virtual void setAnalyticState( const double t ) { 
    _sAnalytic.assign( _nEqn,0 ); }

  // energy (unit mass): type = 0(ODE calc), 1(analytic)
  virtual double Energy( const int type=0 ) const { return 0; }

  // fractional deviation b/w analytic and numeric state
  double Deviation( const double ana, const double num ) const;

  // set the time and the State & Rate vectors
  void setTime( const double t ) { _t = t; }
  void setState( const std::vector<double>& s ) { _s = s; }
  void setRate( const std::vector<double>& dsdt ) { _dsdt = dsdt; }
  void reset();

  // access to protected elements (this is probably bad form)
  int getN() { return _nEqn; }
  double getTime() { return _t; }
  std::vector<double> getState( const int type=0 ){ 
    return( (type == 0 ) ? _s : _sAnalytic ); }

protected:
  int _nEqn;              // number of equations
  double _t0;             // initial time value
  double _t;              // current time value
  std::vector<double> _s0;     // initial state vector: e.g. x(t0),v(t0),a(t0)...
  std::vector<double> _s;      // current state vector: e.g. x(t),v(t),a(t)...
  std::vector<double> _sAnalytic; // current state vector for analytic solution
  std::vector<double> _dsdt;   // current rate vector: ds^(i)/dt
}; // cpODEmodel

//========================================================================
// cpODEmodelFall1D: 1D falling point mass
//   s    = [x,v]
//   s0   = [y0,v0]
//   p[0] = g (accel due to gravity)
//========================================================================
class cpODEmodelFall1D : public cpODEmodel {
public:
  // Constructor: sets initial state and rate
  cpODEmodelFall1D( const std::vector<double>& s0,  const std::vector<double>& p,
		    const double t0=0.0 ) :
    cpODEmodel(s0,p,t0), _g(p[0]) {
    setRate( getRate(t0,s0) );
  }

  std::vector<double> getRate( const double t, const std::vector<double>& s );
  void setAnalyticState( const double t );
  double Energy( const int type=0 ) const;

protected:
  double _g;    // accel due to gravity

}; // cpODEmodelFall1D

//========================================================================
// cpODEmodelFallAir: 1D falling object w/ air resistant
//   s    = [x,v]
//   s0   = [y0,v0]
//   p[0] = g (accel due to gravity)
//   p[1] = v_t (terminal velocity)
// note: no energy calculated
//       analytic state = state w/ *no* air resistance
//========================================================================
class cpODEmodelFallAir : public cpODEmodel {
public:
  // Constructor: sets initial state and rate
  cpODEmodelFallAir( const std::vector<double>& s0,  const std::vector<double>& p,
		    const double t0=0.0 ) :
    cpODEmodel(s0,p,t0), _g(p[0]), _vt(p[1]) {
    setRate( getRate(t0,s0) );
  }

  std::vector<double> getRate( const double t, const std::vector<double>& s );
  void setAnalyticState( const double t );

protected:
  double _g;    // accel due to gravity
  double _vt;   // terminal velocity

}; // cpODEmodelFallAir

//========================================================================
// cpODEmodelSHO: Simple Harmonic Oscillator
//   s    = [x,v]
//   s0   = [Amplitude,0]  (assumed to be released at maximum)
//   p[0] = T (period of oscillation)
//========================================================================
class cpODEmodelSHO : public cpODEmodel {
public:
  // Constructor: sets initial state and rate
  cpODEmodelSHO( const std::vector<double>& s0,  const std::vector<double>& p,
		    const double t0=0.0 ) :
    cpODEmodel(s0,p,t0), _omega(2*cpPI/p[0]) {
    setRate( getRate(t0,s0) );
  }

  std::vector<double> getRate( const double t, const std::vector<double>& s );
  void setAnalyticState( const double t );
  double Energy( const int type=0 ) const;

protected:
  double _omega;    // angular frequency = sqrt(k/m)

}; // cpODEmodelSHO

//========================================================================
// cpODEmodelDecay: Exponential Decay - dN/dt = -lambda N
//   s    = [N]
//   s0   = [N0]
//   p[0] = lambda (1/tau)
// note: no energy calculated
//========================================================================
class cpODEmodelDecay : public cpODEmodel {
public:
  // Constructor: sets initial state and rate
  cpODEmodelDecay( const std::vector<double>& s0,  const std::vector<double>& p,
		    const double t0=0.0 ) :
    cpODEmodel(s0,p,t0), _lambda(p[0]) {
    setRate( getRate(t0,s0) );
  }

  std::vector<double> getRate( const double t, const std::vector<double>& s );
  void setAnalyticState( const double t );

protected:
  double _lambda;    // decay constant

}; // cpODEmodelDecay

#endif
