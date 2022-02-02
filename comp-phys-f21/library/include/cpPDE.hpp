#ifndef CPPDE_H
#define CPPDE_H

//**********************************************************************
// cpPDE.hpp: classes for solving PDEs
// Note - these are not very general. Have to create a daughter class
//   for every problem
//
// Base Method: cpPDEiv - 1+1D Initial Value PDEs
// Flux Conserving Daughter Methods
//   cpPDEivFluxFTCS: FTCS solution, Dirichlet boundary conditions
//   cpPDEivFluxLax:  Lax solution, Dirichlet boundary conditions
// Diffusion Daughter Methods
//   cpPDEivDiffFTCS: FTCS solution
//   cpPDEivDiffRW:   Random Walk solution
//
// Note on time and space grids
//   nStepX/T = total number of points in space/time
//   Spatial grid defined from xLo - xHi
//     x[0] = xLo ; x[nStepX-1] = xHi - dx
//   Time grid defined from tLo - tHi
//     t0 = tLo ; tf = tHi - dt
//
// Modifications
// 06-Sep-19  change use of _nStepX, _nStepT
// 03-Sep-19  add variable _dimU and function xPos to allow for complex state vectors _u
// 05-Jan-11  add _alpha to daughter classes
//            add cpPDEivDiffuseRW
//            make reset() virtual
// 04-Jan-11  create as library
//            cpPDEiv: add bndType,bndVal[2]
//            add boundary conditions functions
// 04-Nov-10  created
//**********************************************************************

#include <iostream>
#include <iomanip>
#include <vector>

#include "cpUtils.hpp"
#include "cpRandom.hpp"

//===========================================================
// Analytic Solutions for Initial Value problems: same type as FuncXT
// GaussWave: Gaussian w/ mean = v*t; sigma = p[1]; v = p[0]
// GaussDiffuse: Gaussian w/ time-dependent sigma: D = p[0]; sigma_0 = p[1]
//===========================================================
typedef double (*FuncXT)( const double, const double, const std::vector<double>& );
double GaussWave( const double x, const double t, const std::vector<double>& p );
double GaussDiffuse( const double x, const double t, const std::vector<double>& p );

//===========================================================
// Boundary Conditions: same type as BoundX
// o j    : =0 (lower bound); =J (upper bound)
// o uN   : solution vector at t=n
// o val  : value associated with boundary condition
//          = u_0 or u_J for Dirichlet
//          = du_0/dx or du_J/dx for Neumann
//
// Return u_0 or u_J
//   Bnd1dDirichlet: Dirichlet boundary conditions
//   Bnd1dCyclic:    simple Cyclic boundary condition
// Return u_{-1} or u_{J+1}
//   Bnd1dNeumann:   Neumann boundary conditions (2nd order approx)
//===========================================================
const int BNDdir = 0;  // Boundary condition types
const int BNDcyc = 1;
const int BNDneu = 2;
typedef double (*BoundX)( const int, const std::vector<double>&, const double );
double Bnd1dDirichlet( const int j, const std::vector<double>& uN, const double val );
double Bnd1dCyclic( const int j, const std::vector<double>& uN, const double val );
double Bnd1dNeumann( const int j, const std::vector<double>& uN, const double val );

//===========================================================
// cpPDEiv class: base class for Initial Value PDEs
//   p    = parameters of PDE
//   u0   = vector (size nStepX+1) of initial values u_j^0
//   bndType = boundary condition type
//   bndVal  = array of boundary values/derivs
//   xLo,xHi = boundaries in x (last x-point = xHi-dx)
//   nStepX  = number of steps in x
//   tLo,tHi = boundaries in t (last t-point = tHi-dt)
//   nStepT  = number of steps in t
//===========================================================
class cpPDEiv {
  // output method
  friend std::ostream& operator <<( std::ostream&, const cpPDEiv& );

public:
  // Constructor
  cpPDEiv( const std::vector<double>& p, const std::vector<double>& u0,
	   const int bndType, double bndVal[2],
	   const double xLo=0., const double xHi=1., const int nStepX=50,
	   const double tLo=0., const double tHi=1., const int nStepT=10 ) :
    _bndType(bndType), _bndVal(bndVal),
    _xLo(xLo),_xHi(xHi),_nStepX(nStepX), _tLo(tLo),_tHi(tHi),_nStepT(nStepT), 
    _u0(u0), _u(u0)
  { _dx = (_xHi -_xLo)/_nStepX; _dt = (_tHi - _tLo)/nStepT; _dimU = nStepX;
    switch( bndType ){
    case BNDdir: _BndFunc = Bnd1dDirichlet; break;
    case BNDcyc: _BndFunc = Bnd1dCyclic;    break;
    case BNDneu: _BndFunc = Bnd1dNeumann;   break;
    }
  }

  virtual void update() { _t += _dt; }  // make time step - impose boundaries
  virtual void reset(){ _u = _u0; _t = _tLo; }  // reset to initial state
  double xPos( int j ) const { return( (j < _nStepX) ? _xLo + j*_dx : _xLo + (j-_nStepX)*_dx ); }

protected:
  int _bndType;             // boundary type
  double* _bndVal;          // lower [0] and upper [1] boundary values or derivs
  BoundX _BndFunc;          // function for boundary calculation

  double _xLo,_xHi;         // spatial region
  int _nStepX;
  double _dx;

  double _tLo,_tHi;         // time region
  int _nStepT;
  double _dt;

  double _t;                // current time
  int _dimU;                // dimension of state vector (allows for complex _u)
  std::vector<double> _u0;       // u_j^0 (initial state)
  std::vector<double> _u;        // u_j^n (current state)
}; // cpPDEiv


//==============================================================
// cpPDEivFluxFTCS: FTCS solution of flux-conserving PDE w/ F = v u
//   du/dt = -v du/dx
//     p[0] = v
//   analytic solutions: u = f(x - vt)
//==============================================================
class cpPDEivFluxFTCS : public cpPDEiv {
public:
  // Constructor
  cpPDEivFluxFTCS( const std::vector<double>& p, const std::vector<double>& u0,
	       const int bndType, double bndVal[2],
	       const double xLo=0., const double xHi=1., const int nStepX=50,
	       const double tLo=0., const double tHi=1., const int nStepT=10 ) :
    cpPDEiv( p,u0,bndType,bndVal,xLo,xHi,nStepX,tLo,tHi,nStepT ), _v(p[0])
  { _alpha = 0.5 * _v * _dt / _dx; }

  void update();

protected:
  double _v;
  double _alpha;
}; // cpPDEivFluxFTCS


//==============================================================
// cpPDEivFluxLax: Lax solution of flux-conserving PDE w/ F = v u
//   du/dt = -v du/dx
//     p[0] = v
//   analytic solutions: u = f(x - vt)
//==============================================================
class cpPDEivFluxLax : public cpPDEiv {
public:
  // Constructor
  cpPDEivFluxLax( const std::vector<double>& p, const std::vector<double>& u0,
   	      const int bndType, double bndVal[2],
	      const double xLo=0., const double xHi=1., const int nStepX=50,
	      const double tLo=0., const double tHi=1., const int nStepT=10 ) :
    cpPDEiv( p,u0,bndType,bndVal,xLo,xHi,nStepX,tLo,tHi,nStepT ), _v(p[0]) 
  { _alpha = 0.5 * _v * _dt / _dx; }

  void update();

protected:
  double _v;
  double _alpha;
}; // cpPDEivFluxLax 


//==============================================================
// cpPDEivDiffuseFTCS: FTCS solution of 1D diffusion equation
//   du/dt = D d^2u/dx^2
//    p[0] = D
//==============================================================
class cpPDEivDiffuseFTCS : public cpPDEiv {
public:
  // Constructor
  cpPDEivDiffuseFTCS( const std::vector<double>& p, const std::vector<double>& u0,
   	       const int bndType, double bndVal[2],
	       const double xLo=0., const double xHi=1., const int nStepX=50,
	       const double tLo=0., const double tHi=1., const int nStepT=10 ) :
    cpPDEiv( p,u0,bndType,bndVal,xLo,xHi,nStepX,tLo,tHi,nStepT ), _D(p[0]) 
  { _alpha = _D * _dt / (_dx*_dx); }

  void update();

protected:
  double _D;
  double _alpha;
}; // cpPDEivDiffuseFTCS


//==============================================================
// cpPDEivDiffuseRW: Random Walk solution of 1D diffusion equation
//   du/dt = D d^2u/dx^2
//    p[0] = D
//    p[1] = N(walkers) - initial
//    p[2] = seed for random number generator
// Note: make sure that grid is set up s.t.
//    dx = sqrt( 2 D dt )
//==============================================================
class cpPDEivDiffuseRW : public cpPDEiv {
public:
  // Constructor
  cpPDEivDiffuseRW( const std::vector<double>& p, const std::vector<double>& u0,
   	       const int bndType, double bndVal[2],
	       const double xLo=0., const double xHi=1., const int nStepX=50,
	       const double tLo=0., const double tHi=1., const int nStepT=10 ) :
    cpPDEiv( p,u0,bndType,bndVal,xLo,xHi,nStepX,tLo,tHi,nStepT ),
    _D(p[0]), _nW0((int)p[1]), _seed((Ullong)p[2])
  { if( _dx != sqrt( 2.*_D*_dt ) ) throw("cpPDEivDiffuseRW: use dx^2 = 2 D dt");
    _norm = Integrate( _u0 ); _nW = _nW0; _xW0 = Func2Walker( _u0 ); _xW = _xW0; 
    _Rand = new cpRandom( _seed ); }

  void update();
  // reset recaluclates _nW since it may have changed during evolution
  void reset(){ _u = _u0; _t = _tLo, _nW = (int)_xW0.size(); _xW = _xW0; }

  double Integrate( const std::vector<double>& u );        // Sum u_j
  int FindIndex( const double x );                    // index of x in vector

  std::vector<double> Func2Walker( const std::vector<double>& u ); // walker x's from u_j
  std::vector<double> Walker2Func( const std::vector<double>& w ); // u_j from walker x's
  // note: use Func2Walker with care since it can change _nW

protected:
  double _D;             // diffusion coefficient

  int _nW0;              // requested N(walkers)
  int _nW;               // used N(walkers) [may be different than request]
  double _norm;          // integral of initial distribution
  std::vector<double> _xW0;   // initial x-positions of walkers
  std::vector<double> _xW;    // current x-positions of walkers


  cpRandom* _Rand;   // random number generator
  Ullong _seed;      // random seed
}; // cpPDEivDiffuseRW


#endif
