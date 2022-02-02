#include <cmath>
#include <string>

#include "cpPDE.hpp"
using namespace std;

//====================================================================
// Analytic Solutions
//====================================================================
// GaussWave ---------------------------------------------------------
double GaussWave( const double x, const double t, const vector<double>& p )
{
  double v   = p[0];
  double sig = p[1];
  double nrm = sqrt( 2.*cpPI ) * sig;
  double arg = ( x - v*t ) / sig;
  return( exp( -0.5*arg*arg ) / nrm );
} // GaussWave

// GaussDiffuse -------------------------------------------------------
double GaussDiffuse( const double x, const double t, const vector<double>& p )
{
  double D    = p[0];
  double sig0 = p[1];
  double sig2 = 2.0*D*t + sig0*sig0;    // sig^2 = 2Dt + sig_0^2
  double nrm  = sqrt( 2.*cpPI * sig2 );
  double arg2 = x*x / sig2;
  return( exp( -0.5*arg2 ) / nrm );
} // GaussDiffuse


//===================================================================
// Boundary Conditions
//===================================================================
// Bnd1dDirichlet: returns u_0 or u_j -------------------------------
double Bnd1dDirichlet( const int j, const vector<double>& uN, const double val )
{
  return( val );
} // Bnd1dDirichlet

// Bnd1dCyclic: returns u_0 or u_J ----------------------------------
double Bnd1dCyclic( const int j, const vector<double>& uN, const double val )
{
  int J = uN.size() - 1;             // index of upper bound
  if( j == 0 )      return uN[J];    // u_0^n+1 = u_J^n
  else if( j == J ) return uN[J-1];  // u_J^n+1 = u_J-1^n
  else              return 0.0;      // bad index
} // Bnd1dCyclic

// Bnd1dNeumann: returns u_-1 or u_J+1 -------------------------------
// note: must use val = 2 Delta(x) du/dx_0,J when calling
double Bnd1dNeumann( const int j, const vector<double>& uN, const double val )
{
  int J = uN.size() - 1;                     // index of upper bound
  if( j == 0 ) return( uN[1] - val );        // u_{-1}  = u_1     - 2 dx du_0/dx
  else if( j == J ) return( uN[J-1] + val ); // u_{J+1} = u_{J-1} + 2 dx du_J/dx
  else              return 0.0;              // bad index
} // Bnd1dNeumann



//====================================================================
// << operator for cpPDEiv
//=====================================================================
ostream& operator <<( ostream& os, const cpPDEiv& cls ){
  int prec = 4;
  os << "Time: " << cls._t << endl;

  for( int j = 0; j < cls._dimU; ++j )
    cout << "  " << "x: " << showpos << scientific << setprecision(prec)
	 << cls.xPos(j) << "  Numeric: " << cls._u[j] << endl;

  return os;
}; // <<

//====================================================================
// cpPDEivFluxFTCS: FTCS solution of PDE
//====================================================================
void cpPDEivFluxFTCS::update()
{
  //------------------
  // state at t = n+1 - for j = 1,...,N-1
  //------------------
  vector<double> u_np1( _nStepX,0 );
  for( int j = 1; j < _nStepX-1; ++j )
    u_np1[j] = _u[j] - _alpha * ( _u[j+1] - _u[j-1] );

  //------------------
  // update state at boundaries
  //------------------
  if( _bndType == BNDneu ){            // Neumann Boundary Condition
    double uM1  = _BndFunc( 0,         _u, 2.0*_dx*_bndVal[0] ); // u_{-1}
    double uJp1 = _BndFunc( _nStepX-1, _u, 2.0*_dx*_bndVal[1] ); // u_{J+1}
    u_np1[0]         = _u[0]         - _alpha * ( _u[1] - uM1 );
    u_np1[_nStepX-1] = _u[_nStepX-1] - _alpha * ( uJp1 - _u[_nStepX-2] );
  }
  else{                                // Dirichlet or Cyclic
    u_np1[0]       = _BndFunc( 0,         _u, _bndVal[0] );
    u_np1[_nStepX] = _BndFunc( _nStepX-1, _u, _bndVal[1] );
  }

  // update time & state
  _t += _dt;
  _u  = u_np1;
} // cpPDEivFluxFTCS::update

//====================================================================
// cpPDEivFluxLax: LAX solution of Flux PDE
//====================================================================
void cpPDEivFluxLax::update()
{
  //------------------
  // state at t = n+1 - for j = 1,...,N-1
  //------------------
  vector<double> u_np1( _nStepX,0 );
  for( int j = 1; j < _nStepX-1; ++j )
    u_np1[j] = 0.5 * ( _u[j+1] + _u[j-1] ) - _alpha * ( _u[j+1] - _u[j-1] );

  //------------------
  // update state at boundaries
  //------------------
  if( _bndType == BNDneu ){            // Neumann Boundary Condition
    double uM1  = _BndFunc( 0,         _u, 2.0*_dx*_bndVal[0] ); // u_{-1}
    double uJp1 = _BndFunc( _nStepX-1, _u, 2.0*_dx*_bndVal[1] ); // u_{J+1}
    u_np1[0] = 0.5 * ( _u[1] + uM1 ) 
      - _alpha * ( _u[1] - uM1 );
    u_np1[_nStepX-1] = 0.5 * ( uJp1 + _u[_nStepX-2] ) 
      - _alpha * ( uJp1 - _u[_nStepX-2] );
  }
  else{                                // Dirichlet or Cyclic
    u_np1[0]         = _BndFunc( 0,         _u, _bndVal[0] );
    u_np1[_nStepX-1] = _BndFunc( _nStepX-1, _u, _bndVal[1] );
  }

  // update time & state
  _t += _dt;
  _u  = u_np1;
} // cpPDEivFluxLax::update

//====================================================================
// cpPDEivDiffuseFTCS: FTCS solution of Diffusion PDE
//====================================================================
void cpPDEivDiffuseFTCS::update()
{
  //------------------
  // state at t = n+1 - for j = 1,...,N-1
  //------------------
  vector<double> u_np1( _nStepX,0 );
  for( int j = 1; j < _nStepX-1; ++j )
    u_np1[j] = _u[j] + _alpha * ( _u[j+1] - 2.*_u[j] + _u[j-1] );

  //------------------
  // update state at boundaries
  //------------------
  if( _bndType == BNDneu ){            // Neumann Boundary Condition
    double uM1  = _BndFunc( 0,         _u, 2.0*_dx*_bndVal[0] ); // u_{-1}
    double uJp1 = _BndFunc( _nStepX-1, _u, 2.0*_dx*_bndVal[1] ); // u_{J+1}
    u_np1[0]       = _u[0]       
      + _alpha * ( _u[1] - 2.*_u[0] + uM1 );
    u_np1[_nStepX-1] = _u[_nStepX-1] 
      + _alpha * ( uJp1 - 2.*_u[_nStepX-1] + _u[_nStepX-2] );
  }
  else{                                // Dirichlet or Cyclic
    u_np1[0]         = _BndFunc( 0,       _u, _bndVal[0] );
    u_np1[_nStepX-1] = _BndFunc( _nStepX-1, _u, _bndVal[1] );
  }

  // update time & state
  _t += _dt;
  _u  = u_np1;
} // cpPDEivDiffuseFTCS::update

//====================================================================
// cpPDEivDiffuseRW: Random Walk solution of Diffusion PDE
//====================================================================
void cpPDEivDiffuseRW::update()
{
  //------------------
  // Update walker positions: 
  // move each walker by +-dx with equal probability
  // note: this allows walkers to fall out of the region
  //   be careful when translating to a distribution
  //------------------
  for( int i = 0; i < _nW; ++i ){
    if( _Rand->randDouble() < 0.5 ) _xW[i] -= _dx;
    else                            _xW[i] += _dx;
  } // i = walkers

  // update time & state
  _t += _dt;
  _u = Walker2Func( _xW );
} // cpPDEivDiffuseRW::update

double cpPDEivDiffuseRW::Integrate( const vector<double>& u )
{
  double sum = 0;
  for( int j = 0; j < _nStepX; ++j ) sum += u[j];
  return( sum );
} // cpPDEivDiffuseRW::Integrate

//---------------------------------------------------------------------
vector<double> cpPDEivDiffuseRW::Func2Walker( const vector<double>& u )
{
  vector<double> xW;
  double scl = _nW / _norm;     // scales _u array to _nW entries

  //--------------------
  // Fill walker vector based on values of _u[j] scaled to number of walkers
  // note: we will end this loop with a vector that does not (generally) 
  //   have _nW entries
  //--------------------
  for( int j = 0; j < _nStepX; ++j ){
    //--------------
    // loop over number of walkers in each bin
    //   calculated from scaled value of u_j^n (rounded to nearest int)
    //   note: this will generally *not* give _nW entries in the vector
    //    b/c of rounding issues
    // load x-position for each walker into vector
    //--------------
    for( int i = 0; i < cpRound(scl*_u[j]); ++i )
      xW.push_back( _xLo + j*_dx );  // use x_j for each walker
  } // j = x-positions

  _nW = (int)xW.size();  // N(walkers) updated to number of entries in vector
  return xW;
} // cpPDEivDiffuseRW::Func2Walker

//------------------------------------------------------------------------
vector<double> cpPDEivDiffuseRW::Walker2Func( const vector<double>& w )
{
  vector<double> u( _nStepX,0 );

  // Loop over walkers: fill u like a histogram with each walker's position
  int j;
  for( int i = 0; i < _nW; ++i )
    if( (j = FindIndex(_xW[i])) >= 0 ) u[j] += 1;

  // normalize u to same sum as initial distribution
  double scl = _norm / Integrate( u );
  for( j = 0; j < _nStepX; ++j ) u[j] *= scl;

  return u;
} // cpPDEivDiffuseRW::Walker2Func

//-----------------------------------------------------------------------
int cpPDEivDiffuseRW::FindIndex( const double x )
{
  if( x < _xLo ) return -1;                 // x below lower boundary

  for( int j = 0; j < _nStepX; ++j )
    if( x < (_xLo + (j+1)*_dx) ) return j;  // returns index j if x < x_{j+1}

  return -1;                                // x above upper bound
} // cpPDEivDiffuseRW::FindIndex
