//***********************************************************************
// tests: runs ODE tests
//
// $ ./decay <test> <stepper> <model> [params...]  :  see USAGE below
//
// Modifications
// 22-Sep-10  change output format for easier reading
// 22-Jan-10  add test=4: testFallAir
//            change DefaultParam(), TimeScale(), initModel(): vector of params
//            chages to testXXX() to use vector of params
// 13-Jan-10  modify writeStep
//            move deviation to cpODEmodel
//            add DefaultParam, TimeScale, initModel functions
//            add model parameter: modifies tests 1,2
// 22-Dec-09  add testFallAir, testConvvsNstep, cpTimer
//            change output of all tests to be consistent
// 18-Dec-09  add RK4 stepper
// 15-Dec-09  created
//***********************************************************************
using namespace std;

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

#include "cpUtils.hpp"
#include "cpODEmodel.hpp"
#include "cpODEstep.hpp"
#include "cpInterpolate.hpp"

cpODEstep* getStepper( int stepID, double h, cpODEmodel* ode,
		       double conv = 1.0e-14 );
void writeStep( const int nStep, cpODEstep* stepper );
vector<double> DefaultParam( int model );
double TimeScale( int model, double param );
cpODEmodel* initModel( int model, vector<double>& param );
int testModelXvsT( int stepID, int argc, char* argv[] );
int testModelDvsH( int stepID, int argc, char* argv[] );
int testConvvsNstep( int stepID, int argc, char* argv[] );
int testFallAir( int stepID, int argc, char* argv[] );

//-----------------------------
// Constants
//-----------------------------
const double g   = 9.8;      // acceleration of gravity near Earth [m/s]
const double dum = 0;        // dummy value

// IDs for various stepping algorithms
const int idEuler = 1;       // Euler
const int idRK2   = 2;       // Runge-Kutta 2nd Order
const int idRK4   = 3;       // Runge-Kutta 4th order
const int idDP5   = 4;       // Dormand-Prince stepper

// IDs for various models
const int idDecay   = 1;     // radioactive decay
const int idFall    = 2;     // 1D falling w/out air resist
const int idSHO     = 3;     // simple harmonic oscillator
const int idFallAir = 4;     // 1D falling w/ air resist

//================================================================
int main( int argc, char* argv[] )
{
  //---------------------
  // Read command line inputs: test number
  //---------------------
  if( argc < 4 ){
    cout << "--------------------------------------" << endl;
    cout << "USAGE: decay <test> <stepper> <model> [params]" << endl;
    cout << "stepper    " << idEuler << "=Euler " << idRK2 << "=RK2 "
	 << idRK4 << "=RK4 " << idDP5 << "=DP5" << endl;
    cout << "model      " << idDecay << "=Decay " << idFall << "=Falling "
	 << idSHO << "=SHO " << idFallAir << "=AirResist" << endl;
    cout << endl;
    cout << "Test=1  -  X(t_i) vs t_i" << endl;
    cout << "  param 1  step size in natural units (def=0.1)" << endl;
    cout << "  param 2  t_max in natural units (def=5)" << endl;
    cout << "  param 3  lambda/height/period (def=1/50/1)" << endl;
    cout << endl;
    cout << "Test=2  -  Deviation vs h (at t_max)" << endl;
    cout << "  param 1  h_min in natural units (def=1.0e-3)" << endl;
    cout << "  param 2  h_max in natural units (def=1.0e-2)" << endl;
    cout << "  param 3  t_max in natural units (def=5)" << endl;
    cout << endl;
    cout << "Test=3  -  SHO Convergence vs N(steps)" << endl;
    cout << "  param 1  min N(step) (def=1.0e2)" << endl;
    cout << "  param 2  max N(step) (def=1.0e7)" << endl;
    cout << "  param 3  DP5 convergence (def=1.0e-14)" << endl;
    cout << endl;
    cout << "Test=4  -  Falling w/ Air Resistance" << endl;
    cout << "  param 1  step size in natural units (def=0.01)" << endl;
    cout << "  param 2  terminal velocity (def=21 m/s)" << endl;
    cout << "  param 3  initial height (def=20 m)" << endl;
    cout << "--------------------------------------" << endl;
    return -1;
  } // argc < 4

  int test_no = atoi( argv[1] );   // test ID
  int stepID  = atoi( argv[2] );   // stepper ID

  //---------------------
  // The Tests (each is a separate function)
  //---------------------
  int stat;
  switch( test_no ){

  case 1:                       // Decay: N(t) vs t
    stat = testModelXvsT( stepID, argc,argv );  break;

  case 2:                       // Decay: Deviation vs h
    stat = testModelDvsH( stepID, argc,argv );  break;

  case 3:                       // Convergence vs Nstep
    stat = testConvvsNstep( stepID, argc,argv );     break;

  case 4:                       // Falling w/ air resistance
    stat = testFallAir( stepID, argc,argv );  break;

  default:                      // bad test number
    cout << "Invalid test number: " << test_no << endl;
    break;

  } // switch test_no

  return 0;
} // main

//===============================================================
// getStepper: returns correct stepper based on stepID and ODEmodel
//===============================================================
cpODEstep* getStepper( int stepID, double h, cpODEmodel* ode, double conv )
{
  cpODEstep* stepper;

  switch( stepID ){
  case idEuler:                                     // Euler method
    stepper = new cpODEstepEuler( h,ode );  break;
  case idRK2:                                       // RK2 method
    stepper = new cpODEstepRK2( h,ode );    break;
  case idRK4:                                       // RK4 method
    stepper = new cpODEstepRK4( h,ode );    break;
  case idDP5:                                       // DP5 method
    stepper = new cpODEstepDP5( h,ode,conv ); break;
  default:                                          // invalid stepper
    cout << "getStepper: invalid stepper ID: " << stepID << endl;
    return (cpODEstep *)0;
    break;
  } // switch stepID

  return stepper;
}

//================================================================
// writeStep: write out info about current step in a consistent way
//================================================================
void writeStep( const int nStep, cpODEstep* stepper )
{
  cout.setf( ios_base::left, ios_base::adjustfield );
  cout << setw(8) << nStep << " ";
  cout.setf( ios_base::right, ios_base::adjustfield );
  cout << *stepper;
} // writeStep

//================================================================
// initModel: creates ODE model with initial conditions
// parameters
//   model = model ID
//   param = vector of model params
//================================================================
cpODEmodel* initModel( int model, vector<double>& param )
{
  cpODEmodel* ode;
  vector<double> s0;
  vector<double> pars;

  switch( model ){
  case idDecay:                 // Decay
    s0.push_back( 1000 );       //   N_0 
    pars.push_back( param[0] ); //   lambda
    ode = new cpODEmodelDecay( s0,pars );
    break;

  case idFall:                  // Falling
    s0.push_back( param[0] );   //   height
    s0.push_back( 0 );          //   v_0 = 0
    pars.push_back( g );        //   gravity
    ode = new cpODEmodelFall1D( s0,pars );
    break;

  case idSHO:                   // SHO
    s0.push_back( 1 );          //   x_0 = A
    s0.push_back( 0 );          //   v_0 = 0
    pars.push_back( param[0] ); //   period (T)
    ode = new cpODEmodelSHO( s0,pars );
    break;

  case idFallAir:               // Falling w/ air resist
    s0.push_back( param[0] );   //   height
    s0.push_back( 0 );          //   v_0 = 0
    pars.push_back( g );        //   gravity
    pars.push_back( param[1] ); //   terminal velocity
    ode = new cpODEmodelFallAir( s0,pars );
    break;
  } // switch

  return ode;
} // initModel


//================================================================
// DefaultParam: returns default time-related parameter for a model
// parameters
//   model = model ID
//================================================================
vector<double> DefaultParam( int model )
{
  vector<double> param;

  switch( model ){
  case idDecay:             // Decay: lambda
    param.push_back( 1 );   break;
  case idFall:              // Falling: height
    param.push_back( 50 );  break;
  case idSHO:               // SHO: period
    param.push_back( 1 );   break;
  case idFallAir:           // Falling w/ air resist: height, v_t
    param.push_back( 20 ); param.push_back( 21 );  break;
  } // switch

  return param;
} // DefaultParam

//================================================================
// TimeScale: sets the natural timescale for a model
// parameters
//   model = model ID
//   param = time parameter of model
//================================================================
double TimeScale( int model, double param )
{
  double tscale = 1.0;                          // default

  switch( model ){
  case idDecay:                                 // Decay: 1/lambda
    tscale = 1 / param;  break;

  case idFall:                                  // Falling: t_fall
    tscale = sqrt( 2.0*param/g );  break;       //  (for v_0 = 0)

  case idSHO:                                   // SHO: period
    tscale = param;  break;

  case idFallAir:                               // Air Resist: t_fall(no air)
    tscale = sqrt( 2.0*param/g );  break;       //  (for v_0 = 0)
  } // switch

  return tscale;
} // TimeScale

//================================================================
// testModelXvsT: X(t) vs t
// parameters
//   argv[4] = step size [h]
//   argv[5] = t_max
//   argv[6] = lambda/height/period
//================================================================
int testModelXvsT( int stepID, int argc, char* argv[] )
{
  //---------------------
  // Which model are we using
  //---------------------
  int model = atoi( argv[3] );

  //---------------------
  // Read parameters from command line
  // param = lambda         for Decay
  //       = initial height for Falling
  //       = period         for SHO
  //---------------------
  vector<double> param = DefaultParam( model );  // model parameter
  if( argc > 6 ) param[0] = atod( argv[6] );
  double tscale = TimeScale( model,param[0] );   // tscale for model

  double tmax = 5 * tscale;                      // evolve from 0 - t_max
  if( argc > 5 ) tmax   = atod( argv[5] ) * tscale;

  double h = 0.1 * tscale;                       // step size
  if( argc > 4 ) h      = atod( argv[4] ) * tscale;

  //---------------------
  // Set up the ODE model and the stepper
  //--------------------
  cpODEmodel* ode = initModel( model, param );
  cpODEstep* stepper = getStepper( stepID, h, ode );

  //----------------------
  // Step until t = t_max
  //----------------------
  int nStep  = 0;
  writeStep( nStep,stepper );              // print out initial state
  while( ode->getTime() < tmax ){
    nStep++;                               // update step count
    stepper->step();                       // perform step
    writeStep( nStep,stepper );            // print out new state
  } // while t < tmax

  return 0;
} // testModelXvsT

//================================================================
// testModelDvsH: Exponential Decay - (N_num(t)-N_ana(t))/N_ana(t) vs h
// parameters
//   argv[4] = h_max
//   argv[5] = h_min
//   argv[6] = t_max
//================================================================
int testModelDvsH( int stepID, int argc, char* argv[] )
{
  //---------------------
  // Which model are we using
  //---------------------
  int model = atoi( argv[3] );

  //---------------------
  // Read parameters from command line
  // Set timescale and time at which to evaluate function for deviation
  // and step size limits
  //---------------------
  vector<double> param  = DefaultParam( model ); // set up params for model
  double tscale = TimeScale( model,param[0] );   // timescale for model
  double tmax   = 5 * tscale;                // evolve from 0 - t_max
  if( argc > 6 ) tmax = atod( argv[6] ) * tscale;

  double hmax = 0.01, hmin = 0.001;
  double h = hmin;
  if( argc > 4 ) hmin = atod( argv[4] ) * tscale;
  if( argc > 5 ) hmax = atod( argv[5] ) * tscale;

  //---------------------
  // Set up the ODE model and the stepper
  //--------------------
  cpODEmodel* ode = initModel( model, param );
  cpODEstep* stepper = getStepper( stepID, h, ode );

  //---------------------
  // Loop over steps from hmin to hmax
  // Choose equal spacing in log(h)
  //---------------------
  int nVals      = 50;
  double loghmin = log10( hmin );
  double loghmax = log10( hmax );
  double dlogh   = ( loghmax - loghmin ) / nVals;

  for( int i = 0; i < nVals+1; ++i ){   // inclusive: [hmin,hmax]
    h = pow( 10, (loghmin + i*dlogh) ); // the new step size
    ode->reset();                       // reset state to initial conditions
    stepper->setStep( h );              // set step size for this loop

    //----------------------
    // Step until t = t_max
    // Print out final values
    //----------------------
    int nStep = 0;
    while( ode->getTime() <= tmax ){
      stepper->step();
      nStep++;
    } // while t <= tmax
    writeStep( nStep,stepper ); // write out final results

  } // i loop

  return 0;
} // testModelDvsH


//================================================================
// testConvvsNstep: Convergence vs Nstep
//   uses SHO with T=1, A=1 over 5 cycles
// parameters
//   argv[4] = smallest number of steps to try
//   argv[5] = largest number of steps to try
//   argv[6] = convergence param for DP5 method
//================================================================
int testConvvsNstep( int stepID, int argc, char* argv[] )
{
  //---------------------
  // Which model are we using
  //---------------------
  int model = atoi( argv[3] );

  //---------------------
  // Read parameters from command line
  //---------------------
  double nStepLo = 1.0e2;                      // range for Nsteps
  if( argc > 4 ) nStepLo = atod( argv[4] );
  double nStepHi = 1.0e7;
  if( argc > 5 ) nStepHi = atod( argv[5] );

  double conv = 1.0e-14;                       // convergence
  if( argc > 6 ) conv    = atod( argv[6] );    //   DP5 stepper only

  // number of points to try: distributed evenly over log plot
  int nPts = 25;
  double lognSL = log10( nStepLo );
  double lognSH = log10( nStepHi );
  double dlogN  = ( lognSH - lognSL ) / nPts;

  //---------------------
  // Set up the SHO and the stepper
  //--------------------
  vector<double> param   = DefaultParam( model ); // set up params for model
  double tscale  = TimeScale( model,param[0] );   // timescale for model
  double nscales = 5;                         // no. of time scales to evolve
  if( model == idFall ) nscales = 1;
  double tmax   = nscales * tscale;           // evolve from 0 -> tmax
  double h      = tmax / nStepLo;             // intial step size

  cpODEmodel* ode = initModel( model, param );
  cpODEstep* stepper = getStepper( stepID, h, ode );

  // State vectors at end of steps over cycles
  double tStep = 0;
  int nStep    = 0;
  double algT;
  cpTimer* algTime = new cpTimer(); // Timing container

  //----------------------
  // Adaptive step size stepper is special
  //   figures out how many steps to take to reach target convergence
  //   w/in the stepper itself
  //----------------------
  if( stepID == idDP5 ){
    algTime->reset();                    // start timer
    while( tStep < tmax ){
      //cout << stepper->getStep() << endl;
      stepper->step();
      tStep = ode->getTime();
      nStep++;
    } // while t < nT
    algT = algTime->elapsed();           // time since reset
    writeStep( nStep,stepper );          // write out results of this step trial
    cout << "  clock-time: " << algT << endl;

    return 0;
  } // idDP5

  //----------------------
  // For all other steppers: Loop over step size trials
  //----------------------
  double lognStep = lognSL;
  while( lognStep <= lognSH ){
    ode->reset();                        // reset model to initial conditions

    nStep = (int)pow( 10,lognStep );     // step size for this trial
    h     = tmax / nStep;
    stepper->setStep( h );

    algTime->reset();                    // reset timer
    for( int i = 1; i <= nStep; ++i )    // Step thru cycles
      stepper->step();
    algT = algTime->elapsed();           // time since reset

    writeStep( nStep,stepper );          // write out results of this step trial
    cout << "  clock-time: " << algT << endl;
    lognStep += dlogN;                   // update to next number of steps
  } // while nStep <= nStepHi

  return 0;
} // testConvvsNstep

//================================================================
// testFallAir: 1D falling (w/ air resistance)
// parameters
//   argv[4] = step size [s]
//   argv[5] = terminal velocity [m/s]
//   argv[6] = initial height [m]
//================================================================
int testFallAir( int stepID, int argc, char* argv[] )
{
  //---------------------
  // Which model are we using
  //---------------------
  int model = idFallAir;

  //---------------------
  // Read parameters from command line
  // param = height, v_t
  //---------------------
  vector<double> param = DefaultParam( model );    // model parameters
  double y_0 = param[0];                           // initial height
  if( argc > 6 ) y_0 = atod( argv[6] );
  double tscale = TimeScale( model,param[0] );     // tscale for model

  double v_t = param[1];                           // terminal veloc
  if( argc > 5 ) v_t = atod( argv[5] );

  double h = 0.01 * tscale;                        // step size
  if( argc > 4 ) h = atod( argv[4] ) * tscale;

  //---------------------
  // Set up the ODE model and the stepper
  //--------------------
  cpODEmodel* ode = initModel( model, param );
  cpODEstep* stepper = getStepper( stepID, h, ode );

  //----------------------
  // Step until x <= 0
  //----------------------
  vector<double> tim,pos,vel;                // will contain t,x,v for each step
  vector<double> sNum = ode->getState( 0 );  // current numeric state
  writeStep( 0,stepper );                    // print out initial state
  tim.push_back( ode->getTime() );
  pos.push_back( sNum[0] );
  vel.push_back( sNum[1] );

  int nStep = 1;                       // initial state is 1st step
  while( sNum[0] > 0.0 ){
    stepper->step();                   // perform step
    sNum = ode->getState( 0 );         // get numeric state
    writeStep( nStep,stepper );        // print out new state
    tim.push_back( ode->getTime() );   // save t,x,v
    pos.push_back( sNum[0] );
    vel.push_back( sNum[1] );
    nStep++;                           // update step count
  } // while t < tmax

  //-----------------------
  // Now interpolate to find t at x=0
  // Use this time to find v(x=0)
  //----------------------
  double dT,dV;
  double t_x0 = PolyInterp( pos,tim, nStep-3,3, 0.0, dT );
  double v_x0 = PolyInterp( tim,vel, nStep-3,3, t_x0,dV );
  double vAna = -sqrt( 2.0 * g * y_0 );
  cout << -1 << " Velocity at the Ground:" << endl;
  cout << "v(x=0)       : " << v_x0 << " m/s (numerical)       error = " 
       << dV << endl;
  cout << "v(x=0)       : " << vAna << " m/s (no air resist)" << endl;
  cout << "Air/NoAir -1 : " << v_x0/vAna - 1.0 << endl;

  return 0;
} // testFallAir
