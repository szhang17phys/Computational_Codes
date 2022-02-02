#ifndef CPODESTEP_H
#define CPODESTEP_H

//************************************************************************
// cpODEstep.hpp : performs steps for numerical solution of 
//                 Ordinary Differential Equations using a variety of methods
//
// Classes
//   cpODEstep : base class (all other inherit from this)
//   cpODEstepEuler: implements Euler algorithm
//   cpODEstepRK2:   implements 2nd Order Runge-Kutta algorithm
//   cpODEstepRK4:   implements 4th Order Runge-Kutta algorithm
//
// Methods
//   cpODEstep( initStep,ODEmodel,tol )
//     initStep = initial step size
//     ODEmodel = pointer to an already-created cpODEmodelXXX class
//     tol      = convergence condition
//
// Modifications
// 20-Jan-10  changes to cpODEstepDP5
//            o add function tolerance() includes both absolute
//              and relative error tolerances
//            o errorCalc() uses tolerance() rather than returning
//              a relative error
//            o modify needNewStep() to include errorCalc() 
//              also change inputs
//            o remove errorCalc() from step()
// 21-Dec-09  add cpODEstepDP5, getStep()
//            add _tol to cpODEstep and all other algo's
// 18-Dec-09  add cpODEstepRK4
//            add t explicitly to step() functions
//            modify update()
// 17-Dec-09  add setStep()
//************************************************************************

#include <iostream>
#include <vector>

#include "cpODEmodel.hpp"

const double DEFTOL = 1.0e-8;  // default tolerance
const int MAXITER   = 10;      // max number of iterations in DP5 stepper

//========================================================================
// cpODEstep: base class
//========================================================================
class cpODEstep {
  // output method
  friend std::ostream& operator <<( std::ostream&, const cpODEstep& );

public:
  // Constructor
  cpODEstep( const double initStep, cpODEmodel* ODEmodel, 
	     const double tol = DEFTOL ) : 
    _h(initStep), _ode(ODEmodel), _tol(tol) {}

  // the stepping function: returns the step size at end of step
  virtual double step() { return _h; }

  // set the step size to a specific value
  void setStep( const double h ){ _h = h; }
  double getStep(){ return _h; }

  // actually change the time, state, and rate vectors in ODEmodel
  void update( const std::vector<double>& s );

protected:
  double _h;            // the current step size
  cpODEmodel* _ode;     // the ODE being evaluated
  double _tol;          // convergence condition (DP5 only)
}; // cpODEstep

//========================================================================
// cpODEstepEuler: Euler step method
//========================================================================
class cpODEstepEuler : public cpODEstep {
public:
  // Constructor
  cpODEstepEuler( const double initStep, cpODEmodel* ODEmodel,
		  const double tol = DEFTOL ) : 
    cpODEstep(initStep,ODEmodel,tol) {}

  double step();              // performs one step of size h

}; // cpODEstepEuler

//========================================================================
// cpODEstepRK2: 2nd order Runge-Kutta step method
//========================================================================
class cpODEstepRK2 : public cpODEstep {
public:
  // Constructor
  cpODEstepRK2( const double initStep, cpODEmodel* ODEmodel,
		const double tol = DEFTOL ) : 
    cpODEstep(initStep,ODEmodel,tol) {}

  double step();             // performs one step of size h

}; // cpODEstepRK2

//========================================================================
// cpODEstepRK4: 4th order Runge-Kutta step method
//========================================================================
class cpODEstepRK4 : public cpODEstep {
public:
  // Constructor
  cpODEstepRK4( const double initStep, cpODEmodel* ODEmodel,
		const double tol = DEFTOL ) : 
    cpODEstep(initStep,ODEmodel,tol) {}

  double step();             // performs one step of size h

}; // cpODEstepRK4

//========================================================================
// cpODEstepDP5: Dormand-Prince adaptive step scheme
//   use 4th vs 5th order Runge-Kutte methods
//   see NR 17.2 and CSM DormandPrince45.java
//========================================================================
class cpODEstepDP5 : public cpODEstep {
public:
  // Constructor
  cpODEstepDP5( const double initStep, cpODEmodel* ODEmodel, 
		const double tol = DEFTOL ) : 
    cpODEstep(initStep,ODEmodel,tol) {}

  double step();             // performs one step of size h

  // trial step
  void trialStep( const double h_try,
		  const std::vector<double>& sN, const std::vector<double>& dsdtN,
		  std::vector<double>& sN1, std::vector<double>& err );
  // is a new step needed
  bool needNewStep( const double h,
		    const std::vector<double>& err, 
		    const std::vector<double>& sN,
		    const std::vector<double>& sN1,
		    double& h_next );

  // calculates the total error from individual component errors
  double errorCalc( const std::vector<double>& err, const std::vector<double>& sN,
		    const std::vector<double>& sN1 );
  // tolerance calculation for equation i
  double tolerance( const double sNi, const double sN1i );

}; // cpODEstepDP5

#endif
