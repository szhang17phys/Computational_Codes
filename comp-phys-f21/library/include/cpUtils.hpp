#ifndef CPUTILS_H
#define CPUTILS_H

//************************************************************************
// cpUtils: collection of generally useful functions
// Modifications
// 30-Aug-18  cpTimer returns time in micro-seconds
// 05-Jan-11  add cpRound
// 01-Sep-10  create cpUtils.cpp b/c of problems that gmake/make has creating
//            a library using only a header file
// 09-Apr-10  add cpSwap
// 11-Mar-10  add atoT
// 06-Feb-10  add cpHUGE
// 13-Oct-09  add throw and cpTimer
//************************************************************************

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>

//-----------------------------------
// Constants
//-----------------------------------
const double cpPI      = acos( -1.0 );         // pi
const double cpSqrt2PI = sqrt( 2.0*cpPI );     // sqrt( 2pi )
const double cpHUGE    = 1.0e90;               // big number test

//------------------------
// Error Handling: redefine throw using a preprocessor macro
//------------------------
#define throw( message ){ \
  std::cout << "ERROR: " << message << std::endl; \
  std::cout << "  in FILE: " << __FILE__ << " at LINE: " << __LINE__ <<std::endl; \
  exit(1); }

//------------------------
// Timing Class
// o start timer by instantiating a member
// o print out elapsed time from start using << operator
//------------------------
class cpTimer {
  // printout elapsed time:
  //   use form <string> <double> <string> for easy reading
  friend std::ostream& operator <<( std::ostream& os, const cpTimer& cls ){
    os << "cpTimer(Elapsed-time): " << cls.elapsed() << std::endl; 
    return os; }

public:
  cpTimer() : _t0(clock()) {};         // start clock when instantiated
  void reset(){ _t0 = clock(); }       // reset clock
  double elapsed() const {             // elapsed time [micro-seconds]
    return( 1.0e6 * (double)(clock()-_t0) / (double)CLOCKS_PER_SEC ); }

private:
  clock_t _t0;                // start time
}; // cpTimer

//--------------------------
// Other Functions
// note atoT must always be invoked using e.g.:
//   x = atoT<double>( ch ) if x is a double, or
//   x = atoT<unsigned int>( ch ) if x is an unsigned int...
//   since T is not used in the function arguments
//--------------------------
template <class T>               // general char --> type T transform
inline T atoT( char* ch ){
  T d; std::istringstream iss( (std::string)ch ); iss >> d; return d; }


inline double atod( char* ch ){  // char --> double
  return( atoT<double>(ch) ); }  //  for backward compatibility

//-------------------------
// cpSwap: swaps a <--> b
//-------------------------
template <class T>
inline void cpSwap( T& a, T& b ){ T tmp = a; a = b; b = tmp; }

//-------------------------
// cpRound: rounds a double to nearest integer
//-------------------------
inline int cpRound( double d ){ return( (int)floor(d + 0.5) ); }

#endif
