//*****************************************************************
// Run tests of computer numerics
//
// Test 0: numeric_limits printout
//     $ ./testNumerics 0
// Test 1: representation of floating point numbers
//     $ ./testNumerics 1 <number> [double]
// Test 2: accuracy/roundoff tests
//     $ ./testNumerics 2
//
// Modifications
// 29-Jun-18  created from a merger of DoubleRep and Accuracy
//*****************************************************************
using namespace std;

#include <iostream>
#include <iomanip>
#include <bitset>
#include <limits>
#include <cstring>

#include "cpUtils.hpp"
#include "cpDoubleRep.hpp"

//----------------------------------------------
// Function Prototypes
//----------------------------------------------
int DoubleRep( int argc, char* argv[] );
int Accuracy( int argc, char* argv[] );
pair<double,double> QuadSolve( const double a, const double b, const double c,
			       int opt=0 );
double Dev( const double x, const double y );


//==============================================================
int main( int argc, char* argv[] ){
  if( argc < 2 ){
    cout << "Usage: testNumerics <test> [other params]" << endl;
    cout << "    test = 0  ==> numeric_limits printout" << endl;
    cout << "    test = 1  ==> float/double representation" << endl;
    cout << "    test = 2  ==> accuracy/roundoff" << endl;
    return -1;
  } // argc < 2

  //------------------------------------------------
  // Determine test to perform from first argument
  //------------------------------------------------
  int test = atoT<int>( argv[1] );

  int stat = 0;
  switch( test ){
  // Test 1: representation of float or double
  case 1:
    stat = DoubleRep( argc, argv );
    break;

    // Test 2: accuracy and roundoff error
  case 2:
    stat = Accuracy( argc, argv );
    break;

  // Test 0: numeric_limits
  default:
    cout << "Max: unsigned int           = " 
	 << numeric_limits<unsigned int>::max() << endl;
    cout << "Max: unsigned long long int = " 
	 << numeric_limits<unsigned long long int>::max() << endl;
    cout << "Max: double                 = " 
	 << numeric_limits<double>::max() << endl;
    cout << "Min: double                 = " 
	 << numeric_limits<double>::min() << endl;
    cout << "Epsilon: double             = " 
	 << scientific << setprecision(20)
	 << numeric_limits<double>::epsilon() << endl;
    cout << "1 + Epsilon: double         = " 
	 << 1.0 + numeric_limits<double>::epsilon() << endl;
    cout << "1 + 1e-16: double           = " 
	 << 1.0 + 1.0e-16 << endl;
    cout << "Max: float                  = " 
	 << numeric_limits<float>::max() << endl;
    cout << "Min: float                  = " 
	 << numeric_limits<float>::min() << endl;
    cout << "Epsilon: float              = " 
	 << numeric_limits<float>::epsilon() << endl;
    break;
  } // switch test

  return 0;
} // main


//================================================================
// DoubleRep: computer representation of float/double numbers
//=================================================================

// FloatToUint and DoubleToUint copy the bits from a floating point
// number to the appropriately sized unsigned integer without
// invoking undefined behavior
uint32_t FloatToUint(float x) {
  uint32_t result;
  std::memcpy(&result, &x, sizeof(uint32_t));
  return result;
}

uint64_t DoubleToUint(double x) {
  uint64_t result;
  std::memcpy(&result, &x, sizeof(uint64_t));
  return result;
}

int DoubleRep( int argc, char* argv[] )
{
  if( argc < 3 ){
    cout << "Usage: testNumerics 1 <number> [add anothe argument for double precision]"<< endl;
    return -1;
  } // argc < 3

  // Single Precision input (default)
  if( argc == 3 ){
    float d = atoT<float>( argv[2] );            // read input from command line
    Ullong ull = (Ullong) FloatToUint(d);                // convert to 64-bit for PrintRep
    bitset<MAXBITS> bits( BitSetULL(ull) );  // bit representation
    cout << d << " = " << PrintRep( F_EXPBITS,F_MANBITS,bits ) << endl;
  } // argc == 2 (float requested)

  // Double Precision requested
  else{
    double d = atoT<double>( argv[2] );           // read input from command line
    bitset<MAXBITS> bits( BitSetULL( DoubleToUint(d) ) );  // bit representation
    cout << d << " = " << PrintRep( D_EXPBITS,D_MANBITS,bits ) << endl;
  } // argc > 2 (double requested)

  return 0;
} // DoubleRep


//=================================================================
// Accuracy: tests of numerical accuracy and roundoff
//=================================================================
int Accuracy( int argc, char* argv[] ){

  double a=1.0, b=1.0, c;
  pair<double,double> x0,x1,x2;
  cout << "a: " << a << "  b: " << b << endl;

  // find solutions to equation using different methods
  for( int i = -1; i >= -20; --i ){
    c  = pow( 10.0,i );         // this takes lots of CPU time !
    x0 = QuadSolve( a,b,c,0 );  // numerically "exact" solution
    x1 = QuadSolve( a,b,c,1 );  // "normal" solution
    x2 = QuadSolve( a,b,c,2 );  // "inverted" solution

    // write out results
    cout << showpos << scientific << setprecision(5) 
	 << "c: " << c 
	 << "  x1: " << Dev(x1.first,x0.first)   
	 << " " << Dev(x2.first,x0.first) 
	 << "\t x2: " << Dev(x1.second,x0.second) 
	 << " " << Dev(x2.second,x0.second)
	 << endl; 
  } // i loop

  return 0;
} // Accuracy

//=============================================================
// QuadSolve: returns solutions to quadratic equation
//   ax^2 + bx + c = 0 ==> x = [-b +- sqrt(b^2 - 4ac)] / 2a  opt=1
//                           = -2c / [b +- sqrt(b^2 - 4ac)]  opt=2
//                           = q/a,c/q                       opt=0 (default)
//                             ( q = -1/2 [b + sgn(b)sqrt(b^2 - 4ac)] )
// The opt=0 solution is the most numerically stable
//   see NR sect 5.6
//=============================================================
pair<double,double> QuadSolve( const double a, const double b, const double c,
			       int opt )
{
  pair<double,double> x;
  double q;
  double sgnb = ( b >= 0 ) ? 1.0 : -1.0;
  double discrim = b*b - 4.0*a*c;
  discrim = ( discrim > 0 ) ? sqrt( discrim ) : 0;

  switch( opt ){

  // normal solution
  case 1:
    x.first  = ( -b + discrim ) / ( 2.0*a );
    x.second = ( -b - discrim ) / ( 2.0*a );
    break;

  // alternate solution
  case 2:
    x.first  = ( -2.0*c ) / ( b + discrim );
    x.second = ( -2.0*c ) / ( b - discrim );
    break;

  // numerically stable solution
  default:
    q = -0.5 * ( b + sgnb*discrim );
    x.first  = c / q;
    x.second = q / a;
    break;

  } // switch

  return x;
} // QuadSolve

//==============================================================
// Dev: calculates (x-y)/y with various safeguards
//==============================================================
double Dev( const double x, const double y )
{
  if( y == 0 )               return( -numeric_limits<double>::max() );
  if( isinf(x) || isnan(x) ) return(  numeric_limits<double>::max() );
  return( x/y - 1.0 );
} // Dev
