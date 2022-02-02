#include <cmath>

#include "cpInterpolate.hpp"
#include "cpUtils.hpp"
using namespace std;

//======================================================================
// PolyInterp: polynomial interpolation using Neville's algorithm
//======================================================================
double PolyInterp( const vector<double>& xv, const vector<double>& yv,
		   const int ist, const int n, const double x, double& dy )
{
  double y;
  if( (ist+n) > (int)xv.size() || (ist+n) > (int)yv.size() )
    throw( "PolyInterp vector out of bounds" );

  //----------------------
  // The parts of vectors xv,yv that we will use
  // o initialize using iterators on interval [ ist, (ist+n) )
  //----------------------
  vector<double> xa( xv.begin()+ist, xv.begin()+(ist+n) );
  vector<double> ya( yv.begin()+ist, yv.begin()+(ist+n) );

  //----------------------
  // The tableau for n=3
  // y_0=P_0
  //          P_{01}
  // y_1=P_1          P_{012}
  //          P_{12}
  // y_2=P_2
  //----------------------
  vector<double> c( n,0 );  // C_{m,i} = P_{i...(i+m)} - P_{i...(i+m-1)}
  vector<double> d( n,0 );  // D_{m,i} = P_{i...(i+m)} - P_{(i+1)...(i+m)}

  //-----------------------
  // Find closest xa[ns] value to x
  //-----------------------
  int ns = 0;
  double dift;
  double dif = abs( x - xa[0] );
  for( int i = 0; i < n; ++i ){
    if( (dift = abs( x - xa[i] )) < dif ){
      ns  = i;
      dif = dift;
    } // dift < dif
    c[i] = ya[i];      // O(0) values of c's and d's are jsut y_i
    d[i] = ya[i];
  } // i loop

  //-----------------------
  // initial approx of y: 
  // o closest to true P(x) so that corrections are small
  // o decrement ns in preparation for navigation through tableau
  //-----------------------
  y = ya[ns--];

  //--------------------------
  // Construct Neville's Algorithm Tableau
  // N.R. sect. 3.2
  //--------------------------
  double ho,hp,w,den;
  for( int m = 1; m < n; ++m ){      // loop over columns of tableau
    for( int i = 0; i < n-m; ++i ){  // loop over elements in each (new) column
      ho = xa[i]   - x;              // new position differences: column m
      hp = xa[i+m] - x;
      w  = c[i+1]  - d[i];           // previous: C_{(m-1),(i+1)} - D_{(m-1),i}
      if( (den = ho - hp ) == 0.0 ) throw( "PolyInterp error" );
      den  = w / den;
      d[i] = hp * den;               // new: D_{m,i}
      c[i] = ho * den;               // new: C_{m,i}
    } // i loop

    //-------------------------
    // update y with the deviation from this column
    // o starting value is ya[ns] (at x[ns] closest to x)
    // o ns = current position in column
    // o dy = e.g. P_{12} - P_{1}
    // o choose C_{m,j} or D_{m,j} to stay away from tableau boundaries
    //-------------------------
    y += ( dy = ( 2*(ns+1) < (n-m) ? c[ns+1] : d[ns--] ) );
  } // m loop

  return y;
} // PolyInterp
