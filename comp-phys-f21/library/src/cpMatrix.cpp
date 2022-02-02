#include <iomanip>

#include "cpMatrix.hpp"
#include "cpUtils.hpp"
using namespace std;

//------------------------------------------------------------------
// cpMatrixInit
//------------------------------------------------------------------
void cpMatrixInit( cpMatrix& A, int nRow, int nCol, double val ){
  vector<double> row( nCol, val );
  for( int i = 0; i < nRow; ++i )
    A.push_back( row );
} // cpmatrixInit

//------------------------------------------------------------------
// cpMatrixPrint
//------------------------------------------------------------------
void cpMatrixPrint( const cpMatrix& A )
{
  ios_base::fmtflags origFormat = cout.flags();  // save current format
  int n = A.size();
  for( int i = 0; i < n; ++i ){
    for( int j = 0; j < n; ++j )
      cout << "  " << showpos << scientific << setprecision(3) << A[i][j];
    cout << endl;
  } // i loop
  cout.flags( origFormat );            // restore old format
} // printMatrix

//------------------------------------------------------------------
// cpMatrixInvGJ: Matrix Inversion by Gauss-Jordan Elimination
//   solves matrix equation: A Y = 1
//------------------------------------------------------------------
int cpMatrixInvGJ( cpMatrix& A )
{
  int n = A.size();    // size returns the number of rows in the matrix

  //-----------------
  // Bookkeeping arrays and variables
  //-----------------
  vector<int> ipiv( n,0 );     // has the pivot for col icol been used
  vector<int> indxr( n,0 );    // row in which pivot element i was orig located
  vector<int> indxc( n,0 );    // col in which pivot element i was orig located
  double big = 0;              // for choosing largest element as pivot
  int icol,irow;               // column,row of current pivot
  double pivinv;               // 1 / pivot element
  double prv;

  //--------------------------------
  // Main Loop over columns (i) to be reduced
  //--------------------------------
  for( int i = 0; i < n; ++i ){

    //----------------------
    // Search for the pivot element
    //----------------------
    big = 0;
    for( int j = 0; j < n; ++j ){       // search for next unused row j
      if( ipiv[j] != 1 ){
	for( int k = 0; k < n; ++k ){   // loop over col k in this row
	  if( ipiv[k] == 0 ){
	    if( abs(A[j][k]) >= big ){  // check for biggest element
	      big  = abs(A[j][k]);
	      irow = j;                 // save row and col
	      icol = k;
	    } // A[j][k] >= big
	  } // ipiv[k] == 0
	} // k loop
      } // ipiv[j] == 0
    } // j loop

    // update pivot use count
    // this can be > 1 because of row swaps
    ++ipiv[icol];

    //-------------------------------------
    // We now have the pivot row and col
    // Interchange rows in A to put the pivot element on the diagonal
    // This implies a column swap in Y
    //--------------------------------------
    if( irow != icol ){
      for( int l = 0; l < n; ++l ) cpSwap( A[irow][l],A[icol][l] );
    } // irow != icol

    // keep track of the original rows and col of this pivot element
    indxr[i] = irow;
    indxc[i] = icol;

    //--------------------
    // Set pivot element in A to 1
    //--------------------
    if( A[icol][icol] == 0.0 ) return( 1 );   // singular matrix
    pivinv        = 1.0 / A[icol][icol];
    A[icol][icol] = 1.0;

    //--------------------
    // Replace the pivot row of A with what will end up 
    //   on the right-hand-side of the equation
    //   A[icol][icol] --> 1 / pivot   (pivot = prev val of A[icol][icol])
    //   A[icol][l]    --> A[icol][l] / pivot
    //---------------------
    for( int l = 0; l < n; ++l ) A[icol][l] *= pivinv;

    //---------------------
    // Reduce the other rows (expect for pivot row)
    //   such that pivot column elements are zero
    // We don't touch the pivot row b/c it already contains the inverse
    // The pivot column elements from the right-hand-side matrix are
    //   then copied into the pivot column of A
    //---------------------
    for( int ll = 0; ll < n; ++ll ){
      if( ll != icol ){       // skip pivot row
	prv = A[ll][icol];    // save prev. value of elem to become zero
	A[ll][icol] = 0.0;    // set col to zero
	for( int l = 0; l < n; ++l )     // adjust all elem's in this row
	  A[ll][l] -= A[icol][l] * prv;
      } // ll != icol
    } // ll loop
  } // i = main loop over columns

  //-------------------------------------------
  // We're now done with the inversion
  // Unscramble the transformed matrix A by interchanging columns
  //   in the reverse order that they were built up
  //-------------------------------------------
  for( int l = n-1; l >= 0; --l ){
    if( indxr[l] != indxc[l] ){
      for( int k = 0; k < n; ++k )
	cpSwap( A[k][indxr[l]], A[k][indxc[l]] );
    } // indxr != indxc
  } // l loop

  return 0;
} // cpMatrixInvGJ

//------------------------------------------------------------------
// cpMatrixInvGJlong: Matrix Inversion by Gauss-Jordan Elimination
//   solves matrix equation: A Y = 1
// Debug version - does not replace A as we go along
//------------------------------------------------------------------
int cpMatrixInvGJlong( cpMatrix& A, const bool debug )
{
  int n = A.size();    // size returns the number of rows in the matrix

  //---------------------
  // Create the right-hand-side matrix that will
  //   initiall = I
  //---------------------
  cpMatrix I( n, vector<double>(n,0) );
  for( int l = 0; l < n; ++l ) I[l][l] = 1;

  string separator("-------------------------------------------------------" );
  if( debug ){
    cout << separator << endl;
    cout << "Original Matrix A" << endl;
    cpMatrixPrint( A );
    cout << "Original Matrix I" << endl;
    cpMatrixPrint( I );
  } // debug

  //-----------------
  // Bookkeeping arrays and variables
  //-----------------
  vector<int> ipiv( n,0 );     // has the pivot for col icol been used
  vector<int> indxr( n,0 );    // row in which pivot element i was orig located
  vector<int> indxc( n,0 );    // col in which pivot element i was orig located
  double big = 0;              // for choosing largest element as pivot
  int icol,irow;               // column,row of current pivot
  double pivinv;               // 1 / pivot element
  double prv;

  //--------------------------------
  // Main Loop over columns (i) to be reduced
  //--------------------------------
  for( int i = 0; i < n; ++i ){

    //----------------------
    // Search for the pivot element
    //----------------------
    big = 0;
    for( int j = 0; j < n; ++j ){       // search for next unused row j
      if( ipiv[j] != 1 ){
	for( int k = 0; k < n; ++k ){   // loop over col k in this row
	  if( ipiv[k] == 0 ){
	    if( abs(A[j][k]) >= big ){  // check for biggest element
	      big  = abs(A[j][k]);
	      irow = j;                 // save row and col
	      icol = k;
	    } // A[j][k] >= big
	  } // ipiv[k] == 0
	} // k loop
      } // ipiv[j] == 0
    } // j loop

    // update pivot use count
    // this can be > 1 because of row swaps
    ++ipiv[icol];

    if( debug ){
      cout << endl << separator << endl;
      cout << "PIVOT: " << i 
	   << " : row = " << irow << " col = " << icol << endl;
      cout << "A = " << endl;  cpMatrixPrint( A );
      cout << "I = " << endl;  cpMatrixPrint( I );
    } // debug

    //-------------------------------------
    // We now have the pivot row and col
    // Interchange rows in A to put the pivot element on the diagonal
    // .....
    //--------------------------------------
    if( irow != icol ){
      for( int l = 0; l < n; ++l ){
	cpSwap( A[irow][l],A[icol][l] );
	cpSwap( I[irow][l],I[icol][l] );
      } // l loop
    } // irow != icol

    if( debug ){
      cout << endl << "SWAP pivot row so that pivot is on the diagonal: " 
	   << endl;
      if( irow != icol ){
	cout << "A = " << endl;  cpMatrixPrint( A );
	cout << "I = " << endl;  cpMatrixPrint( I );
      } // irow != icol
      else cout << "  already on diagonal" << endl;
    } // debug


    // keep track of the original rows and col of this pivot element
    indxr[i] = irow;
    indxc[i] = icol;

    //--------------------
    // Divide the pivot rows of A,I by the pivot
    //---------------------
    if( A[icol][icol] == 0.0 ) return( 1 );   // singular matrix
    pivinv = 1.0 / A[icol][icol];
    for( int l = 0; l < n; ++l ){
      A[icol][l] *= pivinv;
      I[icol][l] *= pivinv;
    } // l loop

    if( debug ){
      cout << endl << "Adjust Pivot row (" << icol << ") --> pivot=1" << endl;
      cout << "A = " << endl;  cpMatrixPrint( A );
      cout << "I = " << endl;  cpMatrixPrint( I );
    } // debug

    //---------------------
    // Reduce the other rows (expect for pivot row)
    //   such that pivot column elements are zero
    //---------------------
    for( int ll = 0; ll < n; ++ll ){
      if( ll != icol ){       // skip pivot row
	prv = A[ll][icol];    // save prev. value of elem to become zero
	for( int l = 0; l < n; ++l ){
	  A[ll][l] -= A[icol][l] * prv;
	  I[ll][l] -= I[icol][l] * prv;
	} // l loop
      } // ll != icol
    } // ll loop

    if( debug ){
      cout << endl << "Set pivot column (" << icol << ") elements to 0" << endl;
      cout << "A = " << endl;  cpMatrixPrint( A );
      cout << "I = " << endl;  cpMatrixPrint( I );
    } // debug

  } // i = main loop over columns

  if( debug ) cout << separator << endl << endl;

  //---------------------------
  // Copy adjusted I into A for return
  //---------------------------
  for( int i = 0; i < n; ++i ){
    for( int j = 0; j < n; ++j ){
      A[i][j] = I[i][j];
    } // j
  } // i

  return 0;
} // cpMatrixInvGJlong

//-------------------------------------------
// cpCholesky
//-------------------------------------------
void cpCholesky( cpMatrix& A, cpMatrix& L, int n )
{
  for (int r = 0; r < n; r++){
    for (int c = 0; c <= r; c++){

      // diagonal elements
      if (c == r){
	double sum = 0;
	for (int j = 0; j < c; j++)
	  sum += L[c][j] * L[c][j];
	L[c][c] = sqrt( A[c][c] - sum );
      } // c == r

      // off-diagonal elements
      else{
	double sum = 0;
	for (int j = 0; j < c; j++)
	  sum += L[r][j] * L[c][j];
	L[r][c] = ( A[r][c] - sum ) / L[c][c];
      } // c != r
      
    } // c-loop
  } // r-loop
} // cpCholesky

//------------------------------------------------------------------
// cpTriDiag
//------------------------------------------------------------------
int cpTriDiag( vector<double>& a, vector<double>& b, vector<double>& c,
	       vector<double>& r, vector<double>& u )
{
  int N = a.size();                // dimension of matrix
  double bet;
  vector<double> gam(N);

  if( b[0] == 0.0 ) return -1;     // need to eliminate one equation

  //-----------------------
  // Decomposition and Forward Substitution
  //-----------------------
  u[0] = r[0] / ( bet=b[0] );
  for( int j = 1; j < N; ++j ){
    gam[j] = c[j-1] / bet;
    bet    = b[j] - a[j]*gam[j];
    if( bet == 0.0 ) return -1;    // almost never happens in practice
    u[j]   = ( r[j] - a[j]*u[j-1] ) / bet;
  } // j loop (forward subst)

  //------------------------
  // Backward Substitution
  //------------------------
  for( int j = (N-2); j >= 0; --j )
    u[j] -= gam[j+1]*u[j+1];

  return 0;
} // cpTriDiag

//------------------------------------------------------------------
// cpSymSolve
//------------------------------------------------------------------
int cpSymSolve( cpMatrix& A, vector<double>& b, vector<double>& x )
{
  int N = b.size();
  
  //-------------------
  // Cholesky decomp: A = L L^T
  //-------------------
  cpMatrix L;
  cpMatrixInit( L,N,N,0 );
  cpCholesky( A,L,N );

  //------------------
  // Forward substitution to find solution to L.y = b
  //------------------
  vector<double> y( N,0 );
  y[0] = b[0] / A[0][0];
  for( int i = 1; i < N; ++i ){
    y[i] = b[i];
    for( int j = 0; j < i; ++j )
      y[i] -= A[i][j] * y[j];
    y[i] = y[i] / A[i][i];
  } // i loop

  //------------------
  // Back substitution to find solution to L^T.x = y
  //------------------
  x[N-1] = y[N-1] / A[N-1][N-1];
  for( int i = N-2; i >= 0; --i ){
    x[i] = y[i];
    for( int j = i+1; j < N; ++j )
      x[i] -= A[j][i] * x[j];
    x[i] = x[i] / A[i][i];
  } // i loop
  
  return 0;
} // cpSymSolve
