#ifndef CPMATRIX_H
#define CPMATRIX_H

//**********************************************************************
// cpMatrix: a few useful matrix manipulation functions
//
// Typedefs
//   cpMatrix   currently just vector<vector<double> >
//
// Functions
//   cpMatrixInit   initializes an nRowxnCol matrix to <val=0>
//   cpMatrixInvGJ  matrix inversion by Gauss-Jordan eleimination
//   cpTriDiag      solves system of eqn's governed by tridiagonal matrix
//
// Modifications
// 30-Aug-19  add Cholesky (moved from cpGenSample)
//            add cpSymSolve
// 25-Jan-11  add cpTriDiag
//**********************************************************************

#include <cmath>
#include <vector>
#include <string>


//-------------------------------------------------------------------
// cpMatrix typedef
// This should really be a class w/ all matrix operations defined
//   but we'll do things simply to start with
//-------------------------------------------------------------------
typedef std::vector<std::vector<double> > cpMatrix;

//------------------------------------------------------------------
// cpMatrixInit
//------------------------------------------------------------------
void cpMatrixInit( cpMatrix& A, int nRow, int nCol, double val=0 );

//-------------------------------------------------------------------
// cpMatrixPrint: prints out matrix to cout
//-------------------------------------------------------------------
void cpMatrixPrint( const cpMatrix& A );

//-------------------------------------------------------------------
// cpMatrixInvGJ: Gauss-Jordan matrix inversion using Full Pivoting
//    simplified version of NR function gaussj
// Inputs:
//   A     matrix to be inverted
//   debug prints out intermediate steps
// Outputs:
//   A     input matrix is changed on output to be A^{-1}
// Returns: 0=success, 1=singular
//------------------------------------------------------------------
int cpMatrixInvGJ( cpMatrix& A );
int cpMatrixInvGJlong( cpMatrix& A, const bool debug=true );  // debug version

//------------------------------------------------------------------
// Cholesky: does Cholesky decomposition on a symmetric matrix
//   A = covariance matrix (n x n)
//   L = Cholesky decomposition (n x n, lower triangular)
//   A = L L^T
// Uses algorithm from: https://rosettacode.org/wiki/Cholesky_decomposition#C.23
//   l_{kk} = sqrt( a_{kk} - sum_{j=1}^{k-1} l^2_{kj} )
//   l_{ik} = ( a_{ik} - sum_{j=1}^{k-1} l_{ij} l_{kj} ) / l_{kk}
//
// Inputs:
//   A  symmetric matrix
//   n  dimension of matrix
// Outputs:
//   L  lower-diagonal matrix
//------------------------------------------------------------------
void cpCholesky( cpMatrix& A, cpMatrix& L, int n );

//------------------------------------------------------------------
// cpTriDiag: solves A.u = r for tridiagonal A
// (from NR function tridiag - sect 2.4)
//
// Inputs:
//   a,b,c  vectors containing diag-1, diag, diag+1 components of matrix
//          note: a[0] and c[N-1] are unused
//   r      rhs of equation
// Outputs:
//   u      solution vector
// Returns: 0 = success
//------------------------------------------------------------------
int cpTriDiag( std::vector<double>& a, std::vector<double>& b, std::vector<double>& c,
	       std::vector<double>& r, std::vector<double>& u );

//------------------------------------------------------------------
// cpSymSolve: solves A.x = b for symmetric A
//     uses Cholesky decomposition: A = L L^T  (L = lower diagonal)
//
// Inputs:
//   A  symmetric matrix (NxN)
//   b  vector (N) - rhs of equation
// Outputs:
//   x  solution vector (N)
// Returns: 0 = success
//------------------------------------------------------------------
int cpSymSolve( cpMatrix& A, std::vector<double>& b, std::vector<double>& x );

#endif
