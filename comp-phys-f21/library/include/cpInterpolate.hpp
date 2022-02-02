#ifndef CPINTERPOLATE_H
#define CPINTERPOLATE_H

//************************************************************************
// cpInterpolate: interpolation functions
//
// Modifications
// 14-Oct-09  created
//************************************************************************

#include <vector>

//-------------------------------
// PolyInterp: O(n-1) polynomial interpolation using n data points
//   xv  = vector of x values
//   yv  = vector of y values corresponding to xa
//   ist = starting index in above vectors
//         will use points (x,y)[ist]...(x,y)[ist+n-1]
//   n   = number of points to consider
//         interpolating polynomial will be O(n-1)
//   x   = point at which to evaluate interpolated function
//   dy  = (OUTPUT) error est on y from O(n-1) - O(n-2) est's of y
// returns: y (value of interpolated function at x)
//-------------------------------
double PolyInterp( const std::vector<double>& xv, const std::vector<double>& yv,
		   const int ist, const int n, const double x, double& dy );

#endif
