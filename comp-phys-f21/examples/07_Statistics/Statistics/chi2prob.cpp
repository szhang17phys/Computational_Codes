//********************************************************************
// chi2prob.cpp: returns chi^2 probability or inverse chi^2 probability
//  chi^2 prob     = P(chi^2; Ndof)
//  inv chi^2 prob ==> returns x such that P(x; Ndof) = p
//
// Modifications
// 25-Sep-18  merge chi2prob and inv_chi2prob
// 27-Oct-10  created
//********************************************************************
#include <limits>

#include "cpUtils.hpp"
#include "cpStatistics.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  if( argc < 3 ){
    cout << "chi2prob usage: ./chi2prob <chi^2 or prob> <N_dof> [inv]" << endl;
    return -1;
  }

  double arg = atod( argv[1] );  // either chi^2 (default) or prob (if 4 argument present) 
  double Ndof = atod( argv[2] );

  // inverse chi^ probability
  if( argc == 4 ){
    cout << "Chi^2 Prob = " << arg
	 << " (" << Ndof << " dof)"
	 << "  ==>  x = " << cpInvProbChi2( arg,Ndof ) << endl;
  } // inv chi^2 prob

  // chi^2 probability
  else{
    cout << "Chi^2/dof = " << arg/Ndof
	 << "  (" << arg << " / " << Ndof << ")"
	 << "  ==>  Prob = " << cpProbChi2( arg,Ndof ) << endl;
  } // chi^2 prob

  return 0;
} // main
