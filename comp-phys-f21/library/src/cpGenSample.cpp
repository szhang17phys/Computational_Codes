#include <stdlib.h>
#include <fstream>
#include <cmath>

#include "cpGenSample.hpp"
using namespace std;

//====================================================================
// cpGSVarPars classes
//====================================================================

//---------------------------------------
// << operator
//---------------------------------------
ostream& operator <<( ostream& os, const cpGSVarPars& cls ){
  os << cls._name;
  for( int i = 0; i < cls._pars.size(); ++i )
    os << " " << cls._pars[i];

  return os;
}; // <<


//====================================================================
// cpGenSample classes
//====================================================================

//---------------------------------------
// << operator
//---------------------------------------
ostream& operator <<( ostream& os, const cpGenSample& cls ){
  os << "Seed: " << cls._seed << endl << endl;
  
  for( int i = 0; i < cls._nVar; ++i )
    os << "Variable " << i << " : " << " " << *(cls._varPars[i]) << endl;
  os << "Gaussian Variables: " << cls._nGau << endl;

  if( cls._nCor > 0 ){
    os << endl << "Correlation Matrix" << endl;
    cpMatrixPrint( cls._mCor );
  } // _nCor > 0

  return os;
}; // <<

//-------------------------------------------
// SampleInit: reads parameters of the sample to generate from a file
//-------------------------------------------
int cpGenSample::SampleInit()
{
  //-------------------------------------
  // Initialize counters
  //-------------------------------------
  _nEvt = 0;
  _nVar = 0;
  _nGau = 0;
  _nCor = 0;
  _seed = -1;  // default RNG seed

  //-------------------------------------
  // Open Input File
  //-------------------------------------
  ifstream ifl( _inFile.data() );
  if( !ifl ){
    cout << "cpGenSample: could not find file: " << _inFile << endl;
    exit( EXIT_FAILURE );
  } // !ifl

  //------------------------------------
  // Read data from file
  //------------------------------------
  string arg,line,name;
  while( ifl >> arg ){

    // PRNG seed
    if( arg == "SEED" ){
      ifl >> _seed;
      getline( ifl,line );
      
      if( _seed == -1 ) _seed = (Ullong)clock();
    } // arg == SEED
    
    // GAU = Gaussian variable
    else if( arg == "GAU" ){
      double mean,sigma;
      vector<double> pars;
      ifl >> name >> mean >> sigma;
      pars.push_back( mean );
      pars.push_back( sigma );
      _varPars[_nVar] = new cpGSVarPars( IDGAU,name,pars );
      _indGau.push_back( _nVar );  // index in variable list for this Gaussian variable
      _nVar++;
      _nGau++;
      getline( ifl,line );
    } // arg == GAU
    
    // UNI = Uniform variable
    else if( arg == "UNI" ){
      double hi,lo;
      vector<double> pars;
      ifl >> name >> lo >> hi;
      pars.push_back( lo );
      pars.push_back( hi );
      _varPars[_nVar] = new cpGSVarPars( IDUNI,name,pars );
      _nVar++;
      getline( ifl,line );
    } // arg == UNI
    
    // EXP = Exponential variable
    else if( arg == "EXP" ){
      double lambda;
      vector<double> pars;
      ifl >> name >> lambda;
      pars.push_back( lambda );
      _varPars[_nVar] = new cpGSVarPars( IDEXP,name,pars );
      _nVar++;
      getline( ifl,line );
    } // arg == UNI
    
    // CORR = correlation coefficient: ind1 ind2 correlation
    // ind1,2 = indices of the two variables in the variable list
    // correlation = [-1,1]
    else if( arg == "CORR" ){
      int ind1,ind2;
      double corr;
      ifl >> ind1 >> ind2 >> corr;
      _ind1.push_back( ind1 );
      _ind2.push_back( ind2 );
      _corr.push_back( corr );
      _nCor++;
      getline( ifl,line );
    } // arg == UNI

    // anything not recognized is treated as a comment
    else{
      getline( ifl,line );
    } // else
  } // while

  //------------------------------------
  // Set up Covariance & Cholesky matrices
  //------------------------------------
  if( _nCor > 0 ){
    // set up correlation matrix (nVar x nVar)
    cpMatrixInit( _mCor, _nVar, _nVar, 0 );
    for( int i = 0; i < _nCor; ++i ){  // load correlations into nVar x nVar matrix
      _mCor[_ind1[i]][_ind2[i]] = _corr[i];
      _mCor[_ind2[i]][_ind1[i]] = _corr[i];
    } // i = _nCor loop

    // set up covariance matrix (nGau x nGau)
    cpMatrixInit( _mCov, _nGau, _nGau, 0 );
    for( int r = 0; r < _nGau; ++r ){ // load covariance matrix
      _mCov[r][r] = _varPars[_indGau[r]]->getPar(1) * _varPars[_indGau[r]]->getPar(1); // diagonal
      for( int c = 0; c < r; ++c ){
	_mCov[r][c] = _mCor[_indGau[r]][_indGau[c]]                                    // off-diag
	  * _varPars[_indGau[r]]->getPar(1) * _varPars[_indGau[c]]->getPar(1);
	_mCov[c][r] = _mCov[r][c];
      } // c = col loop
    } // r = row loop

    // set up Cholesky decomposition matrix L(nGau x nGau)
    // \vec{y} = \vec{mu} + L \vec{x}
    //    L L^T = Cov
    //    \vec{x} = uncorrelated Gaussian var's with mean=0, var=1
    cpMatrixInit( _mCho, _nGau, _nGau, 0 );
    cpCholesky( _mCov, _mCho, _nGau );

  } // _nCor > 0

  //-----------------------------------
  // Initialize generated variable statistics
  //-----------------------------------
  for( int i = 0; i < MAXVARS; ++i ) _sumX[i] = 0;
  cpMatrixInit( _sumXY, _nVar, _nVar, 0 );
  
  //------------------------------------
  // Initialize Random Number Generator
  //------------------------------------
  _Gen = new cpRandMT( _seed );

  //------------------------------------
  // Open output file
  //------------------------------------
  ifl.close();             // first close input file
  _ofl.open(_outFile);     // now open output file

  return 0;
} // SampleInit


//-------------------------------------------
// TMVA output to file
//-------------------------------------------
// TMVA header from User Guide Code Example 9
void cpGenSample::writeTMVAHeader()
{
  for( int i = 0; i < _nVar; ++i ){
    _ofl << _varPars[i]->getName() << "/F";
    if( i < _nVar-1 ) _ofl << ":";
  } // i = _nVars
  
  _ofl << endl;
} // writeTMVAHeader

// writes out variable values generated for this event in TMVA format
void cpGenSample::writeTMVAEvent()
{
  for( int i = 0; i < _nVar; ++i )
    _ofl << _var[i] << " ";
  _ofl << endl;
} // writeTMVAEvent



//-------------------------------------------
// GenEvent: generates variables for an event
//-------------------------------------------
int cpGenSample::GenEvent()
{
  //----------------------------------
  // loop over variables
  //----------------------------------
  for( int i = 0; i < _nVar; ++i ){
    switch( _varPars[i]->ID() ){
    // Uniform: [lo,hi]
    case IDUNI:
      _var[i] = _Gen->randDouble( _varPars[i]->getPar(0),_varPars[i]->getPar(1) );
      break;
      
    // Exponential: lambda
    case IDEXP:
      _var[i] = _Gen->randDevExp( _varPars[i]->getPar(0) );
      break;
      
    // Gaussian: mean=0, sigma=1
    default:
      _var[i] = _Gen->randDevNrm();
      break;      
    } // switch
  } // i = vars

  //----------------------------------
  // Adjust Gaussian variables for mean/sigma & correlations
  //----------------------------------
  vector<double> y( _nGau,0 );  // temp storage for correlated var's
  for( int i = 0; i < _nGau; ++i ){

    // Correlated Gaussian variables
    // y_i = mu_i + L_{ij} x_j
    if( _nCor > 0 ){
      y[i] = _varPars[_indGau[i]]->getPar(0);
      for( int j = 0 ; j < _nGau; ++j )
	y[i] += _mCho[i][j] * _var[_indGau[j]];
    } // _nCor > 0

    // uncorrelated Gaussian variables
    else{
      y[i] = _varPars[_indGau[i]]->getPar(0)
	+ _var[_indGau[i]] * _varPars[_indGau[i]]->getPar(1);
    } // _nCor = 0
    
  } // i-loop (_nGau)

  // replace generated/uncorrelated Gaussian var's with adjusted/correlated ones
  for( int i = 0; i < _nGau; ++i )
    _var[_indGau[i]] = y[i];

  // update running sums for generated stat's
  for( int i = 0; i < _nVar; ++i ){
    _sumX[i] += _var[i];
    for( int j = 0; j < _nVar; ++j )
      _sumXY[i][j] += _var[i] * _var[j];
  } // i = _nVar

  /*********************************
  cout << "Event " << _nEvt << " :";
  for( int i = 0; i < _nVar; ++i )
    cout << " " << _var[i];
  cout << endl;
  **********************************/

  // update event counter
  _nEvt++;
  return 0;
} // GenEvent
