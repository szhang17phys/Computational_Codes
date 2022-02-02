#ifndef CPGENSAMPLE_H
#define CPGENSAMPLE_H

//******************************************************************************************
// cpGenSample.hpp: define classes & functions for generating TMVA training/test samples
//
// Classes
//
// Functions
//
// Modifications
// 30-Aug-19  move function Cholesky to cpMatrix (rename as cpCholesky)
//******************************************************************************************

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "cpUtils.hpp"
#include "cpRandom.hpp"
#include "cpMatrix.hpp"

//-----------------------
// typedefs as in NR 1.4.1
//-----------------------
typedef unsigned long long int Ullong;

//---------------------------------------------------
// PDF IDs for variable generation
//----------------------------------------------------
const int IDGAU = 0;  // Gaussian
const int IDUNI = 1;  // Uniform
const int IDEXP = 2;  // Exponential

const int MAXVARS = 10;  // max number of variables to generate


//==========================================================================================
// cpGSVarPars class
//==========================================================================================
class cpGSVarPars {
  // output method
  friend std::ostream& operator <<( std::ostream&, const cpGSVarPars& );

public:
  // Constructor: reads from parameter file
  cpGSVarPars( const int id, const std::string name, const std::vector<double>& pars )
    : _id(id), _name(name), _pars(pars) {}

  // getter functions
  int ID() { return _id; }
  std::string getName() { return _name; }
  double getPar( int i ) { return _pars[i]; }

private:
  int _id;              // identifier: Gaussian, Uniform, Exponential
  std::string _name;         // name for variable, e.g. gau0
  std::vector<double> _pars; // parameters
}; // cpGSVarPars


//==========================================================================================
// cpGenSample class
//==========================================================================================
class cpGenSample {
  // output method
  friend std::ostream& operator <<( std::ostream&, const cpGenSample& );

public:
  // Constructor: reads from parameter file
  cpGenSample( std::string inFile, std::string outFile ) : _inFile(inFile), _outFile(outFile) {}

  // Read in sample information
  int SampleInit();

  // Generate an event
  int GenEvent();

  // Stat's of generation
  Ullong GetSeed(){ return _seed; }
  int NumVars(){ return _nVar; }
  int NumEvts(){ return _nEvt; }
  double GenAve( int i ){ return( _sumX[i]/_nEvt ); }   // average for variable i
  // covariance matrix element i,j - not a very stable algorithm
  double GenCov( int i, int j ){ return( _sumXY[i][j]/_nEvt - GenAve(i)*GenAve(j) ); }

  // Output to file
  void writeTMVAHeader();
  void writeTMVAEvent();

  // Close output file
  void Close(){ _ofl.close(); }

private:
  std::string _inFile;         // input file name
  std::string _outFile;        // output file name
  std::ofstream _ofl;          // output file
  int _nEvt;              // number of events generated so far
  int _nVar;              // number of variables to generate per event
  int _nGau;              // number of Gaussian variables in variable list
  std::vector<int> _indGau;    // indices of Gaussian variables in variable list

  Ullong _seed;           // seed for random number generator
  cpRandom* _Gen;         // random number generator

  cpGSVarPars* _varPars[MAXVARS]; // container for variables
  double _var[MAXVARS];           // generated variables
  double _sumX[MAXVARS];          // running sum of variables generated (for average)
  cpMatrix _sumXY;                // for covariance of generated var's

  int _nCor;             // number of correlations read in
  std::vector<int> _ind1;     // correlation indices and values read in
  std::vector<int> _ind2;
  std::vector<double> _corr;
  cpMatrix _mCor;         // correlation matrix: nVar x nVar
  cpMatrix _mCov;         // covariance matrix: nGau x nGau
  cpMatrix _mCho;         // Cholesky decomposition of mCov: nGau x nGAu
}; // cpGenSample

#endif
