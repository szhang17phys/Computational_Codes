#ifndef CPRANDOM_H
#define CPRANDOM_H

//************************************************************
// cpRandom: Uniform (and non) Deviates random number generators
//   used as base functions for generation
//
// Random numbers (uniform and non), optionally over lo-hi
// (these are available to all derived classes)
// o Ullong randInt64( [lo],[hi] )   : 64-bit default [0 - 2^64-1]
// o Uint   randInt32( [lo],[hi] )   : 32-bit default [0 - 2^32-1]
// o double randDouble( [lo],[hi] )  : double default [0.0 - 1.0]
// o double randDevNrm( [mean],[sigma] ) : Gaussian: mean[=0], sigma[=1]
// o int randDevPoisson( [lambda] ) : Poisson distrib
//
// Modifications
// 12-Nov-18  fix bug in randDevExp (1/lambda instead of lambda)
// 10-Jul-11  add cpRandMT
//            add idMT
// 02-Apr-10  change randDevNrm to use Ratio-of-Uniforms
//            previous version now called randDevNrmBM
// 30-Mar-10  add randDevExp
// 07-Mar-10  add writeParams(), sPrng, change PRNG id's
//************************************************************

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>

//-----------------------
// typedefs as in NR 1.4.1
//-----------------------
typedef unsigned int           Uint;
typedef unsigned long long int Ullong;

//-----------------------
// Generator IDs (makes testing easier)
//-----------------------
const int NPRNGS   = 5;
const int idD1A1r  = 0;
const int idLocal  = 1;
const int idLCG    = 2;
const int idDevNrm = 3;
const int idMT     = 4;

const std::string sPrng[NPRNGS] = { "D1(A1_r)","rand()","LCG","NormDev","MT19937-64" };

//============================================================
// cpRandom: base class - all others inherit from this
//  contains Numerical Recipes suggested generator (Ranq1)
//  Ranq1 = D1(A1_r) (sect 7.1.3) - period ~ 1.8x10^19
//
// Constructor Inputs
//   Ullong seed : random number seed
//============================================================
class cpRandom {
public:
  // constructors
  cpRandom( Ullong seed ) : _seed(seed), _v( 4101842887655102017LL ) {
    _v ^= seed;
    _v  = randInt64();
    _y2 = 0;
  }

  // print out PRNG info
  virtual inline void writeParams(){
    std::cout << "D1(A1_r): a = 4101842887655102017 ; a1 = 21, a2 = 35, a3 = 4 - seed = " 
	 << getSeed() << std::endl;
  }

  // access to seed
  Ullong getSeed() { return _seed; }

  //--------------------
  // The Basic Generators: all integers or doubles from 0.0-1.0
  //--------------------
  virtual inline Ullong randInt64(){                 // 64-bit Generator
    _v ^= _v >> 21; _v ^= _v << 35; _v ^= _v >> 4; // 64-bit Xor-shift
    return _v * 2685821657736338717LL; }           // MLCG mod 2^64

  virtual inline Uint randInt32() {                  // 32-bit Generator
    return (Uint)randInt64(); }

  virtual inline double randDouble(){                // Double Generator: [0-1]
    return 5.42101086242752217e-20 * randInt64(); } // 1/2^64

  //--------------------
  // Generate a number in an interval: lo <= x <= hi
  //--------------------
  inline Ullong randInt64( Ullong lo, Ullong hi ){
    return( randInt64()%(hi+1-lo) + lo ); }

  inline Uint randInt32( Uint lo, Uint hi ){
    return( (Uint)randInt64( lo,hi ) ); }

  inline double randDouble( double lo, double hi ){
    return( randDouble()*(hi-lo) + lo ); }

  //--------------------------
  // Non-uniform deviates
  //--------------------------
  double randDevNrm();                          // Gaussian( 1,0 )
  double randDevNrm( double mu, double sig ){   // Ratio-of-Uniforms
    return( mu + sig*randDevNrm() ); }

  double randDevNrmBM();                          // Gaussian( 1,0 )
  double randDevNrmBM( double mu, double sig ){   // Box-Muller
    return( mu + sig*randDevNrmBM() ); }

  double randDevExp( double lambda );           // lambda exp( -lambda x )

  int randDevPoisson( double lambda );          // Poisson

protected:
  Ullong _seed;  // the seed used
  Ullong _v;     // previous value (x_i)
  double _y2;    // for randDevNrm

}; // cpRandom


//============================================================
// cpRandLocal: implements local random number generator (rand)
//   for testing purposes ONLY
//   inherits from cpRandom
// Note: DO NOT use 64-bit generator (randInt64) 
// - it only does 32-bit numbers
//============================================================
class cpRandLocal : public cpRandom {
public:
  // constructors
  cpRandLocal( Ullong seed ) : cpRandom( seed ) {  
    srand(seed); }  // set seed

  // print out PRNG info
  inline void writeParams(){
    std::cout << "Local PRNG: no params - seed = " << getSeed() << std::endl; }

  // The Generators
  inline Uint randInt32() {                  // 32-bit generator
    return (Uint)rand(); }  

  inline double randDouble(){                // [0.0 - 1.0]
    return (double)rand() / (double)RAND_MAX; }

  inline Ullong randInt64(){                 // NON-FUNCTIONAL 64-bit Generator
    return (Ullong)rand(); }

}; // cpRandLocal


//============================================================
// cpRandLCG: Linear Congruential Generator
//   x_n = ( a x_(n-1) + c ) mod m
//   inherits from cpRandom
//
// Constructor Inputs (in addition to those for cpRandom)
//   Ullong a,c,m : parameters of LCG as above
//============================================================
class cpRandLCG : public cpRandom {
public:
  // constructors
  cpRandLCG( Ullong seed, Ullong a, Ullong c, Ullong m ) : 
    cpRandom(seed), _v(seed), _a(a), _c(c), _m(m) {}

  // print out PRNG info
  inline void writeParams(){
    std::cout << "LCG: a = " << _a << ", c = " << _c << ", m = " << _m 
	 << " - seed = " << getSeed() << std::endl; }

  // The Generators
  inline Ullong randInt64(){         // 64-bit Generator
    return ( _v = ( _a*_v + _c ) % _m ); }

  inline Uint randInt32() {          // 32-bit generator
    return (Uint)randInt64(); }  

  inline double randDouble(){        // [0.0 - 1.0]
    return( (double)randInt64() / (double)_m ); }

protected:
  Ullong _v;          // previous value
  Ullong _a,_c,_m;    // LCG parameters

}; // cpRandLCG


//============================================================
// cpRandMT: Mersenne Twister
//   inherits from cpRandom
//
// Code based on: ROOT::TRandom3
//   64-bit version from Takuji Nishimura and Makoto Matsumoto
//   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
//============================================================
class cpRandMT : public cpRandom {
public:
  // constructors
  cpRandMT( Ullong seed ) : cpRandom(seed) { init(seed); }

  // print out PRNG info
  inline void writeParams(){
    std::cout << "Mersenne Twister: seed = " << getSeed() << std::endl; }

  // initialization
  void init( Ullong seed );

  // The Generators
  Ullong randInt64();                // 64-bit Generator

  inline Uint randInt32() {          // 32-bit generator
    return (Uint)randInt64(); }  

  inline double randDouble(){        // Double Generator: [0-1]
    return 5.42101086242752217e-20 * randInt64(); } // 1/2^64

protected:
  Ullong _fMt[312];   // state vector array
  int _fCount;        // counter

}; // cpRandMT



#endif
