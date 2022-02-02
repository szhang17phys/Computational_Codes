#ifndef CPDOUBLEREP_H
#define CPDOUBLEREP_H

//**********************************************************
// cpDoubleRep: extract computer representation of a double
//   this could probably be done more generally with templates
//   but let's keep things simple to start with
//
// Modifications
// 16-Aug-10  add float conversion (Ufloat) + changes to PrintRep
// 06-Jan-10  cast on pow to be consistent w/ Mac OS
// 17-Nov-09  modularized
// 25-Sep-09  created
//**********************************************************

#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

typedef unsigned long long int Ullong;

//=========================================================
// Constants
//=========================================================
const int MAXBITS = 64; // maximum number of bits allowed

// Single Precision - 32 bits
const int F_SGNBITS = 1;  // number of bits for sign
const int F_EXPBITS = 8;  // exponent
const int F_MANBITS = 23; // mantissa
const int F_EOFF =
    (int)pow(2.0, F_EXPBITS - 1) - 1; // exp offset (127 for float)

// Double Precision - 64 bits
const int D_SGNBITS = 1;  // number of bits for sign
const int D_EXPBITS = 11; // exponent
const int D_MANBITS = 52; // mantissa
const int D_EOFF =
    (int)pow(2.0, D_EXPBITS - 1) - 1; // exp offset (1023 for double)

//=========================================================
// Udoub/Ufloat: gets the computer representation of a double/float directly
//   using a union
//=========================================================
union Udoub {
  double d; // double
  Ullong i; // 64-bit unsigned int
};          // Udoub

union Ufloat {
  float d;        // float
  unsigned int i; // 32-bit unsigned int
};                // Udoub

//=========================================================
// Functions
//=========================================================

// PrintRep: pretty printout of bit representation of a double
//   nExp = number of bits in exponent
//   nMan = number of bits in mantissa
std::string PrintRep(int nExp, int nMan, std::bitset<MAXBITS> &bits);

// BitSetULL: allows construction of bitset from a 64-bit int
//   standard bitset only allows initialization with a 32-bit int
std::bitset<64> BitSetULL(const Ullong ull);

#endif
