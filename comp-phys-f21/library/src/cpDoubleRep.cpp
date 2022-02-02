#include "cpDoubleRep.hpp"
using namespace std;

//===========================================================
// PrintRep: pretty printout of double bit representation
//===========================================================
string PrintRep( int nExp, int nMan, bitset<MAXBITS>& bits )
{
  int nTot = 1 + nExp + nMan;   // sign bit + exp bits + mantissa bits

  ostringstream os;
  os << bits[nTot-1] << "  ";                                  // sign bit

  for( int i = nTot-2; i > nTot-2-nExp; i-- ) os << bits[i];   // exponent

  os << "  ";
  for( int i = nTot-2-nExp; i >= 0; i-- ) os << bits[i];       // mantissa

  return( os.str() );
} // PrintRep


//=========================================================
// BitSetULL: allows construction of bitset from unsigned long long
//   standard bitset only allows initialization with a 32-bit int
//=========================================================
bitset<64> BitSetULL( const Ullong ull )
{
  bitset<32> lo32( ull );      // low-order bits
  bitset<32> hi32( ull>>32 );  // high-order bits
  bitset<64> b64;
  for( int i = 0; i < 32; ++i ){
    b64[i+32] = hi32[i];
    b64[i]    = lo32[i];
  }
  return b64;
} // BitSetULL

