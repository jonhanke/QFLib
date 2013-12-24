
#if !defined(MISC_PRIMES_H)
#define MISC_PRIMES_H


////////////////////////////////////////////
// Global variables used to handle primes //
////////////////////////////////////////////

  // Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
  // Our largest prime is less than 15 million, so we only need the mpz_class precision 
  // for the square-free numbers. =)


// Declare the external global variables relating to the primes
extern const char PRIME_DIR[14];
extern const string PRIME_DIR_STRING;
extern const string PRIMEFILE_LIST[3];
extern const long NUMBER_OF_PRIMEFILES;
extern long PRIMEFILE_INDEX;
extern vector<long> Big_Prime_List;




void Use_Next_Primefile();

#include "misc_primes.cc"

#endif
