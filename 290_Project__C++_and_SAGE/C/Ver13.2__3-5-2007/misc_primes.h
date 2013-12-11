
#if !defined(MISC_PRIMES_H)
#define MISC_PRIMES_H


////////////////////////////////////////////
// Global variables used to handle primes //
////////////////////////////////////////////

// This gives the default directory for where the primefiles are stored,
// which is overridden by the shell variable QFLIB_PRIME_DIR.
//extern const char PRIME_DIR[] = "/u/jonhanke/QFLIB_using_hg/qflibrary__2012-06-04__hg_repo_on_Google_code/290_Project__C++_and_SAGE/Primes/";
extern const char PRIME_DIR[] = "../../Primes/";

// Define a list of space-delimited text files containing the primes up to some bound 
// in increasing order, with the list of files in order of increasing size.
const string PRIMEFILE_LIST[2] =
  {
    "primes_lt_1000000.txt",
    "primelist_upto_15millionth_prime.txt"
  };
const long NUMBER_OF_PRIMEFILES = 2;

// Define a global index pointing to the current primefile in use from primefile_filelist
long PRIMEFILE_INDEX = -1;






  // Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
  // Our largest prime is less than 15 million, so we only need the mpz_class precision 
  // for the square-free numbers. =)


// The global list of primes used throughout.
vector<long> Big_Prime_List;










void Use_Next_Primefile();

#include "misc_primes.cc"

#endif
