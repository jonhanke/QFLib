
#include<gmpxx.h>


/*  <<====  THIS IS ALREADY DEFINED IN THE EARLIER LOCAL DENSITY PROJECT...

// LCM of two mpz_classes
mpz_class LCM(mpz_class a, mpz_class b); 

*/




// Print the current date and time 
void PrintTime();

string Time();              // <-------- Need to fix this later!!!


// -------------------------------------------------

// Finds the level of a ternary form Q
long QF_Ternary_Level(long Q[6]);


/////////////////////////////////////////////////////////////
// Computes the valuation of an integer m mod some prime p //
/////////////////////////////////////////////////////////////
// Note: Assumes the given p is a prime.
 

//unsigned long Valuation(const long & m, const long & p);            // This is Depricated....
unsigned long Valuation(const mpz_class & m, const long & p);


// Takes a power of a long
long LongPow(long a, unsigned long b);

// Finds the prime divisors of m (in order)
//vector<long> PrimeDivisors(long m);                                   // This is Depricated...
vector<long> PrimeDivisors(mpz_class m);                                  


// Make all subsets of a given set S, containing exactly num elements.
template <class T>
vector< set<T> > MakeSubsets(const set<T> & S, const long & num);




// -----------------------------------------------

#if !defined(MISC_UTILITIES_H)
#define MISC_UTILITIES_H

#include "misc_utilities.cc"

#endif
