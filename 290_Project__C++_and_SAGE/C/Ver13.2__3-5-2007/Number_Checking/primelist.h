

 
// Define the Kronecker symbol (a/p)
//  (assuming p is a prime >= 2)
int KroneckerSymbol(long a, long p);


// =======================================================================================================



// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoExact(long p, long N, int charvalue);



// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoApprox(long p, long N, int charvalue);


// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoUniversal(long p, long N, int charvalue); 



// =======================================================================================================



// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBound(vector<long> & prime_list, vector<double> & F4_value_list, double bound, long N, long chi_top, vector<long> Big_Prime_List);




// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBoundUniversal(vector<long> & prime_list, vector<double> & F4_value_list, double bound, vector<long> Big_Prime_List);




// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBoundOptimized(vector<long> & prime_list, vector<double> & F4_value_list, double bound, long N, long chi_top, vector<long> Big_Prime_List);


// ---------------------------------------------------

#if !defined(PRIMELIST_H)
#define PRIMELIST_H

#include "primelist.cc"

#endif
