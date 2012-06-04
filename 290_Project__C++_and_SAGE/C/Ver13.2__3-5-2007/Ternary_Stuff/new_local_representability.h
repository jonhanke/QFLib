

////////////////////////////////////////////////
// Get the local conditions on a ternary form //
////////////////////////////////////////////////
void FindTernaryLocalConditions_New(long Q[6], vector<long> & local_mod_vector,
				    long local_repn_array[][9], vector<long> & aniso_vector); 




//////////////////////////////////////////////////
// Get the local conditions on a quadratic form //
//////////////////////////////////////////////////
void FindLocalConditions(Matrix_mpz QQ, vector<long> & local_mod_vector,
				    long local_repn_array[][9], vector<long> & aniso_vector);



// ---------------------------------------------

#if !defined(NEW_LOCAL_REPRESENTABILITY_H)
#define NEW_LOCAL_REPRESENTABILITY_H

#include "new_local_representability.cc"

#endif
