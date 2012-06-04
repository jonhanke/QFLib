



/////////////////////////////////////////////////////////
// This finds the ternary sublattice contained in the  //
// orthogonal complement of the upper left entry.      //
/////////////////////////////////////////////////////////

Matrix_mpz FindDecomposableTernarySublatticeComplementof1(Matrix_mpz Q); 


//////////////////////////////////////////////////////////////////////
// Writes the theta function of a ternary form of desired precision //
//////////////////////////////////////////////////////////////////////

//! \deprecated  This has been moved to boolean_ternary_theta::compute().

void MakeTernaryTheta(long QQ[6], unsigned long Ternary_Precision);


//////////////////////////////////////////////////////////////////////
// Checks the representability of the square free numbers of a form //
//////////////////////////////////////////////////////////////////////

void Check4VarRepresentability(long Ternary[6], double B, long N, long char_top, 
			       unsigned long Precision, vector<long> & all_primes, 
			       Matrix_mpz Q, long form_num);






//////////////////////////////////////////////////////////////////
// Identifies the upper-left ternary in each of our 6560 forms, //
// and describes those whose ternary is not regular.            //
//////////////////////////////////////////////////////////////////

void CheckAll6560Forms(vector<long> & all_primes);



// ---------------------------------------------------

#if !defined(SCRAPS4_H)
#define SCRAPS4_H

#include "readmagmaseries.cc"  // This is needed for now because of line 217 in scraps.cc...
#include "maketheta_small.cc"  // This is needed for the ThetaPrecision routine.

// These are the testing routines...
// (they should be moved later...)
#include "local_condition_class.h"
#include "more_local.cc"
#include "more_local__tests.cc"

// This is the real deal
#include "scraps4.cc"

#endif
