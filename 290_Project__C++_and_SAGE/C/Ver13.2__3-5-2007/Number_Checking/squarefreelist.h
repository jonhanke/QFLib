

//////////////////////////////////////////////////////////////////////
// Original routine to make the list of eligible squarefree numbers //
//////////////////////////////////////////////////////////////////////

void SharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP);




////////////////////////////////////////////////////////////////
// Forgetful version which doesn't keep the list in memory =) //
////////////////////////////////////////////////////////////////

void ForgetfulSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP);


///////////////////////////////////////////////////////////////
// Lazy version which doesn't write the list to the disk. =) //
///////////////////////////////////////////////////////////////

void LazySharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP);





////////////////////////////////////////////////////////////////
// An incremental version which checks as it goes along... =) //
////////////////////////////////////////////////////////////////

void NewSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP, 
		  vector<bool> & Ternary_series, long diag_coeff);





////////////////////////////////////////////////////////////////
// An incremental version which checks as it goes along... =) //
////////////////////////////////////////////////////////////////

void NewerSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP, 
		  unsigned long Ternary_series[], unsigned long Precision, long diag_coeff);






// =========================================================================================================





///////////////////////////////////////////////////////////////////////////////////////
// Modified routine to check the bound overflow after multiplying each new factor... //
///////////////////////////////////////////////////////////////////////////////////////

void FasterSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP);





// ---------------------------------------------------

#if !defined(SQUAREFREELIST_H)
#define SQUAREFREELIST_H

#include "squarefreelist.cc"

#endif
