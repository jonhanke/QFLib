

// Check the representability of a vector of eligible numbers.
vector<mpz_class> CheckNumbers(const vector<mpz_class> & num, int depth, const vector<bool> & theta_vec, long diag_coeff);





// Check the representability of a vector of eligible numbers.
void NewCheckNumbers(const vector<mpz_class> & num, int depth, const vector<bool> & theta_vec, long diag_coeff, 
		vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
		vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths);




// Check the representability of a vector of eligible numbers.
void NewerCheckNumbers(const vector<mpz_class> & num, int depth, const unsigned long theta_vec[], 
		       unsigned long Theta_Precision, long diag_coeff, 
		       vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
		       vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths);




// ----------------------------------

#if !defined(CHECKNUMS_H)
#define CHECKNUMS_H
 
#include "checknums.cc"

#endif
