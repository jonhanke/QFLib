






//! Check the representability of a vector of eligible numbers.
void NewestCheckNumbers(const vector<mpz_class> & num, int depth, const unsigned long theta_vec[], 
			unsigned long Theta_Precision, long diag_coeff, 
			vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
			vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths,
			vector<long> & local_mod_vector, long local_repn_array[200][9]);




//! An incremental version which checks as it goes along... =) 

/*! 
  This is the detailed comment... but I probably need to say something more!
  \f$\sqrt{ (x_2-x_1)^2 + (y_2 - y_1)^2 }\f$
  Let's try using the formula \f$\int_a^b f(x) dx \f$. =)
*/

vector<mpz_class> NewestSharpList(vector<mpz_class> & T_list, double bound, //< Bound for F_4(p).
		     const vector<long> & PP, const vector<double> & F4PP, 
		     unsigned long Ternary_series[], unsigned long Precision, 
				  long diag_coeff, Matrix_mpz QQ, long Ternary[6]);


// -------------------------

#if !defined(CHECKING_H)
#define CHECKING_H



#include "checking.cc"

#endif
